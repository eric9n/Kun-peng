use crate::mmscanner::{scan_sequence, MinimizerIterator};
use crate::reader::Reader;
use crate::seq::{Base, SeqFormat};
use crate::{detect_file_format, FastaReader, FastqReader, Meros};
use crossbeam_channel::{bounded, Receiver};
use scoped_threadpool::Pool;
use std::collections::HashMap;
use std::io::Result;
use std::sync::Arc;

pub struct ParallelResult<P>
where
    P: Send,
{
    recv: Receiver<P>,
}

impl<P> ParallelResult<P>
where
    P: Send,
{
    #[inline]
    pub fn next(&mut self) -> Option<P> {
        self.recv.recv().ok()
    }
}

pub fn create_reader(
    file_pair: &[String],
    file_index: usize,
    score: i32,
) -> Result<Box<dyn Reader + Send>> {
    // let mut files_iter = file_pair.iter();
    let paths = crate::OptionPair::from_slice(file_pair);

    match detect_file_format(&file_pair[0])? {
        SeqFormat::Fastq => Ok(Box::new(FastqReader::from_path(paths, file_index, score)?)),
        SeqFormat::Fasta => Ok(Box::new(FastaReader::from_path(&file_pair[0], file_index)?)),
    }
}

pub fn read_parallel<R, W, O, F, Out>(
    reader: &mut R,
    n_threads: usize,
    meros: &Meros,
    work: W,
    func: F,
) -> Result<()>
where
    R: Reader,
    O: Send,
    Out: Send + Default,
    W: Send + Sync + Fn(&mut Vec<Base<MinimizerIterator>>) -> Option<O>,
    F: FnOnce(&mut ParallelResult<Option<O>>) -> Out + Send,
{
    assert!(n_threads > 2);
    let buffer_len = n_threads + 2;
    let (sender, receiver) = bounded::<Vec<Base<Vec<u8>>>>(buffer_len);
    let (done_send, done_recv) = bounded::<Option<O>>(buffer_len);
    let receiver = Arc::new(receiver); // 使用 Arc 来共享 receiver
    let done_send = Arc::new(done_send);
    let mut pool = Pool::new(n_threads as u32);

    let mut parallel_result = ParallelResult { recv: done_recv };

    pool.scoped(|pool_scope| {
        // 生产者线程
        pool_scope.execute(move || {
            while let Ok(Some(seqs)) = reader.next() {
                sender.send(seqs).expect("Failed to send sequences");
            }
        });

        // 消费者线程
        for _ in 0..n_threads - 2 {
            let receiver = Arc::clone(&receiver);
            let work = &work;
            let done_send = Arc::clone(&done_send);
            pool_scope.execute(move || {
                while let Ok(mut seqs) = receiver.recv() {
                    let mut markers: Vec<Base<MinimizerIterator<'_>>> = seqs
                        .iter_mut()
                        .map(|seq| scan_sequence(seq, &meros))
                        .collect();
                    let output = work(&mut markers);
                    done_send.send(output).expect("Failed to send outputs");
                }
            });
        }

        // 引用计数减掉一个,这样都子线程结束时, done_send还能完全释放
        drop(done_send);
        pool_scope.execute(move || {
            let _ = func(&mut parallel_result);
        });

        pool_scope.join_all();
    });

    Ok(())
}

pub fn buffer_read_parallel<R, D, W, O, F, Out>(
    reader: &mut R,
    n_threads: usize,
    buffer_size: usize,
    work: W,
    func: F,
) -> Result<()>
where
    D: Send + Sized + Sync,
    R: std::io::Read + Send,
    O: Send,
    Out: Send + Default,
    W: Send + Sync + Fn(&[D]) -> Option<O>,
    F: FnOnce(&mut ParallelResult<Option<O>>) -> Out + Send,
{
    assert!(n_threads > 2);
    let buffer_len = n_threads + 2;
    let (sender, receiver) = bounded::<&[D]>(buffer_len);
    let (done_send, done_recv) = bounded::<Option<O>>(buffer_len);
    let receiver = Arc::new(receiver); // 使用 Arc 来共享 receiver
    let done_send = Arc::new(done_send);
    let mut pool = Pool::new(n_threads as u32);

    let slot_size = std::mem::size_of::<D>();
    let mut parallel_result = ParallelResult { recv: done_recv };

    pool.scoped(|pool_scope| {
        // 生产者线程
        pool_scope.execute(move || {
            let mut batch_buffer = vec![0u8; slot_size * buffer_size];

            while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
                if bytes_read == 0 {
                    break;
                } // 文件末尾

                let slots_in_batch = bytes_read / slot_size;
                let slots = unsafe {
                    std::slice::from_raw_parts(batch_buffer.as_ptr() as *const D, slots_in_batch)
                };
                sender.send(slots).expect("Failed to send sequences");
            }
        });

        // 消费者线程
        for _ in 0..n_threads - 2 {
            let receiver = Arc::clone(&receiver);
            let work = &work;
            let done_send = Arc::clone(&done_send);
            pool_scope.execute(move || {
                while let Ok(seqs) = receiver.recv() {
                    let output = work(seqs);
                    done_send.send(output).expect("Failed to send outputs");
                }
            });
        }

        // 引用计数减掉一个,这样都子线程结束时, done_send还能完全释放
        drop(done_send);
        pool_scope.execute(move || {
            let _ = func(&mut parallel_result);
        });

        pool_scope.join_all();
    });

    Ok(())
}

pub fn buffer_map_parallel<D, W, O, F, Out>(
    map: &HashMap<u32, Vec<D>>,
    n_threads: usize,
    work: W,
    func: F,
) -> Result<()>
where
    D: Send + Sized + Sync,
    O: Send,
    Out: Send + Default,
    W: Send + Sync + Fn((&u32, &Vec<D>)) -> Option<O>,
    F: FnOnce(&mut ParallelResult<Option<O>>) -> Out + Send,
{
    assert!(n_threads > 2);
    let buffer_len = n_threads + 2;
    let (sender, receiver) = bounded::<(&u32, &Vec<D>)>(buffer_len);
    let (done_send, done_recv) = bounded::<Option<O>>(buffer_len);
    let receiver = Arc::new(receiver); // 使用 Arc 来共享 receiver
    let done_send = Arc::new(done_send);
    let mut pool = Pool::new(n_threads as u32);

    let mut parallel_result = ParallelResult { recv: done_recv };

    pool.scoped(|pool_scope| {
        // 生产者线程
        pool_scope.execute(move || {
            for entry in map {
                sender.send(entry).expect("Failed to send sequences");
            }
        });

        // 消费者线程
        for _ in 0..n_threads - 2 {
            let receiver = Arc::clone(&receiver);
            let work = &work;
            let done_send = Arc::clone(&done_send);
            pool_scope.execute(move || {
                while let Ok(seqs) = receiver.recv() {
                    let output = work(seqs);
                    done_send.send(output).expect("Failed to send outputs");
                }
            });
        }

        // 引用计数减掉一个,这样都子线程结束时, done_send还能完全释放
        drop(done_send);
        pool_scope.execute(move || {
            let _ = func(&mut parallel_result);
        });

        pool_scope.join_all();
    });

    Ok(())
}
