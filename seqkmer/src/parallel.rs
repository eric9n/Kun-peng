use crate::fasta::FastaReader;
use crate::fastq::{FastqPairReader, FastqReader};
use crate::mmscanner::{scan_sequence, MinimizerIterator};
use crate::reader::{detect_file_format, Reader};
use crate::seq::{BaseType, SeqHeader};
use crate::{Meros, SeqFormat};
use crossbeam_channel::{bounded, Receiver};
use scoped_threadpool::Pool;
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
) -> Result<Box<dyn Reader>> {
    let mut files_iter = file_pair.iter();
    let file1 = files_iter.next().cloned().unwrap();
    let file2 = files_iter.next().cloned();
    match detect_file_format(&file_pair[0])? {
        SeqFormat::Fastq => {
            if let Some(file2) = file2 {
                Ok(Box::new(FastqPairReader::from_path(
                    file1, file2, file_index, score,
                )?))
            } else {
                Ok(Box::new(FastqReader::from_path(file1, file_index, score)?))
            }
        }
        SeqFormat::Fasta => Ok(Box::new(FastaReader::from_path(file1, file_index)?)),

        _ => unreachable!(),
    }
}

pub fn read_parallel<R, W, O, F, Out>(
    reader: &mut R,
    n_threads: usize,
    buffer_len: usize,
    meros: &Meros,
    work: W,
    func: F,
) -> Result<()>
where
    R: Reader,
    O: Send,
    Out: Send + Default,
    W: Send + Sync + Fn(&mut Vec<BaseType<SeqHeader, MinimizerIterator>>) -> Option<O>,
    F: FnOnce(&mut ParallelResult<Option<O>>) -> Out + Send,
{
    assert!(n_threads > 2);
    assert!(n_threads <= buffer_len);
    let (sender, receiver) = bounded::<Vec<BaseType<SeqHeader, Vec<u8>>>>(buffer_len);
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
                    let mut markers: Vec<BaseType<SeqHeader, MinimizerIterator<'_>>> = seqs
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
