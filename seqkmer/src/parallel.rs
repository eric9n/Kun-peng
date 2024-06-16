use crate::reader::{Reader, SeqMer};
use crate::seq::Sequence;
use crate::Meros;
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

pub fn read_parallel<R, W, O, F, Out>(
    reader: &mut R,
    n_threads: usize,
    buffer_len: usize,
    meros: Meros,
    work: W,
    func: F,
) -> Result<()>
where
    R: Reader,
    O: Send,
    Out: Send + Default,
    W: Send + Sync + Fn(Vec<SeqMer>) -> Option<O>,
    F: FnOnce(&mut ParallelResult<Option<O>>) -> Out + Send,
{
    assert!(n_threads > 2);
    assert!(n_threads <= buffer_len);
    let (sender, receiver) = bounded::<Vec<Sequence>>(buffer_len);
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
                while let Ok(seqs) = receiver.recv() {
                    let seq_mers: Vec<SeqMer> = seqs
                        .iter()
                        .map(|seq| SeqMer::from_seq(seq, meros))
                        .collect();

                    let output = work(seq_mers);
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
