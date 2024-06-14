use crate::mmscanner::MinimizerScanner;
use crate::reader::Reader;
use crate::seq::Sequence;
use crate::Meros;
use crossbeam_channel::bounded;
use scoped_threadpool::Pool;
use std::fs::File;
use std::io::Read;
use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

pub fn read_parallel<R, W>(
    reader: &mut dyn Reader<R>,
    n_threads: usize,
    buffer_len: usize,
    output_file: Option<&PathBuf>,
    meros: Meros,
    work: W,
) -> Result<()>
where
    R: Read + Send,
    W: Send + Sync + Fn(Vec<u64>, Sequence) -> Option<String>,
{
    assert!(n_threads <= buffer_len);
    let (sender, receiver) = bounded::<Sequence>(buffer_len);
    let receiver = Arc::new(receiver); // 使用 Arc 来共享 receiver
    let mut pool = Pool::new(10);

    let counter = Arc::new(AtomicUsize::new(0));

    let mut writer: Box<dyn Write + Send> = match output_file {
        Some(file_name) => {
            let file = File::create(file_name)?;
            Box::new(BufWriter::new(file)) as Box<dyn Write + Send>
        }
        None => Box::new(io::stdout()) as Box<dyn Write + Send>,
    };

    let _ = pool.scoped(|pool_scope| -> Result<()> {
        // 生产者线程
        pool_scope.execute(move || {
            while let Some(seq) = reader.next().unwrap() {
                sender.send(seq).unwrap();
            }
        });

        // 消费者线程
        for i in 0..n_threads {
            let receiver = Arc::clone(&receiver);
            let counter_clone = Arc::clone(&counter);
            let work = &work;

            let mut temp_writer: Box<dyn Write + Send> = match output_file {
                Some(file_name) => {
                    let parent_dir = file_name.parent().unwrap_or_else(|| Path::new(""));
                    let file_name = file_name.file_name().unwrap().to_str().unwrap();
                    let filename = parent_dir.join(format!("{}.tmp.{}", file_name, i));
                    let file = File::create(filename)?;
                    Box::new(BufWriter::new(file)) as Box<dyn Write + Send>
                }
                None => Box::new(io::stdout()) as Box<dyn Write + Send>,
            };
            pool_scope.execute(move || {
                while let Ok(seq) = receiver.recv() {
                    counter_clone.fetch_add(1, Ordering::Relaxed);
                    let mut ms = MinimizerScanner::new(&seq.seq, meros);
                    let res = ms.iter();
                    if let Some(out) = work(res, seq) {
                        temp_writer
                            .write_all(out.as_bytes())
                            .expect("write data error");
                    }
                }
            });
        }
        pool_scope.join_all();
        Ok(())
    });
    println!("counter {:?}", counter.load(Ordering::Relaxed));
    Ok(())
}
