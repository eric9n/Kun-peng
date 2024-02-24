use seq_io::fastq;
use seq_io::parallel::Reader;

use std::fs::File;
use std::io;
use std::iter;
use std::path::Path;

use seq_io::policy::StdPolicy;

type DefaultBufPolicy = StdPolicy;

pub struct PairReader<R: io::Read, P = DefaultBufPolicy> {
    reader1: fastq::Reader<R, P>,
    reader2: fastq::Reader<R, P>,
}

impl Default for PairRecordSet {
    fn default() -> Self {
        PairRecordSet(fastq::RecordSet::default(), fastq::RecordSet::default())
    }
}

impl PairReader<File, DefaultBufPolicy> {
    /// Creates a reader from a file path.
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path1: P, path2: P) -> io::Result<PairReader<File>> {
        // 分别打开两个文件
        let file1 = File::open(path1)?;
        let file2 = File::open(path2)?;

        // 为每个文件创建一个 fastq::Reader 实例
        let reader1 = fastq::Reader::new(file1);
        let reader2 = fastq::Reader::new(file2);

        // 使用这两个实例构造一个 PairReader 对象
        Ok(PairReader { reader1, reader2 })
    }
}

pub struct PairRecordSet(fastq::RecordSet, fastq::RecordSet);

impl<'a> iter::IntoIterator for &'a PairRecordSet {
    type Item = (fastq::RefRecord<'a>, fastq::RefRecord<'a>);
    type IntoIter = PairRecordSetIter<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        PairRecordSetIter(self.0.into_iter(), self.1.into_iter())
    }
}

pub struct PairRecordSetIter<'a>(fastq::RecordSetIter<'a>, fastq::RecordSetIter<'a>);

impl<'a> Iterator for PairRecordSetIter<'a> {
    type Item = (fastq::RefRecord<'a>, fastq::RefRecord<'a>);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match (self.0.next(), self.1.next()) {
            (Some(record1), Some(record2)) => Some((record1, record2)),
            _ => None, // Return None if either iterator runs out of records
        }
    }
}

impl<R, P> Reader for PairReader<R, P>
where
    R: io::Read,
    P: seq_io::policy::BufPolicy + Send,
{
    type DataSet = PairRecordSet;
    type Err = fastq::Error;

    #[inline]
    fn fill_data(&mut self, rset: &mut PairRecordSet) -> Option<Result<(), Self::Err>> {
        let res1 = self.reader1.read_record_set(&mut rset.0)?.is_err();
        let res2 = self.reader2.read_record_set(&mut rset.1)?.is_err();

        if res1 || res2 {
            return None;
        }

        // If both reads are successful, return Ok(())
        Some(Ok(()))
    }
}
