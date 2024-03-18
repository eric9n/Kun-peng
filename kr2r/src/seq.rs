use crate::mmscanner::MinimizerScanner;
use seq_io::fasta;
use seq_io::fasta::Record as FaRecord;
use seq_io::fastq;
use seq_io::fastq::Record as FqRecord;

use seq_io::parallel::Reader;

use crate::utils::is_gzipped;
use seq_io::policy::StdPolicy;
use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::iter;
use std::path::Path;

use crate::Meros;

type DefaultBufPolicy = StdPolicy;

pub trait SeqX {
    fn seq_x(&self, score: i32) -> Vec<u8>;
}

impl<'a> SeqX for fastq::RefRecord<'a> {
    fn seq_x(&self, score: i32) -> Vec<u8> {
        if score <= 0 {
            return self.seq().to_vec();
        }

        let qual = self.qual();
        self.seq()
            .iter()
            .zip(qual.iter())
            .map(|(&base, &qscore)| {
                if (qscore as i32 - '!' as i32) < score {
                    b'x'
                } else {
                    base
                }
            })
            .collect::<Vec<u8>>()
    }
}

impl SeqX for fastq::OwnedRecord {
    fn seq_x(&self, score: i32) -> Vec<u8> {
        if score <= 0 {
            return self.seq().to_vec();
        }
        let qual = self.qual();
        self.seq()
            .iter()
            .zip(qual.iter())
            .map(|(&base, &qscore)| {
                if (qscore as i32 - '!' as i32) < score {
                    b'x'
                } else {
                    base
                }
            })
            .collect::<Vec<u8>>()
    }
}

impl<'a> SeqX for fasta::RefRecord<'a> {
    #[allow(unused_variables)]
    fn seq_x(&self, score: i32) -> Vec<u8> {
        self.seq().to_vec()
    }
}

#[derive(Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct SeqReads {
    pub dna_id: String,
    pub seq_paired: Vec<Vec<u64>>,
}

pub trait SeqSet {
    fn to_seq_reads(&self, score: i32, meros: Meros) -> HashSet<SeqReads>;
}

pub struct PairFastqReader<P = DefaultBufPolicy> {
    reader1: fastq::Reader<Box<dyn io::Read + Send>, P>,
    reader2: fastq::Reader<Box<dyn io::Read + Send>, P>,
    index: usize, // 新增索引字段
}

impl Default for PairFastqRecordSet {
    fn default() -> Self {
        PairFastqRecordSet(fastq::RecordSet::default(), fastq::RecordSet::default())
    }
}

use flate2::read::GzDecoder;

impl<'a> PairFastqReader<DefaultBufPolicy> {
    /// Creates a reader from a file path.
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path1: P, path2: P) -> io::Result<PairFastqReader> {
        // 分别打开两个文件
        let mut file1 = File::open(&path1)?;
        let mut file2 = File::open(path2)?;

        let read1: Box<dyn io::Read + Send> = if is_gzipped(&mut file1)? {
            Box::new(GzDecoder::new(file1))
        } else {
            Box::new(file1)
        };

        let read2: Box<dyn io::Read + Send> = if is_gzipped(&mut file2)? {
            Box::new(GzDecoder::new(file2))
        } else {
            Box::new(file2)
        };

        // 为每个文件创建一个 fastq::Reader 实例
        let reader1 = fastq::Reader::new(read1);
        let reader2 = fastq::Reader::new(read2);

        // 使用这两个实例构造一个 PairFastqReader 对象
        Ok(PairFastqReader {
            reader1,
            reader2,
            index: 0,
        })
    }

    pub fn next(&mut self) -> Option<PairFastqRecord> {
        let ref_record1 = self
            .reader1
            .next()?
            .expect("fastq file error")
            .to_owned_record();
        let ref_recrod2 = self
            .reader2
            .next()?
            .expect("fastq file error")
            .to_owned_record();

        self.index += 1;
        Some(PairFastqRecord(self.index, ref_record1, ref_recrod2))
    }
}

pub struct PairFastqRecord(pub usize, pub fastq::OwnedRecord, pub fastq::OwnedRecord);

pub struct PairFastqRecordSet(fastq::RecordSet, fastq::RecordSet);

impl<'a> iter::IntoIterator for &'a PairFastqRecordSet {
    type Item = (fastq::RefRecord<'a>, fastq::RefRecord<'a>);
    type IntoIter = PairFastqRecordSetIter<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        PairFastqRecordSetIter(self.0.into_iter(), self.1.into_iter())
    }
}

pub struct PairFastqRecordSetIter<'a>(fastq::RecordSetIter<'a>, fastq::RecordSetIter<'a>);

impl<'a> Iterator for PairFastqRecordSetIter<'a> {
    type Item = (fastq::RefRecord<'a>, fastq::RefRecord<'a>);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match (self.0.next(), self.1.next()) {
            (Some(record1), Some(record2)) => Some((record1, record2)),
            _ => None, // Return None if either iterator runs out of records
        }
    }
}

impl<P> Reader for PairFastqReader<P>
where
    P: seq_io::policy::BufPolicy + Send,
{
    type DataSet = PairFastqRecordSet;
    type Err = fastq::Error;

    #[inline]
    fn fill_data(&mut self, rset: &mut PairFastqRecordSet) -> Option<Result<(), Self::Err>> {
        let res1 = self.reader1.read_record_set(&mut rset.0)?.is_err();
        let res2 = self.reader2.read_record_set(&mut rset.1)?.is_err();

        if res1 || res2 {
            return None;
        }

        // If both reads are successful, return Ok(())
        Some(Ok(()))
    }
}

impl SeqSet for PairFastqRecordSet {
    fn to_seq_reads(&self, score: i32, meros: Meros) -> HashSet<SeqReads> {
        let mut seq_pair_set = HashSet::<SeqReads>::new();

        for records in self.into_iter() {
            let dna_id = records.0.id().unwrap_or_default().to_string();
            let seq1 = records.0.seq_x(score);
            let seq2 = records.1.seq_x(score);

            let kmers1 = MinimizerScanner::new(&seq1, meros).collect();
            let kmers2 = MinimizerScanner::new(&seq2, meros).collect();

            let seq_paired: Vec<Vec<u64>> = vec![kmers1, kmers2];
            seq_pair_set.insert(SeqReads { dna_id, seq_paired });
        }
        seq_pair_set
    }
}

impl SeqSet for fastq::RecordSet {
    fn to_seq_reads(&self, score: i32, meros: Meros) -> HashSet<SeqReads> {
        let mut seq_pair_set = HashSet::<SeqReads>::new();
        for records in self.into_iter() {
            let dna_id = records.id().unwrap_or_default().to_string();
            let seq1 = records.seq_x(score);
            let kmers1: Vec<u64> = MinimizerScanner::new(&seq1, meros).collect();
            let seq_paired: Vec<Vec<u64>> = vec![kmers1];
            seq_pair_set.insert(SeqReads { dna_id, seq_paired });
        }
        seq_pair_set
    }
}

impl SeqSet for fasta::RecordSet {
    fn to_seq_reads(&self, score: i32, meros: Meros) -> HashSet<SeqReads> {
        let mut seq_pair_set = HashSet::<SeqReads>::new();
        for records in self.into_iter() {
            let dna_id = records.id().unwrap_or_default().to_string();
            let seq1 = records.seq_x(score);
            let kmers1: Vec<u64> = MinimizerScanner::new(&seq1, meros).collect();
            let seq_paired: Vec<Vec<u64>> = vec![kmers1];
            seq_pair_set.insert(SeqReads { dna_id, seq_paired });
        }
        seq_pair_set
    }
}

// pub struct PairFastaRecordSet(fasta::RecordSet, fasta::RecordSet);

// impl Default for PairFastaRecordSet {
//     fn default() -> Self {
//         PairFastaRecordSet(fasta::RecordSet::default(), fasta::RecordSet::default())
//     }
// }

// pub struct PairFastaReader<R: io::Read, P = DefaultBufPolicy> {
//     reader1: fasta::Reader<R, P>,
//     reader2: fasta::Reader<R, P>,
// }

// impl PairFastaReader<File, DefaultBufPolicy> {
//     #[inline]
//     pub fn from_path<P: AsRef<Path>>(path1: P, path2: P) -> io::Result<PairFastaReader<File>> {
//         let file1 = File::open(path1)?;
//         let file2 = File::open(path2)?;

//         let reader1 = fasta::Reader::new(file1);
//         let reader2 = fasta::Reader::new(file2);

//         Ok(PairFastaReader { reader1, reader2 })
//     }
// }

// impl<'a> iter::IntoIterator for &'a PairFastaRecordSet {
//     type Item = (fasta::RefRecord<'a>, fasta::RefRecord<'a>);
//     type IntoIter = PairFastaRecordSetIter<'a>;

//     #[inline]
//     fn into_iter(self) -> Self::IntoIter {
//         PairFastaRecordSetIter(self.0.into_iter(), self.1.into_iter())
//     }
// }

// pub struct PairFastaRecordSetIter<'a>(fasta::RecordSetIter<'a>, fasta::RecordSetIter<'a>);

// impl<'a> Iterator for PairFastaRecordSetIter<'a> {
//     type Item = (fasta::RefRecord<'a>, fasta::RefRecord<'a>);

//     #[inline]
//     fn next(&mut self) -> Option<Self::Item> {
//         match (self.0.next(), self.1.next()) {
//             (Some(record1), Some(record2)) => Some((record1, record2)),
//             _ => None,
//         }
//     }
// }

// impl<R, P> Reader for PairFastaReader<R, P>
// where
//     R: io::Read,
//     P: seq_io::policy::BufPolicy + Send,
// {
//     type DataSet = PairFastaRecordSet;
//     type Err = fasta::Error;

//     #[inline]
//     fn fill_data(&mut self, rset: &mut PairFastaRecordSet) -> Option<Result<(), Self::Err>> {
//         let res1 = self.reader1.read_record_set(&mut rset.0)?.is_err();
//         let res2 = self.reader2.read_record_set(&mut rset.1)?.is_err();

//         if res1 || res2 {
//             return None;
//         }

//         // If both reads are successful, return Ok(())
//         Some(Ok(()))
//     }
// }

// impl SeqSet for PairFastaRecordSet {
//     fn to_seq_reads(&self, score: i32, meros: Meros) -> HashSet<SeqReads> {
//         let mut seq_pair_set = HashSet::<SeqReads>::new();

//         for records in self.into_iter() {
//             let dna_id = records.0.id().unwrap_or_default().to_string();
//             let seq1 = records.0.seq_x(score);
//             let seq2 = records.1.seq_x(score);

//             let kmers1 = KmerIterator::new(&seq1, meros).collect();
//             let kmers2 = KmerIterator::new(&seq2, meros).collect();

//             let seq_paired: Vec<Vec<u64>> = vec![kmers1, kmers2];
//             seq_pair_set.insert(SeqReads { dna_id, seq_paired });
//         }
//         seq_pair_set
//     }
// }
