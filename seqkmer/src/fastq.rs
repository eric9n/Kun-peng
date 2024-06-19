use crate::reader::{dyn_reader, trim_end, trim_pair_info, Reader, SeqVecType, BUFSIZE};
use crate::seq::{BaseType, SeqFormat, SeqHeader};
use std::io::{BufRead, BufReader, Read, Result};
use std::path::Path;

type SeqType = BaseType<SeqHeader, Vec<u8>>;

struct QReader<R: Read + Send> {
    reader: BufReader<R>,
    quality_score: i32,

    header: Vec<u8>,
    seq: Vec<u8>,
    plus: Vec<u8>,
    quals: Vec<u8>,
}

impl<R> QReader<R>
where
    R: Read + Send,
{
    pub fn with_capacity(reader: R, capacity: usize, quality_score: i32) -> Self {
        assert!(capacity >= 3);
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            header: Vec::new(),
            seq: Vec::new(),
            plus: Vec::new(),
            quals: Vec::new(),
            quality_score,
        }
    }

    pub fn read_next_entry(&mut self) -> Result<Option<(Vec<u8>, Vec<u8>)>> {
        // 读取fastq文件header部分
        let mut header: Vec<u8> = Vec::new();
        if self.reader.read_until(b'\n', &mut header)? == 0 {
            return Ok(None);
        }
        // 读取fastq文件seq部分
        let mut seq = Vec::new();
        if self.reader.read_until(b'\n', &mut seq)? == 0 {
            return Ok(None);
        }
        trim_end(&mut seq);

        // 读取fastq文件+部分
        self.plus.clear();
        if self.reader.read_until(b'\n', &mut self.plus)? == 0 {
            return Ok(None);
        }

        // 读取fastq文件quals部分
        self.quals.clear();
        if self.reader.read_until(b'\n', &mut self.quals)? == 0 {
            return Ok(None);
        }
        trim_end(&mut self.quals);

        if self.quality_score > 0 {
            for (base, &qscore) in self.seq.iter_mut().zip(self.quals.iter()) {
                if (qscore as i32 - '!' as i32) < self.quality_score {
                    *base = b'x';
                }
            }
        }

        Ok(Some((header, seq)))
    }

    pub fn read_next(&mut self) -> Result<Option<()>> {
        // 读取fastq文件header部分
        self.header.clear();
        if self.reader.read_until(b'\n', &mut self.header)? == 0 {
            return Ok(None);
        }
        // 读取fastq文件seq部分
        self.seq.clear();
        if self.reader.read_until(b'\n', &mut self.seq)? == 0 {
            return Ok(None);
        }
        trim_end(&mut self.seq);

        // 读取fastq文件+部分
        self.plus.clear();
        if self.reader.read_until(b'\n', &mut self.plus)? == 0 {
            return Ok(None);
        }

        // 读取fastq文件quals部分
        self.quals.clear();
        if self.reader.read_until(b'\n', &mut self.quals)? == 0 {
            return Ok(None);
        }
        trim_end(&mut self.quals);

        if self.quality_score > 0 {
            for (base, &qscore) in self.seq.iter_mut().zip(self.quals.iter()) {
                if (qscore as i32 - '!' as i32) < self.quality_score {
                    *base = b'x';
                }
            }
        }

        Ok(Some(()))
    }
}

/// FastqReader
pub struct FastqReader<R: Read + Send> {
    inner: QReader<R>,
    file_index: usize,
    reads_index: usize,
    // 批量读取
    batch_size: usize,
}

impl<R> FastqReader<R>
where
    R: Read + Send,
{
    pub fn new(reader: R, file_index: usize, quality_score: i32) -> Self {
        Self::with_capacity(reader, file_index, BUFSIZE, quality_score, 30)
    }

    pub fn with_capacity<'a>(
        reader: R,
        file_index: usize,
        capacity: usize,
        quality_score: i32,
        batch_size: usize,
    ) -> Self {
        assert!(capacity >= 3);
        Self {
            inner: QReader::with_capacity(reader, capacity, quality_score),
            file_index,
            reads_index: 0,
            batch_size,
        }
    }

    pub fn read_next_entry(&mut self) -> Result<Option<SeqType>> {
        let entry = self.inner.read_next_entry()?;
        if entry.is_none() {
            return Ok(None);
        }

        let (header, seq) = entry.unwrap();

        let seq_id = unsafe {
            let s = std::str::from_utf8_unchecked(&header[1..]);
            let first_space_index = s
                .as_bytes()
                .iter()
                .position(|&c| c == b' ')
                .unwrap_or(s.len());

            // 直接从原始切片创建第一个单词的切片
            &s[..first_space_index]
        };
        self.reads_index += 1;

        let seq_header = SeqHeader {
            file_index: self.file_index,
            reads_index: self.reads_index,
            format: SeqFormat::Fasta,
            id: seq_id.to_owned(),
        };

        let seq = BaseType::Single(seq_header, seq);
        Ok(Some(seq))
    }

    pub fn read_next(&mut self) -> Result<Option<SeqType>> {
        if self.inner.read_next()?.is_none() {
            return Ok(None);
        }

        let seq_id = unsafe {
            let s = std::str::from_utf8_unchecked(&self.inner.header[1..]);
            let first_space_index = s
                .as_bytes()
                .iter()
                .position(|&c| c == b' ')
                .unwrap_or(s.len());

            // 直接从原始切片创建第一个单词的切片
            &s[..first_space_index]
        };
        self.reads_index += 1;

        let seq_header = SeqHeader {
            file_index: self.file_index,
            reads_index: self.reads_index,
            format: SeqFormat::Fasta,
            id: seq_id.to_owned(),
        };

        let seq = BaseType::Single(seq_header, self.inner.seq.to_owned());
        Ok(Some(seq))
    }
}

impl FastqReader<Box<dyn Read + Send>> {
    #[inline]
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        file_index: usize,
        quality_score: i32,
    ) -> Result<Self> {
        let reader = dyn_reader(path)?;
        Ok(Self::new(reader, file_index, quality_score))
    }
}

impl<R> Reader for FastqReader<R>
where
    R: Read + Send,
{
    fn next(&mut self) -> Result<Option<SeqVecType>> {
        let seqs: SeqVecType = (0..self.batch_size)
            .filter_map(|_| self.read_next_entry().transpose()) // 将 Result<Option<_>, _> 转换为 Option<Result<_, _>>
            .collect::<Result<Vec<_>>>()?;

        Ok(Some(seqs).filter(|v| !v.is_empty()))
    }
}

/// FastqPairReader
pub struct FastqPairReader<R: Read + Send> {
    inner1: QReader<R>,
    inner2: QReader<R>,
    file_index: usize,
    reads_index: usize,
    // 批量读取
    batch_size: usize,
}

impl<R> FastqPairReader<R>
where
    R: Read + Send,
{
    pub fn new(reader1: R, reader2: R, file_index: usize, score: i32) -> Self {
        Self::with_capacity(reader1, reader2, file_index, BUFSIZE, score, 30)
    }

    pub fn with_capacity<'a>(
        reader1: R,
        reader2: R,
        file_index: usize,
        capacity: usize,
        score: i32,
        batch_size: usize,
    ) -> Self {
        assert!(capacity >= 3);
        Self {
            inner1: QReader::with_capacity(reader1, capacity, score),
            inner2: QReader::with_capacity(reader2, capacity, score),
            file_index,
            reads_index: 0,
            batch_size,
        }
    }

    pub fn read_next(&mut self) -> Result<Option<SeqType>> {
        if self.inner1.read_next()?.is_none() {
            return Ok(None);
        }

        if self.inner2.read_next()?.is_none() {
            return Ok(None);
        }

        let seq_id = unsafe {
            let s = std::str::from_utf8_unchecked(&self.inner1.header[1..]);
            let first_space_index = s
                .as_bytes()
                .iter()
                .position(|&c| c == b' ')
                .unwrap_or(s.len());

            // 直接从原始切片创建第一个单词的切片
            &s[..first_space_index]
        };
        self.reads_index += 1;

        let seq_header = SeqHeader {
            file_index: self.file_index,
            reads_index: self.reads_index,
            format: SeqFormat::Fasta,
            id: trim_pair_info(seq_id),
        };

        let sequence = BaseType::Pair(
            seq_header,
            self.inner1.seq.to_owned(),
            self.inner2.seq.to_owned(),
        );
        Ok(Some(sequence))
    }
}

impl FastqPairReader<Box<dyn Read + Send>> {
    #[inline]
    pub fn from_path<P: AsRef<Path>>(
        path1: P,
        path2: P,
        file_index: usize,
        quality_score: i32,
    ) -> Result<Self> {
        let reader1 = dyn_reader(path1)?;
        let reader2 = dyn_reader(path2)?;
        Ok(Self::new(reader1, reader2, file_index, quality_score))
    }
}

impl<R> Reader for FastqPairReader<R>
where
    R: Read + Send,
{
    fn next(&mut self) -> Result<Option<SeqVecType>> {
        let seqs: SeqVecType = (0..self.batch_size)
            .filter_map(|_| self.read_next().transpose()) // 将 Result<Option<_>, _> 转换为 Option<Result<_, _>>
            .collect::<Result<Vec<_>>>()?;

        Ok(Some(seqs).filter(|v| !v.is_empty()))
    }
}

use crate::seq::OptionPair;

pub struct FastxPairReader<R: Read + Send> {
    inner: OptionPair<QReader<R>>,
    file_index: usize,
    reads_index: usize,
    // 批量读取
    batch_size: usize,
}

impl<R> FastxPairReader<R>
where
    R: Read + Send,
{
    pub fn new(readers: OptionPair<R>, file_index: usize, quality_score: i32) -> Self {
        Self::with_capacity(readers, file_index, BUFSIZE, quality_score, 30)
    }

    pub fn with_capacity<'a>(
        readers: OptionPair<R>,
        file_index: usize,
        capacity: usize,
        quality_score: i32,
        batch_size: usize,
    ) -> Self {
        assert!(capacity >= 3);
        let inner = match readers {
            OptionPair::Single(reader) => {
                OptionPair::Single(QReader::with_capacity(reader, capacity, quality_score))
            }
            OptionPair::Pair(reader1, reader2) => OptionPair::Pair(
                QReader::with_capacity(reader1, capacity, quality_score),
                QReader::with_capacity(reader2, capacity, quality_score),
            ),
        };
        Self {
            inner,
            file_index,
            reads_index: 0,
            batch_size,
        }
    }

    fn create_seq_header(reader: &QReader<R>, file_index: usize, reads_index: usize) -> SeqHeader {
        let seq_id = unsafe {
            let s = std::str::from_utf8_unchecked(&reader.header[1..]);
            let first_space_index = s
                .as_bytes()
                .iter()
                .position(|&c| c == b' ')
                .unwrap_or(s.len());

            // 直接从原始切片创建第一个单词的切片
            &s[..first_space_index]
        };
        SeqHeader {
            file_index,
            reads_index,
            format: SeqFormat::Fastq,
            id: trim_pair_info(seq_id),
        }
    }

    pub fn read_next(&mut self) -> Result<Option<SeqType>> {
        match &mut self.inner {
            OptionPair::Single(reader) => {
                if reader.read_next()?.is_none() {
                    return Ok(None);
                }

                self.reads_index += 1;

                let seq_header =
                    Self::create_seq_header(&reader, self.file_index, self.reads_index);
                Ok(Some(BaseType::Single(seq_header, reader.seq.to_owned())))
            }
            OptionPair::Pair(reader1, reader2) => {
                if reader1.read_next()?.is_none() {
                    return Ok(None);
                }
                if reader2.read_next()?.is_none() {
                    return Ok(None);
                }

                self.reads_index += 1;
                let seq_header =
                    Self::create_seq_header(&reader1, self.file_index, self.reads_index);

                Ok(Some(BaseType::Pair(
                    seq_header,
                    reader1.seq.to_owned(),
                    reader2.seq.to_owned(),
                )))
            }
        }
    }
}

impl FastxPairReader<Box<dyn Read + Send>> {
    #[inline]
    pub fn from_path<P: AsRef<Path>>(
        paths: OptionPair<P>,
        file_index: usize,
        quality_score: i32,
    ) -> Result<Self> {
        let readers = paths.map(|path| dyn_reader(path))?;
        Ok(Self::new(readers, file_index, quality_score))
    }
}

impl<R> Reader for FastxPairReader<R>
where
    R: Read + Send,
{
    fn next(&mut self) -> Result<Option<SeqVecType>> {
        let seqs: SeqVecType = (0..self.batch_size)
            .filter_map(|_| self.read_next().transpose()) // 将 Result<Option<_>, _> 转换为 Option<Result<_, _>>
            .collect::<Result<Vec<_>>>()?;

        Ok(Some(seqs).filter(|v| !v.is_empty()))
    }
}
