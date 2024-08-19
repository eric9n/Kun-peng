use crate::reader::{dyn_reader, trim_end, trim_pair_info, Reader, BUFSIZE};
use crate::seq::{Base, SeqFormat, SeqHeader};
use crate::utils::OptionPair;
use std::io::{BufRead, BufReader, Read, Result};
use std::path::Path;

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

pub struct FastqReader<R: Read + Send> {
    inner: OptionPair<QReader<R>>,
    file_index: usize,
    reads_index: usize,
    // 批量读取
    batch_size: usize,
}

impl<R> FastqReader<R>
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
                .find(|c: char| c.is_whitespace() || c == '\u{1}')
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

    pub fn read_next(&mut self) -> Result<Option<Base<Vec<u8>>>> {
        match &mut self.inner {
            OptionPair::Single(reader) => {
                if reader.read_next()?.is_none() {
                    return Ok(None);
                }

                self.reads_index += 1;

                let seq_header =
                    Self::create_seq_header(&reader, self.file_index, self.reads_index);
                Ok(Some(Base::new(
                    seq_header,
                    OptionPair::Single(reader.seq.to_owned()),
                )))
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

                Ok(Some(Base::new(
                    seq_header,
                    OptionPair::Pair(reader1.seq.to_owned(), reader2.seq.to_owned()),
                )))
            }
        }
    }
}

impl FastqReader<Box<dyn Read + Send>> {
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

impl<R> Reader for FastqReader<R>
where
    R: Read + Send,
{
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        let seqs: Vec<Base<Vec<u8>>> = (0..self.batch_size)
            .filter_map(|_| self.read_next().transpose())
            .collect::<Result<Vec<_>>>()?;

        Ok(Some(seqs).filter(|v| !v.is_empty()))
    }
}
