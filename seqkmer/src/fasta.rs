use crate::reader::{dyn_reader, trim_end, Reader, BUFSIZE};
use crate::seq::{Base, SeqFormat, SeqHeader};
use crate::utils::OptionPair;
use std::io::{BufRead, BufReader, Read, Result};
use std::path::Path;

const SEQ_LIMIT: u64 = u64::pow(2, 32);

/// FastaReader
pub struct FastaReader<R>
where
    R: Read + Send,
{
    reader: BufReader<R>,
    file_index: usize,
    reads_index: usize,
    header: Vec<u8>,
    seq: Vec<u8>,

    // 批量读取
    batch_size: usize,
}

impl<R> FastaReader<R>
where
    R: Read + Send,
{
    pub fn new(reader: R, file_index: usize) -> Self {
        Self::with_capacity(reader, file_index, BUFSIZE, 30)
    }

    pub fn with_capacity(reader: R, file_index: usize, capacity: usize, batch_size: usize) -> Self {
        assert!(capacity >= 3);
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            file_index,
            reads_index: 0,
            header: Vec::new(),
            seq: Vec::new(),
            batch_size,
        }
    }

    pub fn read_next(&mut self) -> Result<Option<()>> {
        // 读取fastq文件header部分
        self.header.clear();
        if self.reader.read_until(b'\n', &mut self.header)? == 0 {
            return Ok(None);
        }
        // 读取fasta文件seq部分
        self.seq.clear();
        if self.reader.read_until(b'>', &mut self.seq)? == 0 {
            return Ok(None);
        }
        trim_end(&mut self.seq);
        Ok(Some(()))
    }

    pub fn _next(&mut self) -> Result<Option<(usize, Base<Vec<u8>>)>> {
        if self.read_next()?.is_none() {
            return Ok(None);
        }

        let seq_len = self.seq.len();
        // 检查seq的长度是否大于2的32次方
        if seq_len as u64 > SEQ_LIMIT {
            eprintln!("Sequence length exceeds 2^32, which is not handled.");
            return Ok(None);
        }

        let seq_id = unsafe {
            let slice = if self.header.starts_with(b">") {
                &self.header[1..]
            } else {
                &self.header[..]
            };

            let s = std::str::from_utf8_unchecked(slice);
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
        Ok(Some((
            seq_len,
            Base::new(seq_header, OptionPair::Single(self.seq.to_owned())),
        )))
    }
}

impl FastaReader<Box<dyn Read + Send>> {
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path: P, file_index: usize) -> Result<Self> {
        let reader = dyn_reader(path)?;
        Ok(Self::new(reader, file_index))
    }
}

impl<R: Read + Send> Reader for FastaReader<R> {
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        let mut seqs = Vec::new();
        let mut total_bytes = 0;
        let max_bytes = 10 * 1024 * 1024;

        for _ in 0..self.batch_size {
            if let Some((seq_len, seq)) = self._next()? {
                seqs.push(seq);
                total_bytes += seq_len;
                if total_bytes > max_bytes {
                    break;
                }
            } else {
                break;
            }
        }

        Ok(if seqs.is_empty() { None } else { Some(seqs) })
    }
}

/// BufferFastaReader
pub struct BufferFastaReader<R>
where
    R: Read + Send,
{
    reader: BufReader<R>,
    file_index: usize,
    reads_index: usize,
    header: Vec<u8>,
    seq: Vec<u8>,
    line_num: usize,

    // 批量读取
    batch_size: usize,
}

impl<R> BufferFastaReader<R>
where
    R: Read + Send,
{
    pub fn new(reader: R, file_index: usize) -> Self {
        Self::with_capacity(reader, file_index, BUFSIZE, 60)
    }

    pub fn with_capacity(reader: R, file_index: usize, capacity: usize, batch_size: usize) -> Self {
        assert!(capacity >= 3);
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            file_index,
            reads_index: 0,
            line_num: 0,
            header: Vec::new(),
            seq: Vec::new(),
            batch_size,
        }
    }

    pub fn read_next(&mut self) -> Result<Option<()>> {
        // 读取fastq文件header部分
        if self.header.is_empty() {
            if self.reader.read_until(b'\n', &mut self.header)? == 0 {
                return Ok(None);
            }
        }

        if self.reader.read_until(b'\n', &mut self.seq)? == 0 {
            return Ok(None);
        }
        if self.seq.starts_with(&[b'>']) {
            self.header = self.seq.clone();
            self.seq.clear();
            if self.reader.read_until(b'\n', &mut self.seq)? == 0 {
                return Ok(None);
            }
        }
        self.line_num += 1;
        trim_end(&mut self.seq);
        Ok(Some(()))
    }

    pub fn _next(&mut self) -> Result<Option<Base<Vec<u8>>>> {
        self.seq.clear();
        for _ in 0..self.batch_size {
            if self.read_next()?.is_none() {
                return Ok(None);
            }
        }

        let seq_len = self.seq.len();
        // 检查seq的长度是否大于2的32次方
        if seq_len as u64 > SEQ_LIMIT {
            eprintln!("Sequence length exceeds 2^32, which is not handled.");
            return Ok(None);
        }

        let seq_id = unsafe {
            let slice = if self.header.starts_with(b">") {
                &self.header[1..]
            } else {
                &self.header[..]
            };

            let s = std::str::from_utf8_unchecked(slice);
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
        Ok(Some(Base::new(
            seq_header,
            OptionPair::Single(self.seq.to_owned()),
        )))
    }
}

impl BufferFastaReader<Box<dyn Read + Send>> {
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path: P, file_index: usize) -> Result<Self> {
        let reader = dyn_reader(path)?;
        Ok(Self::new(reader, file_index))
    }
}

impl<R: Read + Send> Reader for BufferFastaReader<R> {
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        let mut seqs = Vec::new();
        if let Some(seq) = self._next()? {
            seqs.push(seq);
        }

        Ok(if seqs.is_empty() { None } else { Some(seqs) })
    }
}
