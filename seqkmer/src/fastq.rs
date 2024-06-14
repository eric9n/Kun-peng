use crate::reader::{dyn_reader, trim_end, Reader, BUFSIZE};
use crate::seq::{SeqFormat, Sequence};
use std::io::{BufRead, BufReader, Read, Result};
use std::path::Path;

/// FastqReader
pub struct FastqReader<R: Read + Send> {
    pub reader: BufReader<R>,
    pub file_index: u64,
    pub reads_index: u64,
    pub seq_id: String,

    score: i32,
    header: Vec<u8>,
    seq: Vec<u8>,
    plus: Vec<u8>,
    quals: Vec<u8>,
}

impl<R> FastqReader<R>
where
    R: Read + Send,
{
    pub fn new(reader: R, file_index: u64, score: i32) -> Self {
        Self::with_capacity(reader, file_index, BUFSIZE, score)
    }

    pub fn with_capacity<'a>(reader: R, file_index: u64, capacity: usize, score: i32) -> Self {
        assert!(capacity >= 3);
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            file_index,
            reads_index: 0,
            seq_id: String::new(),
            header: Vec::new(),
            seq: Vec::new(),
            plus: Vec::new(),
            quals: Vec::new(),
            score,
        }
    }

    pub fn read_lines(&mut self) -> Result<Option<Sequence>> {
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

        let seq_id = unsafe {
            let s = std::str::from_utf8_unchecked(&self.header[1..]);
            let first_space_index = s
                .as_bytes()
                .iter()
                .position(|&c| c == b' ')
                .unwrap_or(s.len());

            // 直接从原始切片创建第一个单词的切片
            &s[..first_space_index]
        };
        self.reads_index += 1;

        if self.score > 0 {
            for (base, &qscore) in self.seq.iter_mut().zip(self.quals.iter()) {
                if (qscore as i32 - '!' as i32) < self.score {
                    *base = b'x';
                }
            }
        }

        let sequence = Sequence {
            file_index: self.file_index,
            reads_index: self.reads_index,
            id: seq_id.to_owned(),
            seq: self.seq.to_owned(),
            format: SeqFormat::FASTQ,
        };
        Ok(Some(sequence))
    }
}

impl FastqReader<Box<dyn Read + Send>> {
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path: P, file_index: u64, score: i32) -> Result<Self> {
        let reader = dyn_reader(path)?;
        Ok(Self::new(reader, file_index, score))
    }
}

impl<R> Reader<R> for FastqReader<R>
where
    R: Read + Send,
{
    fn next(&mut self) -> Result<Option<Sequence>> {
        self.read_lines()
    }
}
