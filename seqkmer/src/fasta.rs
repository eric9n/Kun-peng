use crate::reader::{dyn_reader, trim_end, Reader, SeqVecType, BUFSIZE};
use crate::seq::{self, BaseType, SeqFormat, SeqHeader};
use std::io::{BufRead, BufReader, Read, Result};
use std::path::Path;

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
}

impl<R> FastaReader<R>
where
    R: Read + Send,
{
    pub fn new(reader: R, file_index: usize) -> Self {
        Self::with_capacity(reader, file_index, BUFSIZE)
    }

    pub fn with_capacity(reader: R, file_index: usize, capacity: usize) -> Self {
        assert!(capacity >= 3);
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            file_index,
            reads_index: 0,
            header: Vec::new(),
            seq: Vec::new(),
        }
    }

    pub fn read_next_entry(&mut self) -> Result<Option<(Vec<u8>, Vec<u8>)>> {
        // 读取fastq文件header部分
        let mut header = Vec::new();
        if self.reader.read_until(b'\n', &mut header)? == 0 {
            return Ok(None);
        }
        // 读取fasta文件seq部分
        let mut seq = Vec::new();
        if self.reader.read_until(b'>', &mut seq)? == 0 {
            return Ok(None);
        }
        trim_end(&mut self.seq);
        Ok(Some((header, seq)))
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
}

impl FastaReader<Box<dyn Read + Send>> {
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path: P, file_index: usize) -> Result<Self> {
        let reader = dyn_reader(path)?;
        Ok(Self::new(reader, file_index))
    }
}

fn check_sequence_length(seq: &Vec<u8>) -> bool {
    let limit = u64::pow(2, 32);
    // 检查seq的长度是否大于2的32次方
    (seq.len() as u64) > limit
}

impl<R: Read + Send> Reader for FastaReader<R> {
    fn next(&mut self) -> Result<Option<SeqVecType>> {
        // if self.read_next()?.is_none() {
        //     return Ok(None);
        // }

        let entry = self.read_next_entry()?;
        if entry.is_none() {
            return Ok(None);
        }
        let (header, seq) = entry.unwrap();
        if check_sequence_length(&seq) {
            eprintln!("Sequence length exceeds 2^32, which is not handled.");
            return Ok(None);
        }

        let seq_id = unsafe {
            let slice = if header.starts_with(b">") {
                &header[1..]
            } else {
                &header[..]
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
        let seq = BaseType::Single(seq_header, seq);
        Ok(Some(vec![seq]))
    }
}
