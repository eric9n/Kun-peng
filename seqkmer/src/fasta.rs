use crate::reader::{dyn_reader, trim_end, Reader, BUFSIZE};
use crate::seq::{BaseType, SeqFormat, Sequence};
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
    fn next(&mut self) -> Result<Option<Vec<Sequence>>> {
        if self.read_next()?.is_none() {
            return Ok(None);
        }

        if check_sequence_length(&self.seq) {
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

        let sequence = Sequence {
            file_index: self.file_index,
            reads_index: self.reads_index,
            id: seq_id.to_owned(),
            seq: BaseType::Single(self.seq.to_owned()),
            format: SeqFormat::Fasta,
        };
        Ok(Some(vec![sequence]))
    }
}
