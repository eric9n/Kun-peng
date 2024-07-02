use crate::seq::{Base, SeqFormat};
use crate::utils::OptionPair;
use flate2::read::GzDecoder;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Result, Seek};
use std::path::Path;

pub(crate) fn dyn_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn Read + Send>> {
    let mut file = open_file(path)?;
    if is_gzipped(&mut file)? {
        let decoder = GzDecoder::new(file);
        Ok(Box::new(decoder))
    } else {
        Ok(Box::new(file))
    }
}

pub(crate) fn is_gzipped(file: &mut File) -> Result<bool> {
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;
    file.rewind()?; // 重置文件指针到开头
    Ok(buffer == [0x1F, 0x8B])
}

pub fn trim_pair_info(id: &str) -> String {
    let sz = id.len();
    if sz <= 2 {
        return id.to_string();
    }
    if id.ends_with("/1") || id.ends_with("/2") {
        return id[0..sz - 2].to_string();
    }
    id.to_string()
}

pub fn open_file<P: AsRef<Path>>(path: P) -> Result<File> {
    File::open(&path).map_err(|e| {
        if e.kind() == io::ErrorKind::NotFound {
            io::Error::new(e.kind(), format!("File not found: {:?}", path.as_ref()))
        } else {
            e
        }
    })
}

pub(crate) fn detect_file_format<P: AsRef<Path>>(path: P) -> io::Result<SeqFormat> {
    // let mut file = open_file(path)?;
    let read1: Box<dyn io::Read + Send> = dyn_reader(path)?;
    let reader = BufReader::new(read1);
    let mut lines = reader.lines();

    if let Some(first_line) = lines.next() {
        let line = first_line?;

        if line.starts_with('>') {
            return Ok(SeqFormat::Fasta);
        } else if line.starts_with('@') {
            let _ = lines.next();
            if let Some(third_line) = lines.next() {
                let line: String = third_line?;
                if line.starts_with('+') {
                    return Ok(SeqFormat::Fastq);
                }
            }
        } else {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Unrecognized fasta(fastq) file format",
            ));
        }
    }

    Err(io::Error::new(
        io::ErrorKind::Other,
        "Unrecognized fasta(fastq) file format",
    ))
}

pub(crate) fn trim_end(buffer: &mut Vec<u8>) {
    while let Some(&b'\n' | &b'\r' | &b'>' | &b'@') = buffer.last() {
        buffer.pop();
    }
}

pub const BUFSIZE: usize = 16 * 1024 * 1024;

pub trait Reader: Send {
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>>;
}

impl Reader for Box<dyn Reader + Send> {
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        (**self).next()
    }
}

#[derive(Debug)]
pub struct PosData {
    /// 外部 taxonomy id
    pub ext_code: u64,
    /// 连续命中次数
    pub count: usize,
}

impl PosData {
    pub fn new(ext_code: u64, count: usize) -> Self {
        Self { ext_code, count }
    }
}

impl fmt::Display for PosData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.ext_code, self.count)
    }
}

/// 序列的命中分布
#[derive(Debug)]
pub struct SpaceDist {
    pub value: Vec<PosData>,
    /// example: (0, 10], 左开右闭
    pub range: (usize, usize),
    pos: usize,
}

impl SpaceDist {
    pub fn new(range: (usize, usize)) -> Self {
        Self {
            value: Vec::new(),
            range,
            pos: range.0,
        }
    }

    fn fill_with_zeros(&mut self, gap: usize) {
        if gap > 0 {
            self.value.push(PosData::new(0, gap));
        }
    }

    pub fn add(&mut self, ext_code: u64, pos: usize) {
        if pos <= self.pos || pos > self.range.1 {
            return; // 早期返回，不做任何处理
        }
        let gap = pos - self.pos - 1;

        if gap > 0 {
            self.fill_with_zeros(gap);
        }

        if let Some(last) = self.value.last_mut() {
            if last.ext_code == ext_code {
                last.count += 1;
            } else {
                self.value.push(PosData::new(ext_code, 1));
            }
        } else {
            self.value.push(PosData::new(ext_code, 1));
        }
        self.pos = pos;
    }

    /// Fills the end of the distribution with zeros if there is remaining space.
    pub fn fill_tail_with_zeros(&mut self) {
        if self.pos < self.range.1 {
            self.fill_with_zeros(self.range.1 - self.pos);
            self.pos = self.range.1;
        }
    }
}

impl fmt::Display for SpaceDist {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, data) in self.value.iter().enumerate() {
            if i > 0 {
                write!(f, " ")?;
            }
            write!(f, "{}", data)?;
        }
        write!(f, "")
    }
}

impl OptionPair<SpaceDist> {
    pub fn add(&mut self, ext_code: u64, pos: usize) {
        match self {
            OptionPair::Single(sd) => sd.add(ext_code, pos),
            OptionPair::Pair(sd1, sd2) => {
                if pos > sd1.range.1 {
                    sd2.add(ext_code, pos)
                } else {
                    sd1.add(ext_code, pos)
                }
            }
        }
    }

    pub fn fill_tail_with_zeros(&mut self) {
        self.apply_mut(|sd| sd.fill_tail_with_zeros());
    }
}
