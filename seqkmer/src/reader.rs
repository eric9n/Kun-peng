use crate::seq::{BaseType, Marker, SeqFormat, Sequence};
use crate::{mmscanner::MinimizerScanner, Meros};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Result, Seek};
use std::path::Path;

pub fn dyn_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn Read + Send>> {
    let mut file = open_file(path)?;
    if is_gzipped(&mut file)? {
        let decoder = GzDecoder::new(file);
        Ok(Box::new(decoder))
    } else {
        Ok(Box::new(file))
    }
}

pub fn is_gzipped(file: &mut File) -> Result<bool> {
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;
    file.rewind()?; // 重置文件指针到开头
    Ok(buffer == [0x1F, 0x8B])
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

pub fn detect_file_format<P: AsRef<Path>>(path: P) -> io::Result<SeqFormat> {
    let mut file = open_file(path)?;
    let read1: Box<dyn io::Read + Send> = if is_gzipped(&mut file)? {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };

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

pub fn trim_end(buffer: &mut Vec<u8>) {
    while let Some(&b'\n' | &b'\r' | &b'>' | &b'@') = buffer.last() {
        buffer.pop();
    }
}

pub const BUFSIZE: usize = 16 * 1024 * 1024;

pub trait Reader: Send {
    fn next(&mut self) -> Result<Option<Vec<Sequence>>>;
}

impl Reader for Box<dyn Reader> {
    fn next(&mut self) -> Result<Option<Vec<Sequence>>> {
        (**self).next()
    }
}

#[derive(Debug, Clone)]
pub struct SeqMer {
    pub id: String,
    pub file_index: usize,
    pub reads_index: usize,
    pub marker: BaseType<Marker>,
}

impl SeqMer {
    pub fn from_seq(seq: &Sequence, meros: Meros) -> Self {
        let mut ms = MinimizerScanner::new(&seq.seq, meros);
        let marker = ms.iter();
        Self {
            marker,
            id: seq.id.clone(),
            file_index: seq.file_index,
            reads_index: seq.reads_index,
        }
    }

    pub fn cap_str(&self) -> BaseType<String> {
        self.marker.apply(|marker| marker.cap.to_string())
    }

    pub fn total_size(&self) -> usize {
        match &self.marker {
            BaseType::Single(marker1) => marker1.size(),
            BaseType::Pair((marker1, marker2)) => marker1.size() + marker2.size(),
        }
    }

    pub fn fmt_cap(&self) -> String {
        match &self.marker {
            BaseType::Single(marker1) => marker1.cap.to_string(),
            BaseType::Pair((marker1, marker2)) => format!("{}|{}", marker1.cap, marker2.cap),
        }
    }

    pub fn fmt_size(&self) -> String {
        match &self.marker {
            BaseType::Single(marker1) => marker1.size().to_string(),
            BaseType::Pair((marker1, marker2)) => format!("{}|{}", marker1.size(), marker2.size()),
        }
    }
}

#[derive(Debug)]
pub struct HitGroup<T> {
    /// minimizer capacity
    pub cap: usize,
    /// hit value vector
    pub rows: Vec<T>,
    /// pair offset
    pub offset: u32,
}

impl<T> HitGroup<T> {
    pub fn new(cap: usize, rows: Vec<T>, offset: u32) -> Self {
        Self { cap, rows, offset }
    }
}

impl<T> BaseType<HitGroup<T>> {
    /// Synchronizes the offset of the second element of a `Pair` to the `cap` of the first element.
    /// This alignment is only necessary when the `rows` property of the `HitGroup` is in an
    /// increasing order. If `rows` is not increasing, aligning the offset based on `cap` may not
    /// be appropriate or required.
    ///
    /// # Example
    ///
    /// ```
    /// let mut hit_group1 = HitGroup::new(10, vec![1, 2, 3], 0); // Increasing `rows`
    /// let mut hit_group2 = HitGroup::new(20, vec![4, 5, 6], 0);
    ///
    /// let mut pair = BaseType::Pair((hit_group1, hit_group2));
    /// pair.align_offset();
    /// ```
    pub fn align_offset(&mut self) {
        if let BaseType::Pair((ref first, ref mut second)) = self {
            second.offset = first.cap as u32;
        }
    }
}
