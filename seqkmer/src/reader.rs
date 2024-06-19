use crate::seq::{BaseType, SeqFormat, SeqHeader};
use flate2::read::GzDecoder;
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

pub(crate) fn trim_end(buffer: &mut Vec<u8>) {
    while let Some(&b'\n' | &b'\r' | &b'>' | &b'@') = buffer.last() {
        buffer.pop();
    }
}

pub const BUFSIZE: usize = 16 * 1024 * 1024;

pub type SeqVecType = Vec<BaseType<SeqHeader, Vec<u8>>>;

pub trait Reader: Send {
    fn next(&mut self) -> Result<Option<SeqVecType>>;
}

impl Reader for Box<dyn Reader> {
    fn next(&mut self) -> Result<Option<SeqVecType>> {
        (**self).next()
    }
}

#[derive(Debug)]
pub struct HitGroup<T> {
    /// minimizer data size
    pub marker_size: usize,
    /// hit value vector
    pub rows: Vec<T>,
    /// pair offset
    pub offset: u32,
}

impl<T> HitGroup<T> {
    pub fn new(marker_size: usize, rows: Vec<T>, offset: u32) -> Self {
        Self {
            marker_size,
            rows,
            offset,
        }
    }
}

impl<S, T> BaseType<S, HitGroup<T>> {
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
        if let BaseType::Pair(_, ref first, ref mut second) = self {
            second.offset = first.marker_size as u32;
        }
    }

    pub fn total_marker_size(&self) -> usize {
        match &self {
            BaseType::Single(_, hit) => hit.marker_size,
            BaseType::Pair(_, hit1, hit2) => hit1.marker_size + hit2.marker_size,
        }
    }
}
