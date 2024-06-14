use crate::seq::Sequence;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, Read, Result, Seek};
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

pub fn trim_end(buffer: &mut Vec<u8>) {
    while let Some(&b'\n' | &b'\r') = buffer.last() {
        buffer.pop();
    }
}

pub const BUFSIZE: usize = 8 * 1024 * 1024;

pub trait Reader<R: Read + Send>: Send {
    fn next(&mut self) -> Result<Option<Sequence>>;
}
