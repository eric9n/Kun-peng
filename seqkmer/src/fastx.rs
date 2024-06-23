use crate::fasta::FastaReader;
use crate::fastq::FastqReader;
use crate::reader::{detect_file_format, Reader};
use crate::seq::{Base, SeqFormat};
use crate::utils::OptionPair;
use std::io::Result;
use std::path::Path;

pub struct FastxReader<R: Reader> {
    inner: R,
}

impl<R: Reader> FastxReader<R> {
    pub fn new(inner: R) -> Self {
        Self { inner }
    }
}

impl<R: Reader> Reader for FastxReader<R> {
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        self.inner.next()
    }
}
impl FastxReader<Box<dyn Reader + Send>> {
    pub fn from_paths<P: AsRef<Path>>(
        paths: OptionPair<P>,
        file_index: usize,
        quality_score: i32,
    ) -> Result<Self> {
        let file_format = paths.map(|path: &P| detect_file_format(path));

        match file_format? {
            OptionPair::Single(SeqFormat::Fasta) => {
                let reader = FastaReader::from_path(paths.single().unwrap().as_ref(), file_index)?;
                Ok(Self::new(Box::new(reader) as Box<dyn Reader + Send>))
            }
            OptionPair::Single(SeqFormat::Fastq)
            | OptionPair::Pair(SeqFormat::Fastq, SeqFormat::Fastq) => {
                let reader = FastqReader::from_path(paths, file_index, quality_score)?;
                Ok(Self::new(Box::new(reader) as Box<dyn Reader + Send>))
            }
            _ => panic!("Unsupported file format combination"),
        }
    }
}
