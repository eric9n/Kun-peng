use crate::utils::OptionPair;

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum SeqFormat {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqHeader {
    pub id: String,
    pub file_index: usize,
    pub reads_index: usize,
    pub format: SeqFormat,
}

#[derive(Debug)]
pub struct Base<T> {
    pub header: SeqHeader,
    pub body: OptionPair<T>,
}

impl<T> Base<T> {
    pub fn new(header: SeqHeader, body: OptionPair<T>) -> Self {
        Self { header, body }
    }

    pub fn map<U, E, F>(&self, mut f: F) -> Result<Base<U>, E>
    where
        F: FnMut(&T) -> Result<U, E>,
    {
        self.body.map(|t| f(&t)).map(|body| Base {
            header: self.header.clone(),
            body,
        })
    }
}
