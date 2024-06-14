#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum SeqFormat {
    FASTA,
    FASTQ,
}

#[derive(Debug, Clone)]
pub struct Sequence {
    pub file_index: u64,
    pub reads_index: u64,
    pub id: String,
    pub seq: Vec<u8>,
    pub format: SeqFormat,
}
