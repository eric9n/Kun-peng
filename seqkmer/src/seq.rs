#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum SeqFormat {
    Fasta,
    Fastq,
    PairFastq,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaseType<T> {
    Single(T),
    Pair((T, T)),
}

impl<T> BaseType<T> {
    // 泛型方法，根据序列类型执行操作
    pub fn apply<U, F>(&self, mut func: F) -> BaseType<U>
    where
        F: FnMut(&T) -> U,
    {
        match self {
            BaseType::Single(seq) => BaseType::Single(func(seq)),
            BaseType::Pair((seq1, seq2)) => BaseType::Pair((func(seq1), func(seq2))),
        }
    }

    pub fn fold<U, F, V>(&self, init: &mut V, mut func: F) -> BaseType<U>
    where
        F: FnMut(&mut V, &T) -> U,
    {
        match self {
            BaseType::Single(seq) => BaseType::Single(func(init, seq)),
            BaseType::Pair((seq1, seq2)) => {
                let res1 = func(init, seq1);
                let res2 = func(init, seq2);
                BaseType::Pair((res1, res2))
            }
        }
    }

    pub fn modify<F>(&mut self, mut func: F)
    where
        F: FnMut(&mut T),
    {
        match self {
            BaseType::Single(ref mut seq) => func(seq),
            BaseType::Pair((ref mut seq1, ref mut seq2)) => {
                func(seq1);
                func(seq2);
            }
        }
    }
}

impl<U> BaseType<Vec<U>> {
    pub fn len(&self) -> BaseType<usize> {
        self.apply(|seq| seq.len())
    }
}

#[derive(Debug, Clone)]
pub struct Marker {
    pub cap: usize,
    pub minimizer: Vec<u64>,
}

impl Marker {
    pub fn new(cap: usize, minimizer: Vec<u64>) -> Self {
        Self { cap, minimizer }
    }

    pub fn size(&self) -> usize {
        self.minimizer.len()
    }
}

#[derive(Debug, Clone)]
pub struct Sequence {
    pub file_index: usize,
    pub reads_index: usize,
    pub id: String,
    pub seq: BaseType<Vec<u8>>,
    pub format: SeqFormat,
}
