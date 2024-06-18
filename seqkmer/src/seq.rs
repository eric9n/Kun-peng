#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum SeqFormat {
    Fasta,
    Fastq,
    PairFastq,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqHeader {
    pub id: String,
    pub file_index: usize,
    pub reads_index: usize,
    pub format: SeqFormat,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseType<S, T> {
    Single(S, T),
    Pair(S, T, T),
}

impl<S, T> BaseType<S, T> {
    // 泛型方法，根据序列类型执行操作
    pub fn apply<'a, U, F>(&'a self, mut func: F) -> BaseType<(), U>
    where
        F: FnMut(&'a S, &'a T) -> U,
    {
        match self {
            BaseType::Single(prop, seq) => BaseType::Single((), func(prop, seq)),
            BaseType::Pair(prop, seq1, seq2) => {
                BaseType::Pair((), func(prop, seq1), func(prop, seq2))
            }
        }
    }

    pub fn apply_mut<'a, U, F>(&'a mut self, mut func: F) -> BaseType<(), U>
    where
        F: FnMut(&'a S, &'a mut T) -> U,
    {
        match self {
            BaseType::Single(prop, seq) => BaseType::Single((), func(prop, seq)),
            BaseType::Pair(prop, seq1, seq2) => {
                BaseType::Pair((), func(prop, seq1), func(prop, seq2))
            }
        }
    }

    pub fn get_s(&self) -> &S {
        match self {
            BaseType::Single(prop, _) => prop,
            BaseType::Pair(prop, _, _) => prop,
        }
    }

    pub fn transform<'a, U, F, V>(&mut self, mut func: F) -> BaseType<V, U>
    where
        F: for<'b> FnMut(&S, &'b mut T) -> (V, U),
    {
        match self {
            BaseType::Single(prop, seq) => {
                let res1 = func(prop, seq);
                BaseType::Single(res1.0, res1.1)
            }
            BaseType::Pair(prop, seq1, seq2) => {
                let res1 = func(prop, seq1);
                let res2 = func(prop, seq2);
                BaseType::Pair(res1.0, res1.1, res2.1)
            }
        }
    }

    pub fn fold<U, F, V>(&mut self, init: &mut V, mut func: F) -> BaseType<(), U>
    where
        F: FnMut(&mut V, &S, &mut T) -> U,
    {
        match self {
            BaseType::Single(prop, seq) => BaseType::Single((), func(init, prop, seq)),
            BaseType::Pair(prop, ref mut seq1, ref mut seq2) => {
                let res1 = func(init, prop, seq1);
                let res2 = func(init, prop, seq2);
                BaseType::Pair((), res1, res2)
            }
        }
    }

    pub fn modify<F>(&mut self, mut func: F)
    where
        F: FnMut(&S, &mut T),
    {
        match self {
            BaseType::Single(prop, ref mut seq) => func(prop, seq),
            BaseType::Pair(prop, ref mut seq1, ref mut seq2) => {
                func(prop, seq1);
                func(prop, seq2);
            }
        }
    }
}

impl<S, U> BaseType<S, Vec<U>>
where
    S: Copy,
{
    pub fn len(&self) -> BaseType<(), usize> {
        self.apply(|_, seq| seq.len())
    }
}
