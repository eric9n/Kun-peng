#[derive(Debug, Clone)]
pub enum OptionPair<T> {
    Single(T),
    Pair(T, T),
}

impl<T> OptionPair<T> {
    pub fn single(&self) -> Option<&T> {
        match self {
            OptionPair::Single(value) => Some(value),
            _ => None,
        }
    }

    // 它接受一个泛型闭包 F，并返回一个新的 OptionPair<U>
    pub fn map<U, E, F>(&self, mut f: F) -> Result<OptionPair<U>, E>
    where
        F: FnMut(&T) -> Result<U, E>,
    {
        match self {
            OptionPair::Single(t) => f(t).map(OptionPair::Single),
            OptionPair::Pair(t1, t2) => {
                let u1 = f(t1)?;
                let u2 = f(t2)?;
                Ok(OptionPair::Pair(u1, u2))
            }
        }
    }

    // pub fn concat<U, V, F>(&self, init: &mut U, mut f: F) -> V
    // where
    //     F: FnMut(&mut U, &T) -> V,
    // {
    //     match self {
    //         OptionPair::Single(t) => f(init, t),
    //         OptionPair::Pair(t1, t2) => {
    //             f(init, t1);
    //             f(init, t2)
    //         }
    //     }
    // }
    pub fn reduce<U, F>(&self, init: U, mut f: F) -> U
    where
        F: FnMut(U, &T) -> U,
    {
        match self {
            OptionPair::Single(t) => f(init, t),
            OptionPair::Pair(t1, t2) => {
                let result = f(init, t1);
                f(result, t2)
            }
        }
    }

    pub fn reduce_str<F>(&self, sep: &str, mut f: F) -> String
    where
        F: FnMut(&T) -> String,
    {
        self.reduce(String::new(), |acc, t| {
            if acc.is_empty() {
                f(t)
            } else {
                format!("{}{}{}", acc, sep, f(t))
            }
        })
    }

    pub fn apply<U, F>(&self, mut f: F) -> OptionPair<U>
    where
        F: FnMut(&T) -> U,
    {
        match self {
            OptionPair::Single(t) => OptionPair::Single(f(t)),
            OptionPair::Pair(t1, t2) => OptionPair::Pair(f(t1), f(t2)),
        }
    }

    pub fn apply_mut<U, F>(&mut self, mut f: F) -> OptionPair<U>
    where
        F: FnMut(&mut T) -> U,
    {
        match self {
            OptionPair::Single(t) => OptionPair::Single(f(t)),
            OptionPair::Pair(t1, t2) => OptionPair::Pair(f(t1), f(t2)),
        }
    }
}

impl<T: Clone> OptionPair<T> {
    pub fn from_slice(slice: &[T]) -> OptionPair<T> {
        match slice {
            [a, b] => OptionPair::Pair(a.clone(), b.clone()),
            [a] => OptionPair::Single(a.clone()),
            _ => unreachable!(),
        }
    }
}

impl<T> From<(T, Option<T>)> for OptionPair<T> {
    fn from(tuple: (T, Option<T>)) -> Self {
        match tuple {
            (a, Some(b)) => OptionPair::Pair(a, b),
            (a, None) => OptionPair::Single(a),
        }
    }
}
