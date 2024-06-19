// #[derive(Debug, Clone)]
// pub enum OptionPair<S, T> {
//     Single(S, T),
//     Pair(S, T, T),
// }

// impl<S, T> OptionPair<S, T> {
//     // 它接受一个泛型闭包 F，并返回一个新的 OptionPair<U>
//     pub fn map<U, E, F>(self, mut f: F) -> Result<OptionPair<S, U>, E>
//     where
//         F: FnMut(S, T) -> Result<U, E>,
//     {
//         match self {
//             OptionPair::Single(s, t) => f(s, t).map(|u| OptionPair::Single(s.clone(), u)),
//             OptionPair::Pair(s, t1, t2) => {
//                 let u1 = f(s, t1)?;
//                 let u2 = f(s, t2)?;
//                 Ok(OptionPair::Pair(s, u1, u2))
//             }
//         }
//     }
// }

// impl<S, T: Clone> OptionPair<S, T> {
//     pub fn from_slice(s: S, slice: &[T]) -> OptionPair<S, T> {
//         match slice {
//             [a, b] => OptionPair::Pair(s, a.clone(), b.clone()),
//             [a] => OptionPair::Single(s, a.clone()),
//             _ => unreachable!(),
//         }
//     }
// }

// impl<S, T> From<(S, T, Option<T>)> for OptionPair<S, T> {
//     fn from(tuple: (S, T, Option<T>)) -> Self {
//         match tuple {
//             (s, a, Some(b)) => OptionPair::Pair(s, a, b),
//             (s, a, None) => OptionPair::Single(s, a),
//         }
//     }
// }
