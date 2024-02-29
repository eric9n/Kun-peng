// kraken 2 使用的是murmur_hash3 算法的 fmix64作为 hash
use crate::{canonical_representation, char_to_value, Meros, BITS_PER_CHAR};
use std::collections::VecDeque;

#[inline]
fn to_candidate_lmer(meros: &Meros, lmer: u64) -> u64 {
    let mut canonical_lmer = canonical_representation(lmer, meros.l_mer);
    if meros.spaced_seed_mask > 0 {
        canonical_lmer &= meros.spaced_seed_mask;
    }
    canonical_lmer ^ meros.toggle_mask
}

#[derive(Debug)]
struct MinimizerData {
    pos: usize,
    candidate_lmer: u64,
}

impl MinimizerData {
    fn new(candidate_lmer: u64, pos: usize) -> Self {
        Self {
            candidate_lmer,
            pos,
        }
    }
}

pub struct MinimizerWindow {
    queue: VecDeque<MinimizerData>,
    /// 窗口队列的大小
    capacity: usize,
    /// 队列计数
    count: usize,
}

impl MinimizerWindow {
    fn new(capacity: usize) -> Self {
        Self {
            queue: VecDeque::with_capacity(capacity),
            capacity,
            count: 0,
        }
    }

    #[inline]
    fn next(&mut self, candidate_lmer: u64) -> Option<u64> {
        // 无需比较，直接返回
        if self.capacity == 1 {
            return Some(candidate_lmer);
        }

        let data = MinimizerData::new(candidate_lmer, self.count);

        // 移除队列中所有比当前元素大的元素的索引
        // 因为它们不可能是当前窗口的最小值
        while let Some(m_data) = self.queue.back() {
            if m_data.candidate_lmer > candidate_lmer {
                self.queue.pop_back();
            } else {
                break;
            }
        }
        // 将当前元素的索引添加到队列
        self.queue.push_back(data);
        // 确保队列的第一个元素在当前窗口内
        if self.count < self.capacity {
            self.count += 1;
            return None;
        } else if self
            .queue
            .front()
            .map_or(false, |front| front.pos < self.count - self.capacity)
        {
            self.queue.pop_front();
        }

        self.count += 1;
        self.queue.front().map(|front| front.candidate_lmer)
    }

    fn clear(&mut self) {
        self.count = 0;
        self.queue.clear();
    }
}

struct Cursor {
    pos: usize,
    end: usize,
    inner: Vec<u64>,
    capacity: usize,
    value: u64,
    mask: u64,
    window: MinimizerWindow,
}

impl Cursor {
    fn new(meros: &Meros) -> Self {
        Self {
            pos: 0,
            end: 0,
            inner: Vec::with_capacity(meros.l_mer),
            capacity: meros.l_mer,
            value: 0,
            mask: meros.mask,
            window: MinimizerWindow::new(meros.window_size()),
        }
    }

    /// 每次取一个 lmer 值出来，如果为空，表示一直 seq 已处理完成
    /// 遇到换行符,就跳过.
    #[inline]
    fn slide(&mut self, seq: &[u8]) -> Option<u64> {
        while self.pos < self.end {
            let ch = seq[self.pos];
            let code = if ch == b'\n' || ch == b'\r' {
                self.pos += 1;
                char_to_value(seq[self.pos])
            } else {
                char_to_value(ch)
            };
            // let code = char_to_value(seq[self.pos]);
            self.pos += 1;
            if let Some(c) = code {
                if let Some(lmer) = self.next_lmer(c) {
                    return Some(lmer);
                }
            } else {
                self.clear();
            }
        }
        None
    }

    fn next_lmer(&mut self, item: u64) -> Option<u64> {
        self.value <<= BITS_PER_CHAR;
        self.value |= item;
        if self.inner.len() == self.capacity {
            self.inner.remove(0); // 移除最旧的元素
        }
        self.inner.push(item); // 使用 push 方法
        if self.inner.len() >= self.capacity {
            self.value &= self.mask;
            return Some(self.value);
        }

        None
    }

    #[inline]
    fn next_candidate_lmer(&mut self, item: u64) -> Option<u64> {
        self.window.next(item)
    }

    pub fn has_next(&self) -> bool {
        self.pos < self.end
    }

    // 清除元素
    #[inline]
    fn clear(&mut self) {
        self.inner.clear();
        self.value = 0;
        self.window.clear();
    }
}

pub struct MinimizerScanner {
    meros: Meros,
    // l_mer: usize,
    cursor: Cursor,
    /// 存最近一个最小值
    last_minimizer: u64,
}

impl MinimizerScanner {
    pub fn reset(&mut self) {
        self.cursor.clear();
        self.last_minimizer = std::u64::MAX;
    }

    pub fn new(meros: Meros) -> Self {
        Self {
            meros,
            cursor: Cursor::new(&meros),
            last_minimizer: std::u64::MAX,
        }
    }

    pub fn set_seq_end(&mut self, seq: &[u8]) {
        self.cursor.pos = 0;
        self.cursor.end = seq.len();
    }

    fn get_last_minimizer(&mut self) -> Option<u64> {
        None
        // self.cursor
        //     .window
        //     .get_last_minimizer()
        //     .map(|minimizer| minimizer ^ self.toggle_mask)
    }

    /// 在一个序列上滑动一个光标（可能是为了找到下一个有意义的片段或窗口），
    /// 并对滑动得到的片段进行某种转换或处理。如果光标无法继续滑动（例如到达序列的末尾），则返回 None。
    fn next_window(&mut self, seq: &[u8]) -> Option<u64> {
        self.cursor.slide(seq).and_then(|lmer| {
            let candidate_lmer = to_candidate_lmer(&self.meros, lmer);
            self.cursor.next_candidate_lmer(candidate_lmer)
        })
    }

    /// 去除重复的值
    pub fn next_minimizer(&mut self, seq: &[u8]) -> Option<u64> {
        while self.cursor.has_next() {
            if let Some(minimizer) = self.next_window(&seq) {
                if minimizer != self.last_minimizer {
                    self.last_minimizer = minimizer;
                    return Some(minimizer ^ self.meros.toggle_mask);
                }
            }
        }
        // 检查滑动队列中是否存在值
        let last_minimizer = self.get_last_minimizer();
        // 清空所有的值，等下次换取时，必然等于 None
        self.cursor.clear();
        last_minimizer
    }
}

pub struct KmerIterator<'a> {
    seq: &'a [u8],
    meros: Meros,
    // l_mer: usize,
    cursor: Cursor,
    /// 存最近一个最小值
    last_minimizer: u64,
}

impl<'a> KmerIterator<'a> {
    pub fn new(seq: &'a [u8], meros: Meros) -> Self {
        KmerIterator {
            seq,
            meros,
            cursor: Cursor::new(&meros),
            last_minimizer: std::u64::MAX,
        }
    }

    /// 在一个序列上滑动一个光标（可能是为了找到下一个有意义的片段或窗口），
    /// 并对滑动得到的片段进行某种转换或处理。如果光标无法继续滑动（例如到达序列的末尾），则返回 None。
    fn next_window(&mut self) -> Option<u64> {
        self.cursor.slide(self.seq).and_then(|lmer| {
            let candidate_lmer = to_candidate_lmer(&self.meros, lmer);
            self.cursor.next_candidate_lmer(candidate_lmer)
        })
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        while self.cursor.has_next() {
            if let Some(minimizer) = self.next_window() {
                if minimizer != self.last_minimizer {
                    self.last_minimizer = minimizer;
                    return Some(minimizer ^ self.meros.toggle_mask);
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // 编写测试函数
    #[test]
    fn test_minimizer_scanner() {
        // 在这里编写测试代码
        // 使用 assert_eq!、assert!、assert_ne! 等宏来断言测试条件是否为真
        // 如果条件不为真，测试将失败
        let seq: Vec<u8> = b"ACGATCGACGACG".to_vec();
        let meros = Meros::new(10, 5, None, None, None);
        let mut scanner = MinimizerScanner::new(meros);
        scanner.set_seq_end(&seq);
        let m1 = scanner.next_minimizer(&seq);
        let mm1 = format!("{:016x}", m1.unwrap());
        assert_eq!(mm1, "00000000000002d8");
        let m2 = scanner.next_minimizer(&seq);
        let mm2 = format!("{:016x}", m2.unwrap());
        assert_eq!(mm2, "0000000000000218");
    }

    #[test]
    fn test_minimizer() {
        // 1, 2, 3, 4
        let seq: Vec<u64> = vec![1, 2, 3, 4];
        // 窗口大小 = 2 - 0 + 1
        let mut mini: MinimizerWindow = MinimizerWindow::new(1);
        let mut result = vec![];
        for s in seq {
            if let Some(a) = mini.next(s) {
                result.push(a);
            }
        }
        // if let Some(a) = mini.get_last_minimizer() {
        //     result.push(a);
        // }
        assert_eq!(result, [1, 2, 3, 4]);

        let seq: Vec<u64> = vec![4, 3, 5, 2, 6, 2, 1];
        // 窗口大小 = 2 - 0 + 1
        let mut mini = MinimizerWindow::new(2);
        let mut result = vec![];
        for s in seq {
            if let Some(a) = mini.next(s) {
                result.push(a);
            }
        }
        // if let Some(a) = mini.get_last_minimizer() {
        //     result.push(a);
        // }
        assert_eq!(result, [3, 2, 2, 2, 1]);
    }
}
