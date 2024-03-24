// kraken 2 使用的是murmur_hash3 算法的 fmix64作为 hash
use crate::{
    canonical_representation, char_to_value, fmix64 as murmur_hash3, Meros, BITS_PER_CHAR,
};
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
    pub pos: usize,
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
    queue_pos: usize,
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
            queue_pos: 0,
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
        let mut changed = false;

        if (self.queue.is_empty() && self.count >= self.capacity) || self.count == self.capacity {
            changed = true
        }
        // 将当前元素的索引添加到队列
        self.queue.push_back(data);

        while !self.queue.is_empty()
            && self.queue.front().map_or(false, |front| {
                self.count >= self.capacity && front.pos < self.count - self.capacity
            })
        {
            self.queue.pop_front();
            changed = true;
        }

        self.count += 1;
        if changed {
            self.queue.front().map(|front| front.candidate_lmer)
        } else {
            None
        }
    }

    fn clear(&mut self) {
        self.count = 0;
        self.queue_pos = 0;
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
    fn new(meros: &Meros, size: usize) -> Self {
        Self {
            pos: 0,
            end: size,
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

pub struct MinimizerScanner<'a> {
    seq: &'a [u8],
    meros: Meros,
    // l_mer: usize,
    cursor: Cursor,
    // 存最近一个最小值
    // last_minimizer: u64,
}

impl<'a> MinimizerScanner<'a> {
    pub fn new(seq: &'a [u8], meros: Meros) -> Self {
        let size: usize = seq.len();
        MinimizerScanner {
            seq,
            meros,
            cursor: Cursor::new(&meros, size),
            // last_minimizer: std::u64::MAX,
        }
    }

    /// 在一个序列上滑动一个光标（可能是为了找到下一个有意义的片段或窗口），
    /// 并对滑动得到的片段进行某种转换或处理。如果光标无法继续滑动（例如到达序列的末尾），则返回 None。
    fn next_window(&mut self) -> Option<u64> {
        self.cursor.slide(&self.seq).and_then(|lmer| {
            let candidate_lmer: u64 = to_candidate_lmer(&self.meros, lmer);
            self.cursor.next_candidate_lmer(candidate_lmer)
        })
    }
}

impl<'a> Default for MinimizerScanner<'a> {
    fn default() -> Self {
        let meros = Meros::default();
        let seq: &[u8] = &[];
        MinimizerScanner::new(seq, meros)
    }
}

impl<'a> Iterator for MinimizerScanner<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        while self.cursor.has_next() {
            if let Some(minimizer) = self.next_window() {
                return Some(murmur_hash3(minimizer ^ self.meros.toggle_mask));
                // if minimizer != self.last_minimizer {
                //     self.last_minimizer = minimizer;
                //     return Some(murmur_hash3(minimizer ^ self.meros.toggle_mask));
                // }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
