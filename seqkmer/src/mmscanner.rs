// kraken 2 使用的是murmur_hash3 算法的 fmix64作为 hash
use crate::seq::{BaseType, Marker};
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
    capacity: usize,
    value: u64,
    mask: u64,
}

impl Cursor {
    fn new(capacity: usize, mask: u64) -> Self {
        Self {
            pos: 0,
            value: 0,
            capacity,
            mask,
        }
    }

    fn next_lmer(&mut self, item: u64) -> Option<u64> {
        self.value = ((self.value << BITS_PER_CHAR) | item) & self.mask;
        // 更新当前位置
        self.pos += 1;
        // 检查是否达到了容量
        if self.pos >= self.capacity {
            return Some(self.value);
        }
        None
    }

    // 清除元素
    #[inline]
    fn clear(&mut self) {
        self.pos = 0;
        self.value = 0;
    }
}

pub struct MinimizerScanner<'a> {
    seq: &'a BaseType<Vec<u8>>,
    meros: Meros,
    cursor: Cursor,
    window: MinimizerWindow,
}

impl<'a> MinimizerScanner<'a> {
    pub fn new(seq: &'a BaseType<Vec<u8>>, meros: Meros) -> Self {
        MinimizerScanner {
            seq,
            meros,
            cursor: Cursor::new(meros.l_mer, meros.mask),
            window: MinimizerWindow::new(meros.window_size()),
        }
    }

    #[inline]
    fn clear(&mut self) {
        self.cursor.clear();
        self.window.clear();
    }

    fn iter_seq(&mut self, seq: &Vec<u8>) -> Marker {
        let minimizer = seq
            .iter()
            .filter_map(|&ch| {
                if ch == b'\n' || ch == b'\r' {
                    None
                } else {
                    match char_to_value(ch) {
                        Some(code) => self.cursor.next_lmer(code).and_then(|lmer| {
                            let candidate_lmer: u64 = to_candidate_lmer(&self.meros, lmer);
                            self.window
                                .next(candidate_lmer)
                                .map(|minimizer| murmur_hash3(minimizer ^ self.meros.toggle_mask))
                        }),
                        None => {
                            self.clear();
                            None
                        }
                    }
                }
            })
            .collect();

        Marker::new(seq.len(), minimizer)
    }

    pub fn iter(&mut self) -> BaseType<Marker> {
        self.seq.apply(|seq| self.iter_seq(seq))
    }
}
