use std::collections::HashSet;

// -------
// Helper methods for bit operations

const BITMASK: u32 = ((1u32 << (7 - 1)) - 1) << 1;

/// Extract bits from a u32 value using LSB 0 numbering from 7 to 1, including 1,
/// and then shift the extracted bits to the right.
#[inline]
pub fn extract_bits(value: u32) -> u8 {
    // Use the constant bitmask for extracting bits
    let result = value & BITMASK;

    // Shift resulting bits to the right by 1 position
    (result >> 1) as u8
}

/// Extracts the high `hi` bits from a `u64` value.
///
/// # Arguments
/// * `bits` - The `u64` value from which to extract bits.
/// * `hi` - The number of high bits to extract.
///
/// # Returns
/// Returns the extracted high bits as a `u64`.
#[inline]
pub fn extract_high_bits_u64(bits: u64, hi: u8) -> u64 {
    bits >> (64u8 - hi)
}

#[inline]
pub fn extract_high_bits_u32(bits: u32, hi: u8) -> u32 {
    bits >> (32u8 - hi)
}

/// Calculate the rank of a hash value for HyperLogLog.
/// Rank is defined as the number of leading zeros in the hash value plus one,
/// but only considering the top `32 - p` bits.
///
/// `hash_value`: The 32-bit hash value.
/// `p`: The precision parameter of HyperLogLog.
#[inline]
pub fn get_rank_u32(hash_value: u32, p: u8) -> u8 {
    // Shift p bits off the hash value, keeping only the top 32-p bits
    let rank_bits = hash_value << p;

    // Count the number of leading zeros in the remaining bits and add one
    let rank_val = rank_bits.leading_zeros() as u8 + 1;

    // Sanity check to ensure the rank value is within expected bounds
    assert!(rank_val < 32 - p as u8 + 1);
    rank_val
}

#[inline]
pub fn get_rank_u64(hash_value: u64, p: u8) -> u8 {
    // 将哈希值左移 p 位，移除掉低位的 p 位
    let rank_bits = hash_value << p;

    // 计算剩余部分的前导零的数量
    let rank_val = rank_bits.leading_zeros() as u8 + 1;

    // 确保 rank_val 在有效范围内
    assert!(rank_val <= 64 - p + 1);
    rank_val
}

#[inline]
pub fn get_index_u64(hash_value: u64, p: u8) -> u32 {
    (hash_value >> (64 - p)) as u32
}

#[inline]
pub fn get_index_u32(hash_value: u32, p: u8) -> u32 {
    hash_value >> (32 - p)
}

/// Encodes the 64-bit hash code as a 32-bit integer for use in the sparse representation.
///
/// # Arguments
/// * `hash_value` - The 64-bit hash value.
/// * `p_prime` - The precision used in the sparse representation.
/// * `p` - The original precision.
///
/// # Returns
/// Returns the encoded 32-bit hash value.
pub fn encode_hash_in_32bit(hash_value: u64, p_prime: u8, p: u8) -> u32 {
    // 提取前 p_prime 位作为索引
    let idx = (extract_high_bits_u64(hash_value, p_prime) << (32 - p_prime)) as u32;

    // 检查 p 到 p_prime 之间的位是否都是 0
    if idx << p == 0 {
        // 计算额外的排名
        let additional_rank = get_rank_u64(hash_value, p_prime) - (p_prime - p);
        idx | ((additional_rank as u32) << 1) | 1
    } else {
        // 返回 idx，最后一位是 0
        assert!(idx & 1 == 0);
        idx
    }
}

/// Calculate the encoded rank of a hash value for HyperLogLog.
///
/// # Arguments
/// * `encoded_hash_value` - A 32-bit encoded hash value.
/// * `p_prime` - The higher precision parameter.
/// * `p` - The original precision parameter.
#[inline]
pub fn get_encoded_rank(encoded_hash_value: u32, p_prime: u8, p: u8) -> u8 {
    // Check if the least significant bit is 1
    if encoded_hash_value & 1 == 1 {
        // If yes: the hash was stored with higher precision, bits p to pPrime were 0
        let additional_rank = p_prime - p;
        let extracted_bits = extract_bits(encoded_hash_value);
        additional_rank + extracted_bits
    } else {
        // Otherwise, calculate the rank at the original precision
        get_rank_u32(encoded_hash_value, p)
    }
}

/// Ertl - calculation of sigma correction for 0-registers in M.
/// x is the proportion of 0-registers, thus x ∈ [0, 1].
/// sigma := x + sum[from k=1 to Inf] ( x^(2^k) * 2^(k-1) )
///
/// # Tolerance
/// The `tolerance` value is used to determine when the iterative calculation
/// has stabilized. It's a balance between accuracy and preventing infinite loops
/// due to floating-point precision limitations. Typical values range from 1e-12
/// to 1e-15, depending on the required precision of the application.
#[inline]
pub fn sigma(x: f64) -> f64 {
    assert!(x >= 0.0 && x <= 1.0);
    if x == 1.0 {
        return f64::INFINITY;
    }

    let mut prev_sigma_x = 0.0;
    let mut sigma_x = x;
    let mut x_squared = x * x; // x^(2^k) for k=1
    let mut y = 1.0;

    // 定义一个小的容忍度，用于浮点数比较
    let tolerance = 1e-12;

    loop {
        prev_sigma_x = sigma_x;
        sigma_x += x_squared * y;
        x_squared *= x_squared; // x^(2^k) for next k
        y += y; // 2^(k-1) for next k

        // 检查变化量是否低于容忍度
        if (sigma_x - prev_sigma_x).abs() < tolerance {
            break;
        }
    }

    sigma_x
}

/// Ertl - calculation of tau correction for values higher than q in M.
/// `x` is the proportion of registers with a value below q in M, thus x ∈ [0,1].
/// tau := 1/3 * (1 - x - sum[from k=1 to Inf](1 - x^(2^(-k))^2 * 2^(-k) ))
#[inline]
pub fn tau(x: f64) -> f64 {
    assert!(x >= 0.0 && x <= 1.0);
    if x == 1.0 {
        return f64::INFINITY;
    }

    let mut prev_tau_x = 0.0;
    let mut tau_x = 1.0 - x;
    let mut x_pow = x;
    let mut y = 1.0;

    // 定义一个小的容忍度，用于浮点数比较
    let tolerance = 1e-12;

    loop {
        prev_tau_x = tau_x;
        x_pow *= x_pow; // x^(2^(-k))^2 for next k
        y *= 0.5; // 2^(-k) for next k
        let term = (1.0 - x_pow) * y;
        tau_x -= term;

        // 检查变化量是否低于容忍度
        if (tau_x - prev_tau_x).abs() < tolerance {
            break;
        }
    }

    tau_x / 3.0
}

#[inline]
pub fn register_histogram(m: &[u8], q: u8) -> Vec<f64> {
    let mut c = vec![0f64; (q + 2) as usize];
    for &val in m {
        c[val as usize] += 1f64;
    }

    assert_eq!(c.iter().sum::<f64>(), m.len() as f64);
    c
}

/// Calculate the histogram of ranks from a sparse list for HyperLogLog.
///
/// # Arguments
/// * `sparse_list` - A set of encoded hash values.
/// * `p_prime` - Higher precision parameter.
/// * `p` - Original precision parameter.
/// * `q` - The number of bits to consider after the original precision.
#[inline]
pub fn sparse_register_histogram(
    sparse_list: &HashSet<u32>,
    p_prime: u8,
    p: u8,
    q: u8,
) -> Vec<f64> {
    let mut c = vec![0f64; (q + 2) as usize];
    let mut m = 1usize << p_prime;

    for &encoded_hash_value in sparse_list {
        let rank_val = get_encoded_rank(encoded_hash_value, p_prime, p) as usize;
        c[rank_val] += 1f64;
        m -= 1;
    }
    c[0] = m as f64;
    c
}

/// Improved cardinality estimator of Ertl, 2017 (arXiv, section 4).
/// Based on the underlying distribution, the estimator employs correction
/// factors for zero and 'over-subscribed' registers. It does not depend on
/// empirically defined bias correction values or a switch between linear
/// counting and loglog estimation.
///
/// # Formula
/// ```
///                                 alpha_inf * m^2
/// --------------------------------------------------------------------------------------
/// ( m * sigma(C_0/m) + sum[from k=1 to q] C_k * 2^(-k) + m * tau(1-C_(q+1)/m) * 2^(-q)
/// ```
///
#[inline]
pub fn ertl_cardinality(alpha_inf: f64, m: f64, c: &[f64], q: usize) -> f64 {
    let mut est_denominator = m * sigma(c[0] / m);

    for k in 1..=q {
        est_denominator += c[k] * 2f64.powi(-(k as i32));
    }

    est_denominator += m * tau(1.0 - c[q + 1] / m) * 2f64.powi(-(q as i32));

    let m_sq_alpha_inf = alpha_inf * m * m;
    m_sq_alpha_inf / est_denominator
}

pub fn alpha_inf(p: u8) -> f64 {
    match 1usize << p {
        16 => 0.673,
        32 => 0.697,
        64 => 0.709,
        m => 0.7213 / (1.0 + 1.079 / m as f64),
    }
}
