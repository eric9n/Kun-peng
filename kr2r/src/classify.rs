use crate::compact_hash::{Compact, Row};
use crate::readcounts::TaxonCountersDash;
use crate::taxonomy::Taxonomy;
use std::collections::HashMap;

fn generate_hit_string(
    count: u32,
    rows: &Vec<Row>,
    taxonomy: &Taxonomy,
    value_mask: usize,
    offset: u32,
) -> String {
    let mut result = Vec::new();
    let mut last_pos = 0;

    for row in rows {
        if row.kmer_id < offset || row.kmer_id >= offset + count {
            continue;
        }
        let adjusted_pos = row.kmer_id - offset;

        let value = row.value;
        let key = value.right(value_mask);
        let ext_code = taxonomy.nodes[key as usize].external_id;

        if last_pos == 0 && adjusted_pos > 0 {
            result.push((0, adjusted_pos)); // 在开始处添加0
        } else if adjusted_pos - last_pos > 1 {
            result.push((0, adjusted_pos - last_pos - 1)); // 在两个特定位置之间添加0
        }
        if let Some(last) = result.last_mut() {
            if last.0 == ext_code {
                last.1 += 1;
                last_pos = adjusted_pos;
                continue;
            }
        }

        // 添加当前key的计数
        result.push((ext_code, 1));
        last_pos = adjusted_pos;
    }

    // 填充尾随0
    if last_pos < count - 1 {
        if last_pos == 0 {
            result.push((0, count - last_pos));
        } else {
            result.push((0, count - last_pos - 1));
        }
    }

    result
        .iter()
        .map(|i| format!("{}:{}", i.0, i.1))
        .collect::<Vec<String>>()
        .join(" ")
}

pub fn trim_pair_info(id: &str) -> String {
    let sz = id.len();
    if sz <= 2 {
        return id.to_string();
    }
    if id.ends_with("/1") || id.ends_with("/2") {
        return id[0..sz - 2].to_string();
    }
    id.to_string()
}

// &HashMap<u32, u64>,
pub fn resolve_tree(
    hit_counts: &HashMap<u32, u64>,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    confidence_threshold: f64,
) -> u32 {
    let required_score = (confidence_threshold * total_minimizers as f64).ceil() as u64;

    let mut max_taxon = 0u32;
    let mut max_score = 0;

    for (&taxon, _) in hit_counts {
        let mut score = 0;

        for (&taxon2, &count2) in hit_counts {
            if taxonomy.is_a_ancestor_of_b(taxon2, taxon) {
                score += count2;
            }
        }

        if score > max_score {
            max_score = score;
            max_taxon = taxon;
        } else if score == max_score {
            max_taxon = taxonomy.lca(max_taxon, taxon);
        }
    }

    max_score = *hit_counts.get(&max_taxon).unwrap_or(&0);

    while max_taxon != 0 && max_score < required_score {
        max_score = hit_counts
            .iter()
            .filter(|(&taxon, _)| taxonomy.is_a_ancestor_of_b(max_taxon, taxon))
            .map(|(_, &count)| count)
            .sum();

        if max_score >= required_score {
            break;
        }
        max_taxon = taxonomy.nodes[max_taxon as usize].parent_id as u32;
    }

    max_taxon
}

pub fn add_hitlist_string(
    rows: &Vec<Row>,
    value_mask: usize,
    kmer_count1: u32,
    kmer_count2: Option<u32>,
    taxonomy: &Taxonomy,
) -> String {
    let result1 = generate_hit_string(kmer_count1, &rows, taxonomy, value_mask, 0);
    if let Some(count) = kmer_count2 {
        let result2 = generate_hit_string(count, &rows, taxonomy, value_mask, kmer_count1);
        format!("{} |:| {}", result1, result2)
    } else {
        format!("{}", result1)
    }
}

pub fn count_values(
    rows: &Vec<Row>,
    value_mask: usize,
    kmer_count1: u32,
) -> (HashMap<u32, u64>, TaxonCountersDash, usize) {
    let mut counts = HashMap::new();

    let mut hit_count: usize = 0;

    let mut last_row: Row = Row::new(0, 0, 0);
    let cur_taxon_counts = TaxonCountersDash::new();

    for row in rows {
        let value = row.value;
        let key = value.right(value_mask);
        *counts.entry(key).or_insert(0) += 1;

        // 如果切换到第2条seq,就重新计算
        if last_row.kmer_id < kmer_count1 && row.kmer_id > kmer_count1 {
            last_row = Row::new(0, 0, 0);
        }
        if !(last_row.value == value && row.kmer_id - last_row.kmer_id == 1) {
            cur_taxon_counts
                .entry(key as u64)
                .or_default()
                .add_kmer(value as u64);
            hit_count += 1;
        }

        last_row = *row;
    }

    (counts, cur_taxon_counts, hit_count)
}
