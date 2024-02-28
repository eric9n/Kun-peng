use crate::compact_hash::CompactHashTable;
use crate::mmscanner::MinimizerScanner;
use crate::taxonomy::Taxonomy;
use crate::Meros;
use crate::TaxonCounts;
use seq_io::fastq::{Record, RefRecord};
// kraken 2 使用的是murmur_hash3 算法的 fmix64作为 hash
use crate::fmix64 as murmur_hash3;
use std::collections::HashMap;

pub const TAXID_MAX: u32 = u32::MAX - 1;
pub const MATE_PAIR_BORDER_TAXON: u32 = TAXID_MAX;
pub const READING_FRAME_BORDER_TAXON: u32 = TAXID_MAX - 1;
pub const AMBIGUOUS_SPAN_TAXON: u32 = TAXID_MAX - 2;

fn mask_low_quality_bases<'a>(ref_record: &'a RefRecord, minimum_quality_score: i32) -> Vec<u8> {
    let seq = ref_record.seq();
    let qual = ref_record.qual();

    if minimum_quality_score <= 0 {
        return seq.to_vec();
    }
    // // 确保 minimum_quality_score 在 0 到 255 范围内
    // let min_qual = if minimum_quality_score < 0 {
    //     0
    // } else if minimum_quality_score > 255 {
    //     255
    // } else {
    //     minimum_quality_score as u8
    // };

    seq.iter()
        .zip(qual.iter())
        .map(|(&base, &qscore)| {
            if (qscore as i32 - '!' as i32) < minimum_quality_score as i32 {
                b'x'
            } else {
                base
            }
        })
        .collect()
}

pub fn classify_seq<'a>(
    taxonomy: &Taxonomy,
    cht: &CompactHashTable<u32>,
    scanner: &mut MinimizerScanner,
    record_list: &'a Vec<RefRecord>,
    minimum_quality_score: i32,
    meros: Meros,
    confidence_threshold: f64,
    minimum_hit_groups: i32,
    dna_id: String,
) {
    let mut hit_counts = TaxonCounts::new();
    let mut taxa = Vec::<u32>::new();
    let mut minimizer_hit_groups = 0;

    for record in record_list {
        let mut last_minimizer = u64::MAX;
        let mut last_taxon = TAXID_MAX;
        let seq = mask_low_quality_bases(&record, minimum_quality_score);
        scanner.set_seq_end(&seq);
        while let Some(minimizer) = scanner.next_minimizer(&seq) {
            let taxon = if last_minimizer != minimizer {
                let hashed = murmur_hash3(minimizer);
                let mut taxon = 0;
                if meros
                    .min_clear_hash_value
                    .map_or(true, |min_hash| hashed >= min_hash)
                {
                    taxon = cht.get(hashed);
                    last_minimizer = minimizer;
                    last_taxon = taxon;
                    if taxon > 0 {
                        minimizer_hit_groups += 1;
                    }
                }
                taxon
            } else {
                last_taxon
            };
            if taxon > 0 {
                *hit_counts.entry(taxon).or_insert(0) += 1;
            }
            taxa.push(taxon);
        }
        taxa.push(MATE_PAIR_BORDER_TAXON);

        scanner.reset();
    }
    let total_kmers = if record_list.len() > 1 {
        taxa.len() - 2
    } else {
        taxa.len() - 1
    };
    let mut call = resolve_tree(&hit_counts, taxonomy, total_kmers, confidence_threshold);
    if call > 0 && minimizer_hit_groups < minimum_hit_groups {
        call = 0;
    }

    let ext_call = taxonomy.nodes[call as usize].external_id;
    let classify = if call > 0 { "C" } else { "U" };
    println!("{}\t{}\t{}", classify, dna_id, ext_call);
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
        } else {
            max_taxon = taxonomy.nodes[max_taxon as usize].parent_id as u32;
        }
    }

    max_taxon
}

pub fn resolve_tree_optimized(
    hit_counts: &HashMap<u64, u64>,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    confidence_threshold: f64,
) -> u64 {
    let required_score = (confidence_threshold * total_minimizers as f64).ceil() as u64;
    let mut score_cache: HashMap<u64, u64> = HashMap::new();

    // 为每个taxon及其所有祖先累加得分
    for (&taxon, &count) in hit_counts {
        let mut ancestors = vec![taxon];
        let mut current_taxon = taxon;
        // 循环遍历当前taxon的所有祖先
        while current_taxon != 0 {
            let parent_id = taxonomy.nodes[current_taxon as usize].parent_id;
            if parent_id != 0 {
                ancestors.push(parent_id);
                current_taxon = parent_id;
            } else {
                break; // 如果父节点是0，则结束循环
            }
        }

        // 累加得分到score_cache
        for ancestor in ancestors {
            *score_cache.entry(ancestor).or_insert(0) += count;
        }
    }

    // 找到得分最高的taxon，其得分满足所需阈值
    score_cache
        .into_iter()
        .filter(|&(_, score)| score >= required_score)
        .max_by_key(|&(_, score)| score)
        .map(|(taxon, _)| taxon)
        .unwrap_or(0)
}

pub fn add_hitlist_string(taxa: &[u32], taxonomy: &Taxonomy) -> String {
    let mut result = String::new();
    let mut last_code = taxa[0];
    let mut code_count = 1;

    for &code in &taxa[1..] {
        if code == last_code {
            code_count += 1;
        } else {
            match last_code {
                MATE_PAIR_BORDER_TAXON => result.push_str("|:| "),
                READING_FRAME_BORDER_TAXON => result.push_str("-:- "),
                AMBIGUOUS_SPAN_TAXON => result.push_str(&format!("A:{} ", code_count)),
                _ => {
                    let ext_code = taxonomy.nodes[last_code as usize].external_id;
                    result.push_str(&format!("{}:{} ", ext_code, code_count));
                }
            }
            code_count = 1;
            last_code = code;
        }
    }
    match last_code {
        MATE_PAIR_BORDER_TAXON => result.push_str("|:|"),
        READING_FRAME_BORDER_TAXON => result.push_str("-:-"),
        AMBIGUOUS_SPAN_TAXON => result.push_str(&format!("A:{}", code_count)),
        _ => {
            let ext_code = taxonomy.nodes[last_code as usize].external_id;
            result.push_str(&format!("{}:{}", ext_code, code_count));
        }
    }

    result
}
