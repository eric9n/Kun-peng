use crate::taxonomy::Taxonomy;
use std::collections::HashMap;

const TAXID_MAX: u64 = u64::MAX;
const MATE_PAIR_BORDER_TAXON: u64 = TAXID_MAX;
const READING_FRAME_BORDER_TAXON: u64 = TAXID_MAX - 1;
const AMBIGUOUS_SPAN_TAXON: u64 = TAXID_MAX - 2;

fn trim_pair_info(id: &str) -> String {
    let sz = id.len();
    if sz <= 2 {
        return id.to_string();
    }
    if id.ends_with("/1") || id.ends_with("/2") {
        return id[0..sz - 2].to_string();
    }
    id.to_string()
}

fn resolve_tree(
    hit_counts: &HashMap<u64, u32>,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    confidence_threshold: f64,
) -> u64 {
    let mut max_taxon = 0;
    let mut max_score = 0;
    let required_score = (confidence_threshold * total_minimizers as f64).ceil() as u32;

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
            max_taxon = taxonomy.lowest_common_ancestor(max_taxon, taxon);
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
            max_taxon = taxonomy.nodes[max_taxon as usize].parent_id;
        }
    }

    max_taxon
}

fn add_hitlist_string(taxa: &[u64], taxonomy: &Taxonomy) -> String {
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
