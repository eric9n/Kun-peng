use crate::compact_hash::Compact;
use crate::readcounts::TaxonCounters;
use crate::taxonomy::Taxonomy;
use crate::HitGroup;
use seqkmer::SpaceDist;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Resolves the taxonomic classification based on hit counts and taxonomy.
///
/// This function determines the most likely taxonomic classification for a sequence
/// based on the hit counts for different taxa and the taxonomic hierarchy.
///
/// # Arguments
///
/// * `hit_counts` - A HashMap containing the hit counts for each taxon.
/// * `taxonomy` - The Taxonomy object representing the taxonomic hierarchy.
/// * `required_score` - The minimum score required for a classification to be considered valid.
///
/// # Returns
///
/// Returns the taxon ID of the resolved classification.
pub fn resolve_tree(
    hit_counts: &HashMap<u32, u64>,
    taxonomy: &Taxonomy,
    required_score: u64,
) -> u32 {
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

/// Processes hit statistics for a group of hits.
///
/// This function calculates various statistics for a group of hits, including
/// updating hit counts, taxon counters, and computing space distribution.
///
/// # Arguments
///
/// * `hits` - The HitGroup to process.
/// * `counts` - A mutable reference to a HashMap to store hit counts.
/// * `value_mask` - A mask used for processing hit values.
/// * `taxonomy` - The Taxonomy object representing the taxonomic hierarchy.
/// * `cur_taxon_counts` - A mutable reference to TaxonCounters to update.
///
/// # Returns
///
/// Returns a String representing the space distribution of the hits.
fn stat_hits<'a>(
    hits: &HitGroup,
    counts: &mut HashMap<u32, u64>,
    value_mask: usize,
    taxonomy: &Taxonomy,
    cur_taxon_counts: &mut TaxonCounters,
) -> String {
    let mut space_dist = hits.range.apply(|range| SpaceDist::new(*range));
    for row in &hits.rows {
        let value = row.value;
        let key = value.right(value_mask);

        *counts.entry(key).or_insert(0) += 1;

        cur_taxon_counts
            .entry(key as u64)
            .or_default()
            .add_kmer(value as u64);

        let ext_code = taxonomy.nodes[key as usize].external_id;
        let pos = row.kmer_id as usize;
        space_dist.add(ext_code, pos);
    }

    space_dist.fill_tail_with_zeros();
    space_dist.reduce_str(" |:| ", |str| str.to_string())
}

/// Processes a hit group to determine classification and gather statistics.
///
/// This function takes a hit group, processes it to determine the taxonomic
/// classification, and collects various statistics about the hits.
///
/// # Arguments
///
/// * `hits` - The HitGroup to process.
/// * `taxonomy` - The Taxonomy object representing the taxonomic hierarchy.
/// * `classify_counter` - An atomic counter for tracking classifications.
/// * `required_score` - The minimum score required for a classification to be considered valid.
/// * `minimum_hit_groups` - The minimum number of hit groups required for a valid classification.
/// * `value_mask` - A mask used for processing hit values.
///
/// # Returns
///
/// Returns a tuple containing:
/// 1. A String indicating the classification result ("C" for classified, "U" for unclassified).
/// 2. The external ID of the classified taxon.
/// 3. A String representing the hit statistics.
/// 4. The updated TaxonCounters.
pub fn process_hitgroup(
    hits: &HitGroup,
    taxonomy: &Taxonomy,
    classify_counter: &AtomicUsize,
    required_score: u64,
    minimum_hit_groups: usize,
    value_mask: usize,
) -> (String, u64, String, TaxonCounters) {
    let mut cur_taxon_counts = TaxonCounters::new();
    let mut counts = HashMap::new();
    let hit_groups = hits.capacity();
    let hit_string = stat_hits(
        hits,
        &mut counts,
        value_mask,
        taxonomy,
        &mut cur_taxon_counts,
    );

    let mut call = resolve_tree(&counts, taxonomy, required_score);
    if call > 0 && hit_groups < minimum_hit_groups {
        call = 0;
    };

    let ext_call = taxonomy.nodes[call as usize].external_id;
    let clasify = if call > 0 {
        classify_counter.fetch_add(1, Ordering::SeqCst);
        cur_taxon_counts
            .entry(call as u64)
            .or_default()
            .increment_read_count();

        "C"
    } else {
        "U"
    };

    (clasify.to_owned(), ext_call, hit_string, cur_taxon_counts)
}
