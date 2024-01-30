use crate::readcounts::{ReadCounter, TaxonCounters};
use crate::taxonomy::Taxonomy;
use std::collections::HashMap;

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

pub fn get_clade_counts(taxonomy: &Taxonomy, call_counts: &HashMap<u64, u64>) -> HashMap<u64, u64> {
    let mut clade_counts = HashMap::new();

    for (&taxid, &count) in call_counts {
        let mut current_taxid = taxid;
        while current_taxid != 0 {
            *clade_counts.entry(current_taxid).or_insert(0) += count;
            current_taxid = taxonomy.nodes[current_taxid as usize].parent_id;
        }
    }

    clade_counts
}

pub fn get_clade_counters(taxonomy: &Taxonomy, call_counters: &TaxonCounters) -> TaxonCounters {
    let mut clade_counters = TaxonCounters::new();

    for (&taxid, counter) in call_counters {
        let mut current_taxid = taxid;
        while current_taxid != 0 {
            let _ = clade_counters
                .entry(current_taxid)
                .or_insert_with(ReadCounter::default)
                .merge(counter);
            current_taxid = taxonomy.nodes[current_taxid as usize].parent_id;
        }
    }

    clade_counters
}

fn extract_string_from_offset(data: &[u8], offset: usize) -> &str {
    let end = data[offset..]
        .iter()
        .position(|&c| c == b'\0')
        .unwrap_or(data.len() - offset)
        + offset;
    std::str::from_utf8(&data[offset..end]).unwrap_or("")
}

fn print_mpa_style_report_line(
    file: &mut File,
    clade_count: u64,
    taxonomy_line: &str,
) -> io::Result<()> {
    writeln!(file, "{}\t{}", taxonomy_line, clade_count)
}

fn mpa_report_dfs(
    taxid: u64,
    file: &mut File,
    report_zeros: bool,
    taxonomy: &Taxonomy, // 假设你有一个 `Taxonomy` 结构
    clade_counts: &HashMap<u64, u64>,
    taxonomy_names: &mut Vec<String>,
) -> io::Result<()> {
    if !report_zeros && *clade_counts.get(&taxid).unwrap_or(&0) == 0 {
        return Ok(());
    }

    let node = &taxonomy.nodes[taxid as usize]; // 假设 nodes 是 Vec<TaxonomyNode>
    let rank = extract_string_from_offset(&taxonomy.rank_data, node.rank_offset as usize);

    let rank_code = match rank {
        "superkingdom" => 'd',
        "kingdom" => 'k',
        "phylum" => 'p',
        "class" => 'c',
        "order" => 'o',
        "family" => 'f',
        "genus" => 'g',
        "species" => 's',
        _ => '\0',
    };

    if rank_code != '\0' {
        let name_str = extract_string_from_offset(&taxonomy.name_data, node.name_offset as usize);

        let name = format!("{}__{}", rank_code, name_str);
        taxonomy_names.push(name);
        let taxonomy_line = taxonomy_names.join("|");
        print_mpa_style_report_line(
            file,
            *clade_counts.get(&taxid).unwrap_or(&0),
            &taxonomy_line,
        )?;
    }

    let child_count = node.child_count as usize;
    if child_count != 0 {
        let mut children: Vec<u64> = (0..child_count)
            .map(|i| node.first_child + i as u64)
            .collect();

        children.sort_by(|&a, &b| {
            clade_counts
                .get(&b)
                .unwrap_or(&0)
                .cmp(&clade_counts.get(&a).unwrap_or(&0))
        });

        for child in children {
            mpa_report_dfs(
                child,
                file,
                report_zeros,
                taxonomy,
                clade_counts,
                taxonomy_names,
            )?;
        }
    }

    if rank_code != '\0' {
        taxonomy_names.pop();
    }

    Ok(())
}

pub fn report_mpa_style<P: AsRef<Path>>(
    filename: P,
    report_zeros: bool,
    taxonomy: &Taxonomy,
    call_counters: &HashMap<u64, ReadCounter>,
) -> io::Result<()> {
    let call_counts: HashMap<u64, u64> = call_counters
        .iter()
        .map(|(&taxid, counter)| (taxid, counter.read_count()))
        .collect();

    let clade_counts = get_clade_counts(taxonomy, &call_counts);

    let mut file = File::create(filename)?;
    let mut taxonomy_names: Vec<String> = Vec::new();

    mpa_report_dfs(
        1,
        &mut file,
        report_zeros,
        taxonomy,
        &clade_counts,
        &mut taxonomy_names,
    )
}

pub fn print_kraken_style_report_line(
    file: &mut File,
    report_kmer_data: bool,
    total_seqs: u64,
    clade_counter: &mut ReadCounter,
    taxon_counter: &ReadCounter,
    rank_str: &str,
    taxid: u32,
    sci_name: &str,
    depth: usize,
) -> io::Result<()> {
    let pct = 100.0 * clade_counter.read_count() as f64 / total_seqs as f64;
    let pct_str = format!("{:6.2}", pct);

    write!(
        file,
        "{}\t{}\t{}",
        pct_str,
        clade_counter.read_count(),
        taxon_counter.read_count()
    )?;

    if report_kmer_data {
        write!(
            file,
            "\t{}\t{}",
            clade_counter.kmer_count(),
            &clade_counter.distinct_kmer_count()
        )?;
    }

    write!(file, "\t{}\t{}", rank_str, taxid)?;

    for _ in 0..depth {
        write!(file, "  ")?;
    }

    writeln!(file, "{}", sci_name)
}

pub fn kraken_report_dfs(
    taxid: u64,
    file: &mut File,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    clade_counters: &mut HashMap<u64, ReadCounter>,
    call_counters: &HashMap<u64, ReadCounter>,
    total_seqs: u64,
    rank_code: char,
    rank_depth: i32,
    depth: usize,
) -> io::Result<()> {
    if !report_zeros && clade_counters.get(&taxid).map_or(0, |c| c.read_count()) == 0 {
        return Ok(());
    }

    let node = &taxonomy.nodes[taxid as usize];
    let rank = std::str::from_utf8(&taxonomy.rank_data[node.rank_offset as usize..])
        .unwrap_or_default()
        .split('\0')
        .next()
        .unwrap_or("");

    let (new_rank_code, new_rank_depth) = match rank {
        "superkingdom" => ('D', 0),
        "kingdom" => ('K', 0),
        "phylum" => ('P', 0),
        "class" => ('C', 0),
        "order" => ('O', 0),
        "family" => ('F', 0),
        "genus" => ('G', 0),
        "species" => ('S', 0),
        _ => (rank_code, rank_depth + 1),
    };

    let rank_str = if new_rank_depth == 0 {
        new_rank_code.to_string()
    } else {
        format!("{}{}", new_rank_code, new_rank_depth)
    };

    let name = std::str::from_utf8(&taxonomy.name_data[node.name_offset as usize..])
        .unwrap_or_default()
        .split('\0')
        .next()
        .unwrap_or("");

    let mut clade_counter = clade_counters.get_mut(&taxid).unwrap();
    print_kraken_style_report_line(
        file,
        report_kmer_data,
        total_seqs,
        &mut clade_counter,
        call_counters.get(&taxid).unwrap(),
        &rank_str,
        node.external_id as u32,
        name,
        depth,
    )?;

    let mut children: Vec<u64> = (0..node.child_count)
        .map(|i| node.first_child + i)
        .collect();

    children.sort_by_key(|&child_taxid| {
        clade_counters
            .get(&child_taxid)
            .map_or(0, |c| c.read_count())
    });
    children.reverse();

    for child_taxid in children {
        kraken_report_dfs(
            child_taxid,
            file,
            report_zeros,
            report_kmer_data,
            taxonomy,
            clade_counters,
            call_counters,
            total_seqs,
            new_rank_code,
            new_rank_depth,
            depth + 1,
        )?;
    }

    Ok(())
}

pub fn report_kraken_style(
    filename: &str,
    report_zeros: bool,
    report_kmer_data: bool,
    taxonomy: &Taxonomy,
    call_counters: &HashMap<u64, ReadCounter>,
    total_seqs: u64,
    total_unclassified: u64,
) -> io::Result<()> {
    let mut clade_counters = get_clade_counters(taxonomy, call_counters);

    let mut file = File::create(filename)?;

    // 处理未分类序列的特殊情况
    if total_unclassified != 0 || report_zeros {
        let mut rc = ReadCounter::new(total_unclassified, 0);
        let trc = ReadCounter::new(total_unclassified, 0);
        print_kraken_style_report_line(
            &mut file,
            report_kmer_data,
            total_seqs,
            &mut rc,
            &trc,
            "U",
            0,
            "unclassified",
            0,
        )?;
    }

    // 通过 DFS 遍历分类树
    kraken_report_dfs(
        1,
        &mut file,
        report_zeros,
        report_kmer_data,
        taxonomy,
        &mut clade_counters,
        call_counters,
        total_seqs,
        'R',
        -1,
        0,
    )
}
