#[path = "common/mod.rs"]
mod common;

use kun_peng::taxonomy::Taxonomy;
use kun_peng::utils::read_id_to_taxon_map;
use std::collections::BTreeMap;
use std::io;
use std::path::Path;

fn read_cstring(data: &[u8], offset: u64) -> &str {
    let start = offset as usize;
    let end = data[start..]
        .iter()
        .position(|&c| c == b'\0')
        .map(|relative| start + relative)
        .unwrap_or(data.len());
    std::str::from_utf8(&data[start..end]).unwrap_or("")
}

fn lineage_for_taxid(taxonomy: &Taxonomy, taxid: u64) -> Option<String> {
    let internal_id = taxonomy.get_internal_id(taxid);
    if internal_id == 0 {
        return None;
    }

    let path = taxonomy.path_cache.get(&internal_id)?;
    let mut pieces = Vec::with_capacity(path.len());
    for internal in path {
        let node = taxonomy.nodes.get(*internal as usize)?;
        let rank = read_cstring(&taxonomy.rank_data, node.rank_offset);
        let name = read_cstring(&taxonomy.name_data, node.name_offset);
        pieces.push(format!("{}:{}", rank, name));
    }

    Some(pieces.join(" > "))
}

fn require_file(path: &Path, description: &str) -> io::Result<()> {
    if path.exists() {
        return Ok(());
    }
    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!("Required {} `{}` not found.", description, path.display()),
    ))
}

fn main() -> io::Result<()> {
    let workspace_root = common::workspace_root();
    let database_dir = workspace_root.join("test_database");
    require_file(&database_dir, "database directory")?;

    let taxonomy_path = database_dir.join("taxo.k2d");
    let seqid_map_path = database_dir.join("seqid2taxid.map");
    require_file(&taxonomy_path, "taxonomy file")?;
    require_file(&seqid_map_path, "seqid2taxid map")?;

    let taxonomy = Taxonomy::from_file(taxonomy_path)?;
    let id_map = read_id_to_taxon_map(seqid_map_path)?;

    println!("Loaded taxonomy with {} nodes", taxonomy.node_count());
    println!("Showing lineage for a few reference sequences:\n");

    let ordered: BTreeMap<_, _> = id_map.into_iter().collect();
    for (seq_id, taxid) in ordered.iter().take(5) {
        let lineage = lineage_for_taxid(&taxonomy, *taxid)
            .unwrap_or_else(|| "<unknown lineage>".to_string());
        println!("{} (taxid: {}):\n  {}\n", seq_id, taxid, lineage);
    }

    Ok(())
}
