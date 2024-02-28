use kr2r::taxonomy::{NCBITaxonomy, Taxonomy};
use kr2r::utils::read_id_to_taxon_map;
use std::io::Result;

fn main() -> Result<()> {
    // let t = Taxonomy::from_file("taxo.k2d")?;
    // println!("t {:?}", t);
    let mut ncbi = NCBITaxonomy::from_ncbi("lib/taxonomy/nodes.dmp", "lib/taxonomy/names.dmp")?;
    let id_map = read_id_to_taxon_map("lib/seqid2taxid.map")?;
    // let name_map = parse_names_file("lib/taxonomy/names.dmp")?;
    // println!("names {:?}", name_map);

    for (_, id) in id_map.into_iter() {
        ncbi.mark_node(id);
    }
    // println!("ncbi {:?}", ncbi.marked_nodes);
    // println!("name map {:?}", ncbi.name_map);
    let mut taxo = ncbi.convert_to_kraken_taxonomy();
    taxo.generate_external_to_internal_id_map();
    taxo.build_path_cache();
    println!("path_cache {:?}", taxo.path_cache.keys().len());
    println!("node_count {:?}", taxo.node_count());
    // taxo.write_to_disk("taxo.k2d")?;
    Ok(())
}
