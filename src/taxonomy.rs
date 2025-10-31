use crate::utils::open_file;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind, Read, Result, Write};
use std::path::Path;

/// Parse the NCBI taxonomy nodes file
///
/// # Arguments
///
/// * `nodes_filename` - Path to the nodes file
///
/// # Returns
///
/// A tuple containing:
/// - HashMap of node ID to parent ID
/// - HashMap of parent ID to set of child IDs
/// - HashMap of node ID to rank
/// - HashSet of known ranks
pub fn parse_nodes_file<P: AsRef<Path>>(
    nodes_filename: P,
) -> Result<(
    HashMap<u64, u64>,
    HashMap<u64, HashSet<u64>>,
    HashMap<u64, String>,
    HashSet<String>,
)> {
    let nodes_file = open_file(nodes_filename)?;
    let reader = BufReader::new(nodes_file);

    let mut parent_map = HashMap::new();
    let mut child_map = HashMap::new();
    let mut rank_map = HashMap::new();
    let mut known_ranks = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<_> = line.split("\t|\t").collect();
        if fields.len() < 3 {
            continue;
        }

        let node_id = fields[0]
            .parse::<u64>()
            .map_err(|_| Error::new(ErrorKind::InvalidData, "node_id"))?;

        let parent_id = if node_id == 1 {
            0
        } else {
            fields[1]
                .parse::<u64>()
                .map_err(|_| Error::new(ErrorKind::InvalidData, "parent_id"))?
        };

        let rank = fields[2].to_string();

        parent_map.insert(node_id, parent_id);
        child_map
            .entry(parent_id)
            .or_insert_with(HashSet::new)
            .insert(node_id);
        rank_map.insert(node_id, rank.clone());
        known_ranks.insert(rank);
    }

    Ok((parent_map, child_map, rank_map, known_ranks))
}

/// Parse the NCBI taxonomy names file
///
/// # Arguments
///
/// * `names_filename` - Path to the names file
///
/// # Returns
///
/// A HashMap of node ID to scientific name
pub fn parse_names_file<P: AsRef<Path>>(names_filename: P) -> Result<HashMap<u64, String>> {
    let names_file = open_file(names_filename)?;
    let reader = BufReader::new(names_file);

    let mut name_map = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        // Ignore empty lines or comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let line = line.trim_end_matches(|c| c == '\t' || c == '|' || c == '\n');
        // Split the line into fields
        let fields: Vec<_> = line.split("\t|\t").collect();
        if fields.len() < 4 {
            continue; // Skip if the line doesn't have the expected number of fields
        }
        // Parse node ID and name type
        let node_id = fields[0].parse::<u64>().unwrap_or(0);
        let name = fields[1].to_string();
        let name_type = fields[3].to_string();

        // Only add the name to the map if it's a "scientific name"
        if name_type == "scientific name" {
            name_map.insert(node_id, name);
        }
    }

    Ok(name_map)
}

/// Represents a node in the taxonomy
#[derive(Debug)]
pub struct TaxonomyNode {
    pub parent_id: u64,
    pub first_child: u64,
    pub child_count: u64,
    pub name_offset: u64,
    pub rank_offset: u64,
    pub external_id: u64,
    pub godparent_id: u64,
}

impl Default for TaxonomyNode {
    fn default() -> Self {
        Self {
            parent_id: 0,
            first_child: 0,
            child_count: 0,
            name_offset: 0,
            rank_offset: 0,
            external_id: 0,
            godparent_id: 0,
        }
    }
}

// NCBITaxonomy struct definition
pub struct NCBITaxonomy {
    parent_map: HashMap<u64, u64>,
    name_map: HashMap<u64, String>,
    rank_map: HashMap<u64, String>,
    child_map: HashMap<u64, HashSet<u64>>,
    marked_nodes: HashSet<u64>,
    known_ranks: HashSet<String>,
}

impl NCBITaxonomy {
    /// Create a new NCBITaxonomy from NCBI taxonomy files
    ///
    /// # Arguments
    ///
    /// * `nodes_filename` - Path to the nodes file
    /// * `names_filename` - Path to the names file
    ///
    /// # Returns
    ///
    /// A Result containing the new NCBITaxonomy or an error
    pub fn from_ncbi<P: AsRef<Path>>(nodes_filename: P, names_filename: P) -> Result<Self> {
        let mut marked_nodes = HashSet::new();
        let (parent_map, child_map, rank_map, known_ranks) = parse_nodes_file(nodes_filename)?;

        let name_map = parse_names_file(names_filename)?;

        marked_nodes.insert(1); // Mark the root node

        Ok(NCBITaxonomy {
            parent_map,
            name_map,
            rank_map,
            child_map,
            known_ranks,
            marked_nodes,
        })
    }

    /// Mark a node and all its ancestors in the taxonomy
    ///
    /// # Arguments
    ///
    /// * `taxid` - The ID of the node to mark
    pub fn mark_node(&mut self, taxid: u64) {
        let mut current_taxid = taxid;
        while !self.marked_nodes.contains(&current_taxid) {
            self.marked_nodes.insert(current_taxid);
            if let Some(&parent_id) = self.parent_map.get(&current_taxid) {
                current_taxid = parent_id;
            } else {
                break;
            }
        }
    }

    /// Get rank offset data for the taxonomy
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// - HashMap of rank to offset
    /// - String containing all ranks
    pub fn get_rank_offset_data(&self) -> (HashMap<String, u64>, String) {
        let mut rank_data = String::new();
        let mut rank_offsets = HashMap::new();

        let mut known_ranks: Vec<_> = self.known_ranks.iter().collect();
        known_ranks.sort_unstable();

        for rank in known_ranks {
            rank_offsets.insert(rank.clone(), rank_data.len() as u64);
            rank_data.push_str(rank);
            rank_data.push('\0');
        }

        (rank_offsets, rank_data)
    }

    /// Convert the NCBITaxonomy to a Kraken-style Taxonomy
    ///
    /// # Returns
    ///
    /// A new Taxonomy object
    pub fn convert_to_kraken_taxonomy(&self) -> Taxonomy {
        let mut taxo = Taxonomy::default();
        // Preallocate memory
        taxo.nodes.reserve(self.marked_nodes.len() + 1);
        taxo.nodes.push(TaxonomyNode::default());

        let mut name_data = String::new();
        let (rank_offsets, rank_data) = self.get_rank_offset_data();

        let mut bfs_queue = VecDeque::new();
        bfs_queue.push_back(1);
        let mut external_id_map = HashMap::new();
        external_id_map.insert(0, 0);
        external_id_map.insert(1, 1);
        let mut internal_node_id = 0;

        while let Some(external_node_id) = bfs_queue.pop_front() {
            internal_node_id += 1;
            external_id_map.insert(external_node_id, internal_node_id);

            let mut node = TaxonomyNode::default();
            node.parent_id = external_id_map
                .get(&self.parent_map[&external_node_id])
                .unwrap()
                .clone();
            node.external_id = external_node_id;
            node.rank_offset = *rank_offsets.get(&self.rank_map[&external_node_id]).unwrap();
            node.name_offset = name_data.len() as u64;

            let name = self.name_map.get(&external_node_id).unwrap();
            name_data.push_str(name);
            name_data.push('\0');

            node.first_child = internal_node_id + bfs_queue.len() as u64 + 1;

            if let Some(children) = self.child_map.get(&external_node_id) {
                let mut sorted_children: Vec<_> = children.iter().collect();
                sorted_children.sort_unstable();

                for &child_node in sorted_children {
                    if self.marked_nodes.contains(&child_node) {
                        bfs_queue.push_back(child_node);
                        node.child_count += 1;
                    }
                }
            }
            taxo.nodes.push(node);
        }

        taxo.name_data = name_data.into_bytes();
        taxo.rank_data = rank_data.into_bytes();

        taxo
    }
}

// Taxonomy struct definition
#[derive(Debug)]
pub struct Taxonomy {
    pub path_cache: HashMap<u32, Vec<u32>>,
    pub nodes: Vec<TaxonomyNode>,
    pub name_data: Vec<u8>, // String data stored as Vec<u8>
    pub rank_data: Vec<u8>, // String data stored as Vec<u8>
    external_to_internal_id_map: HashMap<u64, u32>,
}

impl Default for Taxonomy {
    fn default() -> Self {
        Taxonomy {
            path_cache: HashMap::new(),
            nodes: Vec::new(),
            name_data: Vec::new(),
            rank_data: Vec::new(),
            external_to_internal_id_map: HashMap::new(),
        }
    }
}

impl Taxonomy {
    const MAGIC: &'static [u8] = b"K2TAXDAT"; // Replace with actual magic bytes

    /// Create a new Taxonomy from a file
    ///
    /// # Arguments
    ///
    /// * `filename` - Path to the taxonomy file
    ///
    /// # Returns
    ///
    /// A Result containing the new Taxonomy or an error
    pub fn from_file<P: AsRef<Path> + Debug>(filename: P) -> Result<Taxonomy> {
        let mut file = open_file(&filename)?;

        let mut magic = vec![0; Self::MAGIC.len()];
        file.read_exact(&mut magic)?;
        if magic != Self::MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Malformed taxonomy file {:?}", &filename),
            ));
        }

        let mut header = [0u8; 24];
        file.read_exact(&mut header)?;
        let node_count = u64::from_le_bytes(header[0..8].try_into().unwrap());
        let name_data_len = u64::from_le_bytes(header[8..16].try_into().unwrap());
        let rank_data_len = u64::from_le_bytes(header[16..24].try_into().unwrap());

        let mut nodes = Vec::with_capacity(node_count as usize);
        for _ in 0..node_count {
            let mut buffer = [0u8; 56];
            file.read_exact(&mut buffer)?;
            let parent_id = u64::from_le_bytes(buffer[0..8].try_into().unwrap());
            let first_child = u64::from_le_bytes(buffer[8..16].try_into().unwrap());
            let child_count = u64::from_le_bytes(buffer[16..24].try_into().unwrap());
            let name_offset = u64::from_le_bytes(buffer[24..32].try_into().unwrap());
            let rank_offset = u64::from_le_bytes(buffer[32..40].try_into().unwrap());
            let external_id = u64::from_le_bytes(buffer[40..48].try_into().unwrap());
            let godparent_id = u64::from_le_bytes(buffer[48..56].try_into().unwrap());

            nodes.push(TaxonomyNode {
                parent_id,
                first_child,
                child_count,
                name_offset,
                rank_offset,
                external_id,
                godparent_id,
            });
        }

        let mut name_data = vec![0; name_data_len as usize];
        file.read_exact(&mut name_data)?;

        let mut rank_data = vec![0; rank_data_len as usize];
        file.read_exact(&mut rank_data)?;

        let mut external_to_internal_id_map = HashMap::new();
        for (internal_id, node) in nodes.iter().enumerate() {
            let external_id = node.external_id;
            external_to_internal_id_map.insert(external_id, internal_id as u32);
        }

        let mut taxo = Taxonomy {
            path_cache: HashMap::new(),
            nodes,
            name_data,
            rank_data,
            external_to_internal_id_map,
        };
        taxo.build_path_cache();
        Ok(taxo)
    }

    /// Check if node A is an ancestor of node B
    ///
    /// # Arguments
    ///
    /// * `a` - The potential ancestor node ID
    /// * `b` - The potential descendant node ID
    ///
    /// # Returns
    ///
    /// A boolean indicating whether A is an ancestor of B
    pub fn is_a_ancestor_of_b(&self, a: u32, b: u32) -> bool {
        if a == 0 || b == 0 {
            return false;
        }

        // Try to get the ancestor path of B from the path cache
        if let Some(path) = self.path_cache.get(&b) {
            // Check if the path contains A
            return path.contains(&a);
        }

        false
    }

    /// Find the lowest common ancestor of two nodes
    ///
    /// # Arguments
    ///
    /// * `a` - The first node ID
    /// * `b` - The second node ID
    ///
    /// # Returns
    ///
    /// The ID of the lowest common ancestor
    pub fn lca(&self, a: u32, b: u32) -> u32 {
        if a == 0 || b == 0 || a == b {
            return if a != 0 { a } else { b };
        }

        let default: Vec<u32> = vec![0];
        let path_a = self.path_cache.get(&a).unwrap_or(&default);
        let path_b = self.path_cache.get(&b).unwrap_or(&default);

        let mut i = 0;
        while i < path_a.len() && i < path_b.len() && path_a[i] == path_b[i] {
            i += 1;
        }

        if i == 0 {
            return 0;
        }

        // Return the last common ancestor
        *path_a.get(i - 1).unwrap_or(&0)
    }

    /// Build the path cache for efficient ancestor lookups
    pub fn build_path_cache(&mut self) {
        let mut cache: HashMap<u32, Vec<u32>> = HashMap::new();
        let root_external_id = 1u64;
        if let Some(&root_internal_id) = self.external_to_internal_id_map.get(&root_external_id) {
            // Start traversing from the root node
            self.build_path_for_node(root_internal_id, &mut cache, Vec::new());
        }
        self.path_cache = cache;
    }

    fn build_path_for_node(
        &self,
        node_id: u32,
        path_cache: &mut HashMap<u32, Vec<u32>>,
        mut current_path: Vec<u32>,
    ) {
        current_path.push(node_id); // Add the current node to the path
                                    // Store the current node's path
        path_cache.insert(node_id, current_path.clone());

        // Get the current node's information
        let node = &self.nodes[node_id as usize];
        let first_child_id = node.first_child as u32;
        let child_count = node.child_count as u32;

        // Traverse all children
        for i in 0..child_count {
            let child_internal_id = first_child_id + i; // Assume child IDs are consecutive
            self.build_path_for_node(child_internal_id, path_cache, current_path.clone());
        }
    }

    /// Get the number of nodes in the taxonomy
    ///
    /// # Returns
    ///
    /// The number of nodes
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Get the internal ID for a given external ID
    ///
    /// # Arguments
    ///
    /// * `external_id` - The external ID to look up
    ///
    /// # Returns
    ///
    /// The corresponding internal ID, or 0 if not found
    pub fn get_internal_id(&self, external_id: u64) -> u32 {
        *self
            .external_to_internal_id_map
            .get(&external_id)
            .unwrap_or(&0)
    }

    /// Generate the mapping from external to internal IDs
    pub fn generate_external_to_internal_id_map(&mut self) {
        self.external_to_internal_id_map.clear();
        self.external_to_internal_id_map.insert(0, 0);

        for (i, node) in self.nodes.iter().enumerate() {
            self.external_to_internal_id_map
                .insert(node.external_id, i as u32);
        }
    }

    /// Write the taxonomy to disk
    ///
    /// # Arguments
    ///
    /// * `filename` - Path to write the taxonomy file
    ///
    /// # Returns
    ///
    /// A Result indicating success or failure
    pub fn write_to_disk<P: AsRef<Path>>(&self, filename: P) -> Result<()> {
        let mut file = File::create(filename)?;

        // Write file magic
        file.write_all(Taxonomy::MAGIC)?;

        // Write node count, name data length, and rank data length
        let node_count = self.nodes.len() as u64;
        let name_data_len = self.name_data.len() as u64;
        let rank_data_len = self.rank_data.len() as u64;
        file.write_all(&node_count.to_le_bytes())?;
        file.write_all(&name_data_len.to_le_bytes())?;
        file.write_all(&rank_data_len.to_le_bytes())?;

        // Write nodes as binary data
        for node in &self.nodes {
            file.write_all(&node.parent_id.to_le_bytes())?;
            file.write_all(&node.first_child.to_le_bytes())?;
            file.write_all(&node.child_count.to_le_bytes())?;
            file.write_all(&node.name_offset.to_le_bytes())?;
            file.write_all(&node.rank_offset.to_le_bytes())?;
            file.write_all(&node.external_id.to_le_bytes())?;
            file.write_all(&node.godparent_id.to_le_bytes())?;
        }

        // Write name data and rank data
        file.write_all(&self.name_data)?;
        file.write_all(&self.rank_data)?;

        Ok(())
    }
}
