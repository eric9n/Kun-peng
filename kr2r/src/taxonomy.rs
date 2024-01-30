use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind, Read, Result, Write};

/// 解析 ncbi 文件的 taxonomy nodes 文件
pub fn parse_nodes_file(
    nodes_filename: &str,
) -> Result<(
    HashMap<u64, u64>,
    HashMap<u64, HashSet<u64>>,
    HashMap<u64, String>,
    HashSet<String>,
)> {
    let nodes_file = File::open(nodes_filename)?;
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

/// 解析 ncbi 文件的 taxonomy names 文件
pub fn parse_names_file(names_filename: &str) -> Result<HashMap<u64, String>> {
    let names_file = File::open(names_filename)?;
    let reader = BufReader::new(names_file);

    let mut name_map = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        // 忽略空行或注释行
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let line = line.trim_end_matches(|c| c == '\t' || c == '|' || c == '\n');
        // 分割行为字段
        let fields: Vec<_> = line.split("\t|\t").collect();
        if fields.len() < 4 {
            continue; // 如果不满足预期的字段数量，则跳过此行
        }
        // 解析节点 ID 和名称类型
        let node_id = fields[0].parse::<u64>().unwrap_or(0);
        let name = fields[1].to_string();
        let name_type = fields[3].to_string();

        // 仅当类型为 "scientific name" 时，将名称添加到 map 中
        if name_type == "scientific name" {
            name_map.insert(node_id, name);
        }
    }

    Ok(name_map)
}

/// 结构体定义
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

// NCBITaxonomy 类型定义
pub struct NCBITaxonomy {
    parent_map: HashMap<u64, u64>,
    name_map: HashMap<u64, String>,
    rank_map: HashMap<u64, String>,
    child_map: HashMap<u64, HashSet<u64>>,
    marked_nodes: HashSet<u64>,
    known_ranks: HashSet<String>,
}

impl NCBITaxonomy {
    // 构造函数等实现
    pub fn from_ncbi(nodes_filename: &str, names_filename: &str) -> Result<Self> {
        let mut marked_nodes = HashSet::new();
        let (parent_map, child_map, rank_map, known_ranks) = parse_nodes_file(nodes_filename)?;

        let name_map = parse_names_file(names_filename)?;

        marked_nodes.insert(1); // 标记根节点

        Ok(NCBITaxonomy {
            parent_map,
            name_map,
            rank_map,
            child_map,
            known_ranks,
            marked_nodes,
        })
    }

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

    pub fn convert_to_kraken_taxonomy(&self) -> Taxonomy {
        let mut taxo = Taxonomy::default();
        // 预分配内存
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

// Taxonomy 类型定义
#[derive(Debug)]
pub struct Taxonomy {
    pub nodes: Vec<TaxonomyNode>,
    pub name_data: Vec<u8>, // 字符串数据以 Vec<u8> 存储
    pub rank_data: Vec<u8>, // 字符串数据以 Vec<u8> 存储
    external_to_internal_id_map: HashMap<u64, u64>,
}

impl Default for Taxonomy {
    fn default() -> Self {
        Taxonomy {
            nodes: Vec::new(),
            name_data: Vec::new(),
            rank_data: Vec::new(),
            external_to_internal_id_map: HashMap::new(),
        }
    }
}

impl Taxonomy {
    const MAGIC: &'static [u8] = b"K2TAXDAT"; // 替换为实际的 magic bytes

    pub fn from_file(filename: &str) -> Result<Taxonomy> {
        let mut file = std::fs::File::open(filename)?;

        let mut magic = vec![0; Self::MAGIC.len()];
        file.read_exact(&mut magic)?;
        if magic != Self::MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Malformed taxonomy file {}", filename),
            ));
        }

        let mut buffer = [0; 24];
        file.read_exact(&mut buffer)?;
        let (node_count, name_data_len, rank_data_len) =
            unsafe { std::mem::transmute::<[u8; 24], (u64, u64, u64)>(buffer) };

        let mut nodes = Vec::with_capacity(node_count as usize);
        for _ in 0..node_count {
            let mut buffer = [0; 56];
            file.read_exact(&mut buffer)?;
            let node = unsafe { std::mem::transmute::<[u8; 56], TaxonomyNode>(buffer) };
            nodes.push(node);
        }

        let mut name_data = vec![0; name_data_len as usize];
        file.read_exact(&mut name_data)?;

        let mut rank_data = vec![0; rank_data_len as usize];
        file.read_exact(&mut rank_data)?;

        let mut external_to_internal_id_map = HashMap::new();
        for (internal_id, node) in nodes.iter().enumerate() {
            let external_id = node.external_id;
            external_to_internal_id_map.insert(external_id, internal_id as u64);
        }

        Ok(Taxonomy {
            nodes,
            name_data,
            rank_data,
            external_to_internal_id_map,
        })
    }

    pub fn is_a_ancestor_of_b(&self, a: u64, b: u64) -> bool {
        if a == 0 || b == 0 {
            return false;
        }

        let mut current = b;

        while current > a {
            current = match self.nodes.get(current as usize) {
                Some(node) => node.parent_id,
                None => return false,
            };
        }

        current == a
    }

    pub fn lowest_common_ancestor(&self, mut a: u64, mut b: u64) -> u64 {
        // 如果任何一个节点是 0，返回另一个节点
        if a == 0 || b == 0 {
            return if a != 0 { a } else { b };
        }

        // 遍历节点直到找到共同的祖先
        while a != b {
            if a > b {
                a = self.nodes.get(a as usize).map_or(0, |node| node.parent_id);
            } else {
                b = self.nodes.get(b as usize).map_or(0, |node| node.parent_id);
            }
        }

        a
    }

    pub fn get_internal_id(&self, external_id: u64) -> u64 {
        self.external_to_internal_id_map
            .get(&external_id)
            .cloned()
            .unwrap_or(0)
    }

    pub fn generate_external_to_internal_id_map(&mut self) {
        self.external_to_internal_id_map.clear();
        self.external_to_internal_id_map.insert(0, 0);

        for (i, node) in self.nodes.iter().enumerate() {
            self.external_to_internal_id_map
                .insert(node.external_id, i as u64);
        }
    }

    pub fn write_to_disk(&self, filename: &str) -> Result<()> {
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
