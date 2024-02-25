pub mod compact_hash;
pub mod iclassify;
pub mod kr2r_data;
mod kv_store;
pub mod mmscanner;
pub mod pair;
pub mod readcounts;
pub mod report;
pub mod taxonomy;
pub mod utils;

mod feat;

pub use feat::constants::*;
pub use feat::*;

pub mod db;
pub use kr2r_data::*;
pub use kv_store::*;
// pub use mmscanner::{Meros, CURRENT_REVCOM_VERSION, DEFAULT_TOGGLE_MASK};
pub use readcounts::TaxonCounts;
