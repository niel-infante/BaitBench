pub mod reader;
pub mod writer;

pub use reader::{parse_fasta, parse_fasta_ids, count_sequences};
pub use writer::{extract_by_ids, concatenate_fastas};
