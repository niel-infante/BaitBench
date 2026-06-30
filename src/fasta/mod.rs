pub mod reader;
pub mod writer;

pub use reader::{parse_fasta, parse_fasta_ids, count_sequences, parse_reads, parse_reads_ids, count_reads};
pub use writer::{extract_by_ids, concatenate_fastas, extract_reads_by_ids};
