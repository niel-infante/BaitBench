pub mod reader;
pub mod writer;

pub use reader::{
    count_reads, count_sequences, parse_fasta, parse_fasta_ids, parse_reads, parse_reads_ids,
};
pub use writer::{concatenate_fastas, extract_by_ids, extract_reads_by_ids};
