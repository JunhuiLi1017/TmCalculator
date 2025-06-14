% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_genomic_ranges.R
\name{vec_to_genomic_ranges}
\alias{vec_to_genomic_ranges}
\title{Convert sequence strings to GenomicRanges object}
\usage{
vec_to_genomic_ranges(input_seq)
}
\arguments{
\item{input_seq}{A character vector of sequences. If named with format "chr2:1-10:[+|-]:[seq_name]" the name will be parsed into GRanges components.}
}
\value{
A GenomicRanges object containing:
  - GRanges information (seqnames, ranges, strand)
  - sequence data
  - Complementary sequences
  - Names from input or auto-generated
}
\description{
This function converts sequence strings to a GenomicRanges object, handling both
named and unnamed sequences. It can also process complementary sequences if provided.
sequence names can be in the format ">chr2:1-10:+:seq2" which will be parsed into
chromosome, position, strand, and name components.
}
\examples{
# Example with named sequences in GRanges format
seqs <- c("ATGCG", "GCTAG")
names(seqs) <- c("chr1:1111-1115:+:seq1", "chr2:1221-1225:+")
gr <- vec_to_genomic_ranges(seqs)

# Example with unnamed sequences
seqs <- c("ATGCG", "GCTAG")
gr <- vec_to_genomic_ranges(seqs)

}
