% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_genomic_ranges.R
\name{fa_to_genomic_ranges}
\alias{fa_to_genomic_ranges}
\title{Convert FASTA file to GenomicRanges object}
\usage{
fa_to_genomic_ranges(input_seq)
}
\arguments{
\item{input_seq}{Path to the input FASTA file}
}
\value{
A GenomicRanges object containing:
  - GRanges information (seqnames, ranges, strand)
  - sequence data from FASTA file
  - Complementary sequences (if provided)
  - Names from FASTA headers
}
\description{
This function reads sequences from a FASTA file and converts them to a GenomicRanges object. If named with format ">chr2:1-10:[+|-]:[seq_name]", the name will be parsed into GRanges components.
}
\examples{
# Example with single FASTA file
input_seq <- system.file("extdata", "example1.fasta", package = "TmCalculator")
gr <- fa_to_genomic_ranges(input_seq)

}
