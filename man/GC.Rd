% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gc.R
\name{gc}
\alias{gc}
\title{Calculate G and C content of nucleotide sequences}
\usage{
gc(input_seq, ambiguous = FALSE)
}
\arguments{
\item{input_seq}{Sequence (5' to 3') of one strand of the nucleic acid duplex. Can be provided as either:
- A character string (e.g., "ATGCG")
- A path to a FASTA file containing the sequence(s)}

\item{ambiguous}{Logical. If TRUE, ambiguous bases are taken into account when computing the G and C content.
The function handles various ambiguous bases (S, W, M, K, R, Y, V, H, D, B) by proportionally
distributing their contribution to GC content based on their possible nucleotide compositions.
For example:
- S (G or C) contributes fully to GC content
- W (A or T) contributes fully to AT content
- M (A or C) contributes proportionally based on the ratio of A to C in the sequence
- And so on for other ambiguous bases}
}
\value{
Content of G and C as a percentage (range from 0 to 100%)
}
\description{
Calculate G and C content of nucleotide sequences. The function calculates the percentage of G and C bases
relative to the total number of A, T, G, and C bases in the sequence.
}
\examples{

# Calculate GC content of a DNA sequence
gc(c("a","t","c","t","g","g","g","c","c","a","g","t","a"))  # 53.85\%

# Calculate GC content including ambiguous bases
gc("GCATSWSYK", ambiguous = TRUE)  # 55.56\%

}
\author{
Junhui Li
}
