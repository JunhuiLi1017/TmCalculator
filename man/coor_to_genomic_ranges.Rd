% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_genomic_ranges.R
\name{coor_to_genomic_ranges}
\alias{coor_to_genomic_ranges}
\title{Convert genomic coordinate strings to GenomicRanges object with sequences}
\usage{
coor_to_genomic_ranges(input_seq)
}
\arguments{
\item{input_seq}{A character vector where each element is a string in the format:
"chr:start-end:strand:species[:name]"
- chr: Chromosome ID (e.g., "chr1", "chrX")
- start: Start position (integer)
- end: End position (integer)
- strand: "+" for positive strand or "-" for negative strand
- species: Species name for reference genome (e.g., "hg38")
- name: (optional) Custom name for the sequence}
}
\value{
A GenomicRanges object containing:
  - GRanges information (seqnames, ranges, strand)
  - sequence data from the reference genome
  - Names either from the optional name parameter or auto-generated as "1", "2", etc.
}
\description{
This function converts genomic coordinate strings in the format "chr:start-end:strand:species[:name]"
to a GenomicRanges object containing the corresponding sequences from the specified reference genome.
}
\examples{
\dontrun{
# Example with multiple coordinates
coords <- c(
  "chr1:1898000-1898050:+:BSgenome.Hsapiens.UCSC.hg38:exon1",
  "chr2:2563000-2563050:-:BSgenome.Hsapiens.UCSC.hg38:exon2"
)
gr <- coor_to_genomic_ranges(coords)
}

}
