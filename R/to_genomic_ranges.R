#' Convert input file into a GenomicRanges Object
#'
#' This function processes a vector of sequences string, a FASTA file, or a character vector with genomic coordinates into a GenomicRanges object,
#' optionally including complementary sequences. sequence names are parsed based on their format:
#' - If names have this pattern "chr:start-end:strand:species[:name]" (e.g., "chr1:1-5:+:seq_1"), parse components into seqnames, ranges, strand, and name.
#' - If names have this pattern "chr:start-end:strand" (e.g., "chr1:1-5:+"), parse components into seqnames, ranges, and strand.
#' - If names have this pattern "chr:start-end" (e.g., "chr1:1-5"), parse components into seqnames and ranges.
#' - If no names are provided, use default values: seqnames = "chr1", start = 1, width = sequence length, strand = "*", name = "1", etc.
#' Complementary sequences are either provided or automatically generated.
#'
#' @param input_seq Input sequence(s) in 5' to 3' direction. Can be provided as either:
#'   - A character string (e.g., c("ATGCG", "GCTAG"))
#'   - A path to a FASTA file containing the sequence(s)
#'   - A character vector where each element is a string in the format "chr:start-end:strand:species" #' (e.g., "chr1:100-200:+:BSgenome.Hsapiens.UCSC.hg38"). Strand is "+" for positive or "-" for negative.
#'     - chr: Chromosome ID
#'     - start: Start position
#'     - end: End position
#'     - strand: positive or negative strand
#'     - species:  Species name for reference genome (e.g., "BSgenome.Hsapiens.UCSC.hg38"), see \code{BSgenome::available.genomes()} for all available genomes. please make sure the genome package is installed, otherwise the function will stop.
#' @param complement_seq Optional complementary sequences. If NULL, complementary sequences will be auto-generated. otherwise, the complementary sequences will be used as metadata. Can be provided as format of input_seq.
#' @return A GenomicRanges object with seqnames, ranges, strand, name, sequence, Complement, and Tm as metadata.
#' @examples
#' # Using a character vector with auto-generated complementary sequences
#' seqs <- c("ATGCG", "GCTAG")
#' names(seqs) <- c("chr1:1-5:+:seq_1", "chr2:1-5:+")
#' gr <- to_genomic_ranges(seqs)
#' gr
#'
#' # Using a character vector with provided complementary sequences
#' seqs <- c("ATGCG", "GCTAG")
#' comp_seqs <- c("TACGC", "CGTA")
#' gr <- to_genomic_ranges(seqs, comp_seqs)
#' gr
#'
#' # Using a FASTA file
#' gr <- to_genomic_ranges(system.file("extdata", "example1.fasta", package = "TmCalculator"))
#' \dontrun{
#' # Using a character vector with genomic coordinates
#' seqs <- c(
#'   "chr1:1898000-1898050:+:BSgenome.Hsapiens.UCSC.hg38",
#'   "chr2:2563000-2563050:-:BSgenome.Hsapiens.UCSC.hg38"
#' )
#' gr <- to_genomic_ranges(seqs)
#' gr
#' }
#' 
#' @encoding UTF-8
#' @author Junhui Li
#' 
#' @export
#' 
#' @importFrom Biostrings getSeq readBStringSet
#' @importFrom GenomicRanges GRanges
#' @importFrom BSgenome available.genomes
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb genome
#' 
#' @export to_genomic_ranges
 
to_genomic_ranges <- function(input_seq, complement_seq = NULL) {
  # Validate input_seq
  if (is.null(input_seq) || length(input_seq) == 0) {
    stop("Input sequence cannot be NULL or empty")
  }
  
  # Process input based on type
  # for fa file
  if (is.character(input_seq) && length(input_seq) == 1 && file.exists(input_seq)) {
    input_gr <- fa_to_genomic_ranges(input_seq)
  # for genomic coordinates
  } else if (is.character(input_seq) && all(grepl(":", input_seq))) {
    input_gr <- coor_to_genomic_ranges(input_seq)
  # for vec of sequence strings
  } else if (is.character(input_seq)) {
    input_gr <- vec_to_genomic_ranges(input_seq)
  } else {
    stop("Input sequence must be a character string (e.g., c('ATGCG', 'GCTAG') ), a character vector of genomic coordinate (e.g., 'chr1:100100-100200:+:BSgenome.Hsapiens.UCSC.hg38'), or a FASTA file")
  }
  
  # Process complementary sequences if provided
  if (!is.null(complement_seq)) {
    # for fa file
    if (is.character(complement_seq) && length(complement_seq) == 1 && file.exists(complement_seq)) {
      input_gr_comp <- fa_to_genomic_ranges(complement_seq)
    # for genomic coordinates
    } else if (is.character(complement_seq) && all(grepl(":", complement_seq))) {
      input_gr_comp <- coor_to_genomic_ranges(complement_seq)
    # for vec of sequence strings
    } else if (is.character(complement_seq)) {
      input_gr_comp <- vec_to_genomic_ranges(complement_seq)
    } else {
      stop("Complementary sequence must be a character string (e.g., c('ATGCG', 'GCTAG') ), a character vector of genomic coordinate (e.g., 'chr1:100-200:+:BSgenome.Hsapiens.UCSC.hg38'), or a FASTA file")
    }
    input_gr$complement <- input_gr_comp$sequence
  } else {
    # Auto-generate complementary sequences
    comp_vector <- generate_complement(as.character(input_gr$sequence))
    input_gr$complement <- comp_vector
  }
  
  return(input_gr)
}

#' Convert sequence strings to GenomicRanges object
#' 
#' This function converts sequence strings to a GenomicRanges object, handling both
#' named and unnamed sequences. It can also process complementary sequences if provided.
#' sequence names can be in the format ">chr2:1-10:+:seq2" which will be parsed into
#' chromosome, position, strand, and name components.
#' 
#' @param input_seq A character vector of sequences. If named with format "chr2:1-10:[+|-]:[seq_name]" the name will be parsed into GRanges components.
#' 
#' @return A GenomicRanges object containing:
#'   - GRanges information (seqnames, ranges, strand)
#'   - sequence data
#'   - Complementary sequences
#'   - Names from input or auto-generated
#' 
#' @examples
#' # Example with named sequences in GRanges format
#' seqs <- c("ATGCG", "GCTAG")
#' names(seqs) <- c("chr1:1111-1115:+:seq1", "chr2:1221-1225:+")
#' gr <- vec_to_genomic_ranges(seqs)
#' 
#' # Example with unnamed sequences
#' seqs <- c("ATGCG", "GCTAG")
#' gr <- vec_to_genomic_ranges(seqs)
#' 
#' @export

vec_to_genomic_ranges <- function(input_seq) {
  if (is.null(input_seq) || length(input_seq) == 0) {
    stop("Input sequence cannot be NULL or empty")
  }
  seq_name <- names(input_seq)
  n <- length(input_seq)

  ## Vectorised construction. The previous implementation built one GRanges
  ## per sequence and combined them with do.call(c, ...), which cost
  ## milliseconds per sequence and dominated genome-scale runs; the parsing
  ## rules below are unchanged.

  ## Defaults for unnamed sequences, or names that match no pattern
  chr    <- rep("chr1", n)
  starts <- rep(1L, n)
  ends   <- as.integer(nchar(input_seq))
  strand <- rep("*", n)
  gen    <- rep(NA_character_, n)

  if (!is.null(seq_name)) {
    nm <- ifelse(is.na(seq_name), "", seq_name)

    ## Patterns are tried from most to least specific, exactly as before.
    p5 <- grepl("^[^:]+:[0-9]+-[0-9]+:[+-\\*]:[^:]+:[^:]$", nm)
    p4 <- !p5 & grepl("^[^:]+:[0-9]+-[0-9]+:[+-\\*]:[^:]+$", nm)
    p3 <- !p5 & !p4 & grepl("^[^:]+:[0-9]+-[0-9]+:[+-\\*]$", nm)
    p2 <- !p5 & !p4 & !p3 & grepl("^[^:]+:[0-9]+-[0-9]+$", nm)
    hit <- p5 | p4 | p3 | p2

    if (any(hit)) {
      idx   <- which(hit)
      parts <- strsplit(nm[idx], ":", fixed = TRUE)
      rng   <- strsplit(vapply(parts, `[`, character(1), 2L), "-", fixed = TRUE)

      chr[idx]    <- vapply(parts, `[`, character(1), 1L)
      starts[idx] <- as.integer(vapply(rng, `[`, character(1), 1L))
      ends[idx]   <- as.integer(vapply(rng, `[`, character(1), 2L))

      ## strand is field 3 for the three patterns that carry it, "*" for p2
      has_strand <- (p5 | p4 | p3)[idx]
      strand[idx][has_strand] <-
        vapply(parts[has_strand], `[`, character(1), 3L)
      strand[idx][!has_strand] <- "*"

      ## genome is taken from field 4 for patterns p5 and p4
      has_gen <- (p5 | p4)[idx]
      if (any(has_gen)) {
        gen[idx][has_gen] <- vapply(
          strsplit(vapply(parts[has_gen], `[`, character(1), 4L), ".",
                   fixed = TRUE),
          function(z) if (length(z) >= 4L) z[4L] else NA_character_,
          character(1))
      }

      bad <- starts[idx] > ends[idx]
      if (any(bad)) {
        stop("Start positions must be less than or equal to end positions")
      }
    }
  }

  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = starts, end = ends),
    strand   = strand
  )
  names(gr) <- as.character(seq_len(n))
  S4Vectors::mcols(gr)$sequence <- input_seq

  ## GenomeInfoDb::genome() is per seqlevel, not per range; set it when the
  ## parsed names agree on a single genome, as the previous per-sequence
  ## code effectively did.
  gen_ok <- gen[!is.na(gen)]
  if (length(gen_ok) && length(unique(gen_ok)) == 1L) {
    suppressWarnings(try(
      GenomeInfoDb::genome(gr) <- unique(gen_ok), silent = TRUE))
  }

  gr
}

#' Convert FASTA file to GenomicRanges object
#' 
#' This function reads sequences from a FASTA file and converts them to a GenomicRanges object. If named with format ">chr2:1-10:[+|-]:[seq_name]", the name will be parsed into GRanges components.
#' 
#' @param input_seq Path to the input FASTA file
#' 
#' @return A GenomicRanges object containing:
#'   - GRanges information (seqnames, ranges, strand)
#'   - sequence data from FASTA file
#'   - Complementary sequences (if provided)
#'   - Names from FASTA headers
#' 
#' @examples
#' # Example with single FASTA file
#' input_seq <- system.file("extdata", "example1.fasta", package = "TmCalculator")
#' gr <- fa_to_genomic_ranges(input_seq)
#' 
#' @export
fa_to_genomic_ranges <- function(input_seq) {
  # Validate input_seq
  if (!file.exists(input_seq)) {
    stop("Input FASTA file does not exist")
  }
  
  # Read with Biostrings, already a hard dependency, rather than seqinr, which
  # was needed for this one call. readBStringSet() is used rather than
  # readDNAStringSet() deliberately: the latter validates against the DNA
  # alphabet and would reject RNA input, but this package ships twelve RNA
  # parameter sets and inst/extdata/example1.fasta contains a U-bearing
  # sequence. BStringSet imposes no alphabet and preserves case, so it matches
  # seqinr::read.fasta(as.string = TRUE, forceDNAtolower = FALSE) exactly.
  dss <- Biostrings::readBStringSet(input_seq)

  if (length(dss) == 0) {
    stop("No sequences found in the FASTA file")
  }

  # seqinr named sequences by the first word of the header; keep that so that
  # downstream region identifiers are unchanged.
  seq_vector <- as.character(dss)
  names(seq_vector) <- sub("\\s.*$", "", names(dss))

  # Create the GenomicRanges object
  gr <- vec_to_genomic_ranges(seq_vector)
  return(gr)
}
