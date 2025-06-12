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
#' @importFrom Biostrings getSeq
#' @importFrom GenomicRanges GRanges
#' @importFrom seqinr read.fasta
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
    stop("Input sequence must be a character string (e.g., c('ATGCG', 'GCTAG') ), a character vector of genomic coordinate (e.g., 'chr1:100-200:+:BSgenome.Hsapiens.UCSC.hg38'), or a FASTA file")
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
  seq_length <- length(input_seq)

  # Check if name matches the pattern "chr:start-end:strand:name"
  suppressWarnings({
    gr_merged <- do.call(c, sapply(seq_along(input_seq), function(x) {
      sub_seq_name <- seq_name[x]
      # if the name is not null, parse the name
      if (!is.null(sub_seq_name)) {        # if the name matches the pattern "chr:start-end:strand:name", parse the name
        if (grepl("^[^:]+:[0-9]+-[0-9]+:[+-\\*]:[^:]+:[^:]$", sub_seq_name)) {
          parts <- strsplit(sub_seq_name, ":")[[1]]
          range_parts <- as.integer(strsplit(parts[2], "-")[[1]])
          # Validate start and end positions
          if (as.integer(range_parts[1]) > as.integer(range_parts[2])) {
            stop("Start positions must be less than or equal to end positions")
          }
          gr <- GenomicRanges::GRanges(
            seqnames = parts[1],
            ranges = IRanges(start = range_parts[1], end = range_parts[2]),
            strand = parts[3]
          )
          names(gr) <- x
          GenomeInfoDb::genome(gr) <- strsplit(parts[4], "\\.")[[1]][4]
        # if the name matches the pattern "chr:start-end:strand", parse the name
        } else if (grepl("^[^:]+:[0-9]+-[0-9]+:[+-\\*]:[^:]+$", sub_seq_name)) {
          parts <- strsplit(sub_seq_name, ":")[[1]]
          range_parts <- as.integer(strsplit(parts[2], "-")[[1]])
          gr <- GenomicRanges::GRanges(
            seqnames = parts[1],
            ranges = IRanges(start = range_parts[1], end = range_parts[2]),
            strand = parts[3]
          )
          names(gr) <- x
          GenomeInfoDb::genome(gr) <- strsplit(parts[4], "\\.")[[1]][4]
        } else if (grepl("^[^:]+:[0-9]+-[0-9]+:[+-\\*]$", sub_seq_name)) {
          parts <- strsplit(sub_seq_name, ":")[[1]]
          range_parts <- as.integer(strsplit(parts[2], "-")[[1]])
          gr <- GenomicRanges::GRanges(
            seqnames = parts[1],
            ranges = IRanges(start = range_parts[1], end = range_parts[2]),
            strand = parts[3]
          )
          names(gr) <- x
        # if the name matches the pattern "chr:start-end", parse the name
        } else if (grepl("^[^:]+:[0-9]+-[0-9]+$", sub_seq_name)) {
          parts <- strsplit(sub_seq_name, ":")[[1]]
          range_parts <- as.integer(strsplit(parts[2], "-")[[1]])
          gr <- GenomicRanges::GRanges(
            seqnames = parts[1],
            ranges = IRanges(start = range_parts[1], end = range_parts[2]),
            strand = "*"
          )
          names(gr) <- x
        # if the name is not matched, use the default values
        } else {
          gr <- GRanges(
            seqnames = "chr1",
            ranges = IRanges(start = 1, end = nchar(input_seq[x])),
            strand = "*"
          )
          names(gr) <- x
        }
      } else {
        gr <- GRanges(
          seqnames = "chr1",
          ranges = IRanges(start = 1, end = nchar(input_seq[x])),
          strand = "*"
        )
        names(gr) <- x
      }
      S4Vectors::mcols(gr)$sequence <- input_seq[x]
      return(gr)
    }))
  })
  return(gr_merged)
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
  
  # Read sequences from FASTA file
  seq_list <- seqinr::read.fasta(input_seq, as.string = TRUE, forceDNAtolower = FALSE)
  if (length(seq_list) == 0) {
    stop("No sequences found in the FASTA file")
  }
  
  # Convert to named character vector
  seq_vector <- unlist(lapply(seq_list, as.character))
  
  # Create the GenomicRanges object
  gr <- vec_to_genomic_ranges(seq_vector)
  return(gr)
}

#' Convert genomic coordinate strings to GenomicRanges object with sequences
#' 
#' This function converts genomic coordinate strings in the format "chr:start-end:strand:species[:name]"
#' to a GenomicRanges object containing the corresponding sequences from the specified reference genome.
#' 
#' @param input_seq A character vector where each element is a string in the format:
#'   "chr:start-end:strand:species[:name]"
#'   - chr: Chromosome ID (e.g., "chr1", "chrX")
#'   - start: Start position (integer)
#'   - end: End position (integer)
#'   - strand: "+" for positive strand or "-" for negative strand
#'   - species: Species name for reference genome (e.g., "hg38")
#'   - name: (optional) Custom name for the sequence
#' 
#' @return A GenomicRanges object containing:
#'   - GRanges information (seqnames, ranges, strand)
#'   - sequence data from the reference genome
#'   - Names either from the optional name parameter or auto-generated as "1", "2", etc.
#' 
#' @examples
#' \dontrun{
#' # Example with multiple coordinates
#' coords <- c(
#'   "chr1:1898000-1898050:+:BSgenome.Hsapiens.UCSC.hg38:exon1",
#'   "chr2:2563000-2563050:-:BSgenome.Hsapiens.UCSC.hg38:exon2"
#' )
#' gr <- coor_to_genomic_ranges(coords)
#' }
#' 
#' @export

coor_to_genomic_ranges <- function(input_seq){
  all_genomes <- BSgenome::available.genomes()
  suppressWarnings({
    genomic_range_object <- do.call(c,sapply(seq_along(input_seq), function(i){
      x <- input_seq[i] # Get the current input string using the index 'i'
      parts <- strsplit(x, ":")[[1]]
      chrom_x <- parts[1]
      position_x <- as.integer(unlist(strsplit(parts[2],"-")))
      strand_x <- parts[3]
      ref_genome_pkg_name <- parts[4]
      
      # Extract name_x if 5 parts, otherwise it will be NA/null for naming purposes
      name_x <- if(length(parts) == 5) parts[5] else i
      
      # --- Validation and loading of genome package ---
      if (!ref_genome_pkg_name %in% all_genomes) {
        stop(sprintf("Genome package %s not found. Please make sure the genome package name is correct. See BSgenome::available.genomes() for all available genomes.", ref_genome_pkg_name))
      }
      if (!requireNamespace(ref_genome_pkg_name, quietly = TRUE)) {
        stop(sprintf("Genome package %s not found. Please install it first using BiocManager::install(\"%s\").", ref_genome_pkg_name, ref_genome_pkg_name))
      }
      #suppressPackageStartupMessages(
        library(ref_genome_pkg_name, character.only = TRUE)
      #)
      
      genome <- get(ref_genome_pkg_name)
      
      # --- Create GRanges object and fetch sequence ---
      sub_genomic_range <- GenomicRanges::GRanges(
        seqnames = chrom_x,
        ranges = IRanges(start = position_x[1], end = position_x[2]),
        strand = strand_x
      )
      names(sub_genomic_range) <- name_x
      
      # Fetch the sequence
      sub_genomic_range$sequence <- Biostrings::getSeq(genome, sub_genomic_range)
      sub_genomic_range$genome <- ref_genome_pkg_name
      
      return(sub_genomic_range)
    }))
  })
  return(genomic_range_object)
}