#' Process sequence input from various formats
#' 
#' This function processes input sequences from different formats (direct string input or FASTA files)
#' and prepares them for Tm calculation. It handles both single sequences and multiple sequences,
#' and preserves sequence names from FASTA files.
#' 
#' @param input_seq Input sequence(s) in one of the following formats:
#'   - A character string (e.g., "ATGCG")
#'   - A vector of character strings (e.g., c("ATGCG", "GCTAG"))
#'   - A path to a FASTA file containing one or more sequences
#' 
#' @return A list of sequences with the following attributes:
#'   - name: Sequence name (from FASTA header or auto-generated)
#'   - Tm: Initialized as NA (will be calculated later)
#' 
#' @details
#' For FASTA files:
#' - Sequences are read using seqinr::read.fasta
#' - Sequence names are preserved from FASTA headers
#' - All sequences are converted to uppercase
#' 
#' For direct input:
#' - Single sequences are converted to a list
#' - Multiple sequences are preserved as a list
#' - Names are auto-generated as "seq_1", "seq_2", etc.
#' 
#' @importFrom seqinr read.fasta
#' 
#' @examples
#' # Single sequence
#' process_seq("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
#' 
#' # Multiple sequences
#' process_seq(c("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC",
#'               "AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC"))
#' 
#' # FASTA file
#' process_seq(system.file("extdata", "example1.fasta", package = "TmCalculator"))
#' 
#' @export process_seq

process_seq <- function(input_seq) {
  #' Add attributes to a sequence
  sequence_attr <- function(seq) {
    seq_vec <- as.character(seq)
    attr(seq_vec, "name") <- attr(seq, "name")
    attr(seq_vec, "Tm") <- NA
    return(seq_vec)
  }
  
  if(is.null(input_seq) || length(input_seq) == 0) {
    stop("Input sequence cannot be NULL or empty")
  }
  
  if(is.character(input_seq) && length(input_seq) == 1 && file.exists(input_seq)) {
    # Read sequences from FASTA file
    seq_list <- seqinr::read.fasta(input_seq, as.string = TRUE, forceDNAtolower = FALSE)
    if(length(seq_list) == 0) {
      stop("No sequences found in the FASTA file")
    }
    # Assign custom class to each sequence
    result <- sequence_attr(seq_list)
    return(result)
  } else if(is.character(input_seq)) {
    # Convert string or vector of strings to list
    seq_list <- as.list(input_seq)
    names(seq_list) <- paste("seq", 1:length(seq_list), sep = "_")
    # Assign custom class to each sequence
    result <- sequence_attr(seq_list)
    return(result)
  } else {
    stop("Input must be either a character string, a vector of character strings, or a path to a FASTA file")
  }
}
