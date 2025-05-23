#' Process sequence input from various formats
#' 
#' @param input_seq Input sequence as a character string or FASTA file path
#' @return A list of sequences with attributes including name, Tm
#' @keywords internal
#' @importFrom seqinr read.fasta
#' @export

process_seq <- function(input_seq) {
  if(is.null(input_seq) || length(input_seq) == 0) {
    stop("Input sequence cannot be NULL or empty")
  }
  
  if(is.character(input_seq) && length(input_seq) == 1 && file.exists(input_seq)) {
    # Read sequences from FASTA file
    seq_list <- seqinr::read.fasta(input_seq)
    if(length(seq_list) == 0) {
      stop("No sequences found in the FASTA file")
    }
    # Assign custom class to each sequence
    for(i in seq_along(seq_list)) {
      class(seq_list[[i]]) <- "TmCalculator"
    }
    return(seq_list)
  } else if(is.character(input_seq)) {
    # Convert string or vector of strings to list
    # Vector of strings input
    seq_list <- as.list(input_seq)
    for(i in seq_along(seq_list)) {
      names(seq_list)[i] <- paste("seq", i, sep = "_")
    }
    # Assign custom class to each sequence
    for(i in seq_along(seq_list)) {
      attr(seq_list[[i]],"name") <- names(seq_list)[i]
      attr(seq_list[[i]],"Tm") <- NA
    }
    return(seq_list)
  } else {
    stop("Input must be either a character string, a vector of character strings, or a path to a FASTA file")
  }
}
