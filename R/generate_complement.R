#' Generate complementary sequence
#' 
#' Generate the complementary sequence of a nucleic acid sequence, with an option to reverse it.
#' 
#' @param input_seq Input sequence(s) in 5' to 3' direction. Must be provided as either:
#'   - A character string (e.g., "ATGCG")
#' 
#' @param reverse Logical. If TRUE, the complementary sequence is reversed (3' to 5').
#'   If FALSE (default), the complementary sequence is in the same direction (5' to 3').
#' 
#' @returns Returns the complementary sequence(s) in the specified direction.
#' 
#' @author Junhui Li
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 
#' @examples
#' 
#' # Generate complementary sequence in same direction (5' to 3')
#' generate_complement("ATGCG", reverse = FALSE)
#' 
#' # Generate complementary sequence in reverse direction (3' to 5')
#' generate_complement("ATGCG", reverse = TRUE)
#' 
#' @export generate_complement

generate_complement <- function(input_seq, reverse = FALSE) {
  # Process input sequence
  raw_seq <- process_seq(input_seq)
  
  # Define complement table
  complement_table <- c(
    "A" = "T", "T" = "A", "G" = "C", "C" = "G",
    "M" = "K", "K" = "M", "R" = "Y", "Y" = "R",
    "W" = "W", "S" = "S", "B" = "V", "V" = "B",
    "D" = "H", "H" = "D", "N" = "N", "I" = "I"
  )
  
  # Process each sequence
  result <- lapply(raw_seq, function(seq) {
    # Convert to character vector
    seq_vec <- s2c(seq)
    
    # Get complement
    comp_vec <- sapply(seq_vec, function(base) {
      if (base %in% names(complement_table)) {
        return(complement_table[base])
      } else {
        return(".")
      }
    })
    
    # Reverse if requested
    if (reverse) {
      comp_vec <- rev(comp_vec)
    }
    
    # Convert back to string and preserve attributes
    result <- c2s(comp_vec)
    attributes(result) <- attributes(seq)
    return(result)
  })
  
  return(result)
}
