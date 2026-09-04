#' Generate complementary sequence
#' 
#' Generate the complementary sequence of a nucleic acid sequence, with an option to reverse it.
#' 
#' @param input_seq Input sequence(s) in 5' to 3' direction. Must be provided as either:
#'   - A character string (e.g., c("ATGCG", "GCTAG"))
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
#' 
generate_complement <- function(input_seq, reverse = FALSE) {
  # Complement map, identical to the lookup table this function used before:
  #   A<->T  G<->C  M<->K  R<->Y  B<->V  D<->H   W, S, N, I self-complement
  # Case-sensitive, as before: lowercase input matched no table entry and so
  # became ".", and that behaviour is preserved rather than quietly improved.
  from <- "ATGCMKRYWSBVDHNI"
  to   <- "TACGKMYRWSVBHDNI"

  x <- as.character(input_seq)
  x <- gsub(paste0("[^", from, "]"), ".", x)   # unmatched characters -> "."
  comp <- chartr(from, to, x)

  if (reverse) {
    comp <- vapply(strsplit(comp, "", fixed = TRUE),
                   function(z) paste(rev(z), collapse = ""),
                   character(1L), USE.NAMES = FALSE)
  }
  # sapply(USE.NAMES = TRUE) used to name the result after the input strings;
  # kept so that callers relying on those names are unaffected.
  names(comp) <- input_seq
  comp
}
