#' Generate complementary sequence
#' 
#' Generate the complementary sequence of a nucleic acid sequence, with an option to reverse it.
#' 
#' @param input_seq Input sequence(s) in 5' to 3' direction. Must be provided as either:
#'   - A character string (e.g., c("ATGCG", "GCTAG"))
#'
#' @param reverse Logical, controlling which of the two ways of writing the
#'   opposite strand is returned.
#'
#'   \code{FALSE} (default) gives the plain complement: base \code{i} of the
#'   result pairs with base \code{i} of the input, so written underneath the
#'   input it runs \strong{3' to 5'}. This is the form
#'   \code{\link{to_genomic_ranges}} expects for \code{complement_seq}, and the
#'   form its auto-generated complements take.
#'
#'   \code{TRUE} gives the reverse complement: the same strand written the
#'   conventional way round, \strong{5' to 3'}. It is what you would order from
#'   a supplier, and it is \emph{not} what \code{complement_seq} wants --
#'   passing it there pairs every position against the wrong base.
#'
#'   \preformatted{
#'   input                 5'-A T G C G-3'
#'   reverse = FALSE          T A C G C     (3' to 5', pairs position by position)
#'   reverse = TRUE        5'-C G C A T-3'  (the same strand, written 5' to 3')
#'   }
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
#' # Plain complement: pairs position by position, reads 3' to 5'
#' generate_complement("ATGCG", reverse = FALSE)
#'
#' # Reverse complement: the same strand written 5' to 3'
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
