#' Calculate G and C content of nucleotide sequences
#' 
#' Calculate G and C content of nucleotide sequences. The function calculates the percentage of G and C bases
#' relative to the total number of A, T, G, and C bases in the sequence.
#' 
#' @param input_seq Sequence (5' to 3') of one strand of the nucleic acid duplex. Can be provided as either:
#'   - A character string (e.g., "ATGCG")
#'   - A path to a FASTA file containing the sequence(s)
#' 
#' @param ambiguous Logical. If TRUE, ambiguous bases are taken into account when computing the G and C content.
#'   The function handles various ambiguous bases (S, W, M, K, R, Y, V, H, D, B) by proportionally
#'   distributing their contribution to GC content based on their possible nucleotide compositions.
#'   For example:
#'   - S (G or C) contributes fully to GC content
#'   - W (A or T) contributes fully to AT content
#'   - M (A or C) contributes proportionally based on the ratio of A to C in the sequence
#'   - And so on for other ambiguous bases
#' 
#' @returns Content of G and C as a percentage (range from 0 to 100%)
#' 
#' @examples 
#' 
#' # Calculate GC content of a DNA sequence
#' gc(c("a","t","c","t","g","g","g","c","c","a","g","t","a"))  # 53.85%
#' 
#' # Calculate GC content including ambiguous bases
#' gc("GCATSWSYK", ambiguous = TRUE)  # 55.56%
#' 
#' @author Junhui Li
#' 
#' @export gc
gc <- function(input_seq, ambiguous = FALSE) {
  if (length(input_seq) == 0) return(NA)
  if (!inherits(input_seq, "character")) {
    stop("sequence must be a character string or vector")
  }
  if (length(input_seq) == 1 && is.na(input_seq)) return(NA)

  # A vector of length > 1 is documented as one sequence supplied as separate
  # characters, e.g. gc(c("a","t","c","g")); collapse it so that a single
  # code path serves both forms. Everything below is .gc_vec(), so gc(),
  # tm_gc(), tm_wallace() and salt_correct() now share one definition of GC
  # and one implementation of the counting.
  if (length(input_seq) > 1) {
    input_seq <- paste0(input_seq, collapse = "")
  }
  .gc_vec(input_seq, ambiguous = ambiguous)
}




#' Vectorized GC percent over a character vector of sequences
#'
#' Mirrors \code{gc()} exactly, including its ambiguity apportioning and its
#' \code{GC/(A+C+G+T)} denominator, but counts every base in a single compiled
#' pass per sequence (\code{cpp_base_counts()}) instead of splitting each
#' sequence into a character vector with \code{s2c()} and scanning it five
#' times. The apportioning arithmetic stays here, vectorised over the count
#' matrix, so the published formula remains the readable one.
#'
#' @param input_seq Character vector of sequences.
#' @param ambiguous Logical; apportion ambiguous IUPAC codes as \code{gc()} does.
#' @return Numeric vector of GC percentages, \code{NA} where a sequence is
#'   \code{NA} or contains no countable base.
#' @keywords internal
.gc_vec <- function(input_seq, ambiguous = FALSE) {
  x <- as.character(input_seq)
  m <- cpp_base_counts(x)          # uppercasing happens in the compiled pass

  # gc() warns once per call; warn once for the whole vector instead.
  if (any(m[, "other"] > 0L, na.rm = TRUE)) {
    warning("Non-nucleic acid bases found in input sequence")
  }

  nA <- m[, "A"]; nC <- m[, "C"]; nG <- m[, "G"]; nT <- m[, "T"]

  if (!isTRUE(ambiguous)) {
    ngc <- nG + nC
    nat <- nA + nT
  } else {
    ngc <- nG + nC + m[, "S"]
    nat <- nA + nT + m[, "W"]
    # gc() skips a code entirely when its denominator is zero, which is the
    # same as adding zero; the mask reproduces that without dividing by 0.
    apportion <- function(ngc, nat, k, denom, gc_part, at_part) {
      ok <- !is.na(denom) & denom != 0
      ngc[ok] <- ngc[ok] + k[ok] * gc_part[ok] / denom[ok]
      nat[ok] <- nat[ok] + k[ok] * at_part[ok] / denom[ok]
      list(ngc = ngc, nat = nat)
    }
    r <- apportion(ngc, nat, m[, "M"], nA + nC,      nC,      nA)      ; ngc <- r$ngc; nat <- r$nat
    r <- apportion(ngc, nat, m[, "K"], nG + nT,      nG,      nT)      ; ngc <- r$ngc; nat <- r$nat
    r <- apportion(ngc, nat, m[, "R"], nG + nA,      nG,      nA)      ; ngc <- r$ngc; nat <- r$nat
    r <- apportion(ngc, nat, m[, "Y"], nC + nT,      nC,      nT)      ; ngc <- r$ngc; nat <- r$nat
    r <- apportion(ngc, nat, m[, "V"], nA + nC + nG, nC + nG, nA)      ; ngc <- r$ngc; nat <- r$nat
    r <- apportion(ngc, nat, m[, "H"], nA + nC + nT, nC,      nA + nT) ; ngc <- r$ngc; nat <- r$nat
    r <- apportion(ngc, nat, m[, "D"], nA + nG + nT, nG,      nA + nT) ; ngc <- r$ngc; nat <- r$nat
    r <- apportion(ngc, nat, m[, "B"], nC + nG + nT, nC + nG, nT)      ; ngc <- r$ngc; nat <- r$nat
  }

  tot <- ngc + nat
  out <- ifelse(!is.na(tot) & tot == 0, NA_real_, 100 * ngc / tot)
  as.numeric(out)
}
