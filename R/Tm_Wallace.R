#' Calculate the melting temperature using the 'Wallace rule'
#' 
#' The Wallace rule is often used as rule of thumb for approximate melting temperature calculations for primers with 14 to 20 nt length.
#'
#' @section Length of validity:
#'
#' The 2 + 4 rule was calibrated on 14 to 20 nt hybridisation probes and
#' carries no length-dependent, salt-dependent or concentration-dependent
#' term. Its error therefore grows without bound as sequences lengthen, and
#' on genome-scale windows it returns a number that is not a melting
#' temperature in any useful sense: a 200 bp window is reported at several
#' hundred degrees Celsius.
#'
#' Sequences longer than 30 nt raise a warning that names
#' \code{\link{tm_nn}} as the appropriate alternative. A warning rather than
#' an error, because the rule stays a legitimate rule of thumb and users who
#' knowingly apply it outside its calibrated range should not be blocked;
#' wrap the call in \code{suppressWarnings()} in that case.
#'
#' @param gr_seq Pre-processed sequence(s) in 5' to 3' direction. This should be the output from
#'   to_genomic_ranges() function.
#'    
#' @param ambiguous Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE.
#'
#' @returns Returns a list of sequences with updated Tm attributes
#' 
#' @export
#' @encoding UTF-8
#'
#' @references
#'
#' Thein S L , Lynch J R , Weatherall D J , et al. DIRECT DETECTION OF HAEMOGLOBIN E WITH SYNTHETIC OLIGONUCLEOTIDES[J]. The Lancet, 1986, 327(8472):93.
#'
#' @author
#' 
#' Junhui Li
#' 
#' @examples
#'
#' input_seq = c('acgtTGCAATGCCGTAWSDBSY','acgtTGCCCCGGCCGCGCCGTAWSDBSY') #for wallace rule
#' gr_seq <- to_genomic_ranges(input_seq)
#' out <- tm_wallace(gr_seq, ambiguous = TRUE)
#' out
#' out$Options
#' 
#' @export tm_wallace

tm_wallace <- function(gr_seq, ambiguous = FALSE) {
  # Filter sequence
  gr_seq$sequence <- check_filter_seq(gr_seq$sequence, method = "tm_wallace")

  # Checked here rather than in .tm_wallace_chunk() so that it fires once
  # per call, on the full input.
  .warn_wallace_length(gr_seq$sequence)

  # Calculate Tm for all sequences in one vectorised pass
  chunk_res <- .tm_wallace_chunk(
    list(sequence = as.character(gr_seq$sequence)),
    ambiguous = ambiguous
  )

  gr_seq$GC <- chunk_res$GC
  gr_seq$Tm <- chunk_res$Tm
  gr_seq <- .normalize_tm_gc_metadata(gr_seq)

  # Create result list with proper structure
  # (result$df is computed lazily via `$.TmCalculator`)
  result_list <- list(
    gr = gr_seq,
    options = list(
      Ambiguous = ambiguous,
      Method = "tm_wallace (Thein & Wallace 1986)"
    )
  )

  # Set class and attributes
  class(result_list) <- c("TmCalculator", "list")
  attr(result_list, "nonhidden") <- "gr"

  return(result_list)
}

# -- Length guard -------------------------------------------------------------
# The Wallace rule has no length term, so there is no length at which it
# degrades gracefully; 30 nt is the point past which the rule is no longer
# defensible even as an approximation, being 50% beyond the upper end of the
# range it was calibrated on.
#
# The message reports the count and the longest sequence rather than listing
# offenders: the genome-wide case has millions of windows and enumerating
# them would be unusable.
#' @keywords internal
.warn_wallace_length <- function(seqs, limit = 30L) {
  n <- nchar(as.character(seqs))
  n <- n[!is.na(n)]
  bad <- sum(n > limit)
  if (bad == 0L) return(invisible(FALSE))

  warning(sprintf(
    paste0("tm_wallace(): %d of %d sequences exceed %d nt (longest %d nt). ",
           "The Wallace rule was calibrated on 14-20 nt oligonucleotides and ",
           "has no length-dependent term, so these values are not meaningful ",
           "melting temperatures. Use tm_nn() for sequences of this length."),
    bad, length(n), limit, max(n)), call. = FALSE)
  invisible(TRUE)
}

# -- Wallace-rule Tm over a block of sequences --------------------------------
# `chunk` is list(sequence=) for the whole input.
#' @keywords internal
.tm_wallace_chunk <- function(chunk, ambiguous) {
  seqs <- chunk$sequence
  m    <- length(seqs)
  if (m == 0L) return(list(Tm = numeric(0), GC = numeric(0)))

  # Was a per-sequence loop calling s2c() twice (once for the length, once
  # inside gc_content()) and scanning the character vector five times. Counting now
  # happens once per sequence in compiled code via .gc_vec(); the arithmetic
  # below is unchanged, including the use of the full sequence length rather
  # than the A+C+G+T count when converting the GC percentage back to a base
  # count, so results are identical.
  n_seq <- nchar(seqs)                       # == length(s2c(x))
  pt_gc <- .gc_vec(seqs, ambiguous = ambiguous)

  n_gc <- n_seq * pt_gc / 100
  n_at <- n_seq - n_gc

  list(Tm = as.numeric(4 * n_gc + 2 * n_at), GC = as.numeric(pt_gc))
}
