#' Calculate the melting temperature using the 'Wallace rule'
#' 
#' The Wallace rule is often used as rule of thumb for approximate melting temperature calculations for primers with 14 to 20 nt length.
#'  
#' @param gr_seq Pre-processed sequence(s) in 5' to 3' direction. This should be the output from
#'   to_genomic_ranges() function.
#'    
#' @param ambiguous Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE.
#'
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   specifying the parallel backend, e.g.
#'   \code{BiocParallel::MulticoreParam(4)} (Unix/macOS) or
#'   \code{BiocParallel::SnowParam(4)} (all platforms). The default,
#'   \code{BiocParallel::SerialParam()}, runs serially.
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

tm_wallace <- function(gr_seq, ambiguous = FALSE,
                       BPPARAM = BiocParallel::SerialParam()) {
  # Filter sequence
  gr_seq$sequence <- check_filter_seq(gr_seq$sequence, method = "tm_wallace")
  # Calculate Tm for each sequence (chunked, optionally in parallel)
  all_seqs <- as.character(gr_seq$sequence)
  chunk_res <- .bp_map_chunks(
    n = length(gr_seq),
    make_chunk = function(idx) list(sequence = all_seqs[idx]),
    worker = .tm_wallace_chunk,
    BPPARAM = BPPARAM,
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

# -- Chunk worker: Wallace-rule Tm over a block of sequences ------------------
# Called by .bp_map_chunks(), either directly (serial) or on a BiocParallel
# worker. `chunk` is list(sequence=) for this worker's block.
#' @keywords internal
.tm_wallace_chunk <- function(chunk, ambiguous) {
  seqs <- chunk$sequence
  m    <- length(seqs)
  if (m == 0L) return(list(Tm = numeric(0), GC = numeric(0)))

  # Was a per-sequence loop calling s2c() twice (once for the length, once
  # inside gc()) and scanning the character vector five times. Counting now
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
