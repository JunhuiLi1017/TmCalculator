#' Calculate the melting temperature using the 'Wallace rule'
#' 
#' The Wallace rule is often used as rule of thumb for approximate melting temperature calculations for primers with 14 to 20 nt length.
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
  # Calculate Tm for each sequence in the list
  seq_tm <- sapply(seq_along(gr_seq), function(i) {
    filter_seq <- gr_seq$sequence[i]
    n_seq <- length(s2c(filter_seq))
    n_gc <- n_seq * gc(filter_seq, ambiguous = ambiguous) / 100
    n_at <- n_seq - n_gc
    tm <- 4 * n_gc + 2 * n_at
    return(tm)
  })
  gr_seq$Tm <- seq_tm

  # Create result list with proper structure
  result_list <- list(
    Tm = gr_seq,
    Options = list(
      Ambiguous = ambiguous,
      Method = "tm_wallace (Thein & Wallace 1986)"
    )
  )

  # Set class and attributes
  class(result_list) <- c("TmCalculator", "list")
  attr(result_list, "nonhidden") <- "Tm"

  return(result_list)
}
