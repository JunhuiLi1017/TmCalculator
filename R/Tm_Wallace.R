#' Calculate the melting temperature using the 'Wallace rule'
#' 
#' The Wallace rule is often used as rule of thumb for approximate melting temperature calculations for primers with 14 to 20 nt length.
#'  
#' @param raw_seq Pre-processed sequence(s) in 5' to 3' direction. This should be the output from
#'   process_seq() function.
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
#' raw_seq <- process_seq(input_seq)
#' out <- tm_wallace(raw_seq, ambiguous = TRUE)
#' out
#' out$Options
#' 
#' @export tm_wallace

tm_wallace <- function(raw_seq, ambiguous = FALSE) {
  # Filter sequence
  seq_checked <- filter_seq(raw_seq, method = "tm_wallace")
  
  # Calculate Tm for each sequence in the list
  seq_tm <- lapply(seq_checked, function(i) {
    filter_seq <- attr(i, "filtered_seq")
    n_seq <- length(s2c(filter_seq))
    n_gc <- n_seq * gc(filter_seq, ambiguous = ambiguous) / 100
    n_at <- n_seq - n_gc
    tm <- 4 * n_gc + 2 * n_at
    
    # Update the Tm attribute
    attr(i, "Tm") <- tm
    return(i)
  })
  
  # Create result list with proper structure
  result_list <- list(
    Tm = seq_tm,
    Options = list(
      Ambiguous = ambiguous,
      Method = "tm_wallace (Thein & Wallace 1986)"
    )
  )
  
  # Set class and attributes
  attr(result_list, "class") <- c("list","TmCalculator")
  attr(result_list, "nonhidden") <- "Tm"
  
  return(result_list)
}
