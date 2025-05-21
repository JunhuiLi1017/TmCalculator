#' Calculate the melting temperature using the 'Wallace rule'
#' 
#' The Wallace rule is often used as rule of thumb for approximate melting temperature calculations for primers with 14 to 20 nt length.
#'  
#' @param input_seq Sequence (5' to 3') of one strand of the DNA nucleic acid duplex
#' as string or vector of characters (\strong{Note:} Non-DNA characters are ignored
#' by this method).
#'    
#' @param ambiguous Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE. 
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
#' input_seq = c('acgtTGCAATGCCGTAWSDBSY') #for wallace rule
#'
#' out <- Tm_Wallace(input_seq,ambiguous = TRUE)
#' out
#' out$Options
#' 
#' @export Tm_Wallace

Tm_Wallace <- function (input_seq, ambiguous = FALSE){
  my_seq <- check_filter(input_seq, method = "Tm_Wallace")
  n_seq <- length(my_seq)
  n_gc <- n_seq * GC(my_seq, ambiguous = ambiguous)/100
  n_at <- n_seq - n_gc
  tm <- 4 * n_gc + 2 * n_at
  result_list <- vector('list',2L)
  names(result_list) <- c("Tm","Options")
  result_list$Tm <- tm
  result_list$Options <- list("Sequence"=input_seq,"Ambiguous"=ambiguous,"Check filter"=c2s(my_seq),Method="Thein & Wallace 1986")
  class(result_list) <- c("TmCalculator","list")
  attr(result_list, "nonhidden") <- "Tm"
  return(result_list)
}
