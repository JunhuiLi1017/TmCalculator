#' Calculate the melting temperature using the 'Wallace rule'
#' 
#' The Wallace rule is often used as rule of thumb for approximate melting temperature calculations for primers with 14 to 20 nt length.
#'  
#' @param ntseq Sequence (5' to 3') of one strand of the DNA nucleic acid duplex
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
#' ntseq = c('acgtTGCAATGCCGTAWSDBSY') #for wallace rule
#'
#' out <- Tm_Wallace(ntseq,ambiguous = TRUE)
#' out
#' out$Options
#' 
#' @export Tm_Wallace

Tm_Wallace <- function (ntseq, ambiguous = FALSE){
  mySeq <- check_filter(ntseq, method = "Tm_Wallace")
  nSeq <- length(mySeq)
  nGC <- nSeq * GC(mySeq, ambiguous = ambiguous)/100
  nAT <- nSeq - nGC
  Tm <- 4 * nGC + 2 * nAT
  resultList <- vector('list',2L)
  names(resultList) <- c("Tm","Options")
  resultList$Tm <- Tm
  resultList$Options <- list("Sequence"=ntseq,"Ambiguous"=ambiguous,"Check filter"=c2s(mySeq),Method="Thein & Wallace 1986")
  class(resultList) <- c("TmCalculator","list")
  attr(resultList, "nonhidden") <- "Tm"
  return(resultList)
}
