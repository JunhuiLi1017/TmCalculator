#' Check and filter invalid base of nucleotide sequences
#' 
#' In general, whitespaces and non-base characters are removed and characters are converted to uppercase in given method.
#'
#' @param ntseq Sequence (5' to 3') of one strand of the DNA nucleic acid duplex
#'   as string or vector of characters
#'   
#' @param method 
#' 
#' TM_Wallace: check and return "A","B","C","D","G","H","I","K","M","N","R","S","T","V","W" and "Y"
#' 
#' TM_GC: check and return "A","B","C","D","G","H","I","K","M","N","R","S","T","V","W", "X" and "Y"
#' 
#' TM_NN: check and return "A","C","G","I" and "T" 
#' 
#' @returns Return a sequence which fullfils the requirements of the given method.
#' 
#' @author Junhui Li
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 
#' @examples
#' 
#' ntseq <- c("ATCGBDHKMNRVYWSqq")
#' check_filter(ntseq,method='Tm_Wallace')
#' check_filter(ntseq,method='Tm_NN')
#' 
#' @export check_filter

check_filter <- function(ntseq,method){
  mySeq <- s2c(ntseq)
  mySeq <- toupper(mySeq)
  if (method == 'Tm_Wallace'){
    baseset <- c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","Y")
  }else if (method == 'Tm_GC'){
    baseset <- c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","X","Y")
  }else if(method == 'Tm_NN'){
    baseset = c('A','C','G','T','I')
  }else{
    stop("Only methods 'Tm_Wallace' or 'Tm_GC' or 'Tm_NN' is allowed")
  }
  finalSeq <- NULL
  #i='A'
  for(i in mySeq){
    if(i %in% baseset){
      finalSeq <- append(finalSeq,i)
    }
  }
  return(finalSeq)
}
