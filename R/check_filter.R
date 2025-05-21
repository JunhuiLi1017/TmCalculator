#' Check and filter invalid base of nucleotide sequences
#' 
#' In general, whitespaces and non-base characters are removed and characters are converted to uppercase in given method.
#'
#' @param input_seq Sequence (5' to 3') of one strand of the DNA nucleic acid duplex
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
#' input_seq <- c("ATCGBDHKMNRVYWSqq")
#' check_filter(input_seq,method='Tm_Wallace')
#' check_filter(input_seq,method='Tm_NN')
#' 
#' @export check_filter

check_filter <- function(input_seq,method){
  my_seq <- s2c(input_seq)
  my_seq <- toupper(my_seq)
  if (method == 'Tm_Wallace'){
    base_set <- c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","Y")
  }else if (method == 'Tm_GC'){
    base_set <- c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","X","Y")
  }else if(method == 'Tm_NN'){
    base_set = c('A','C','G','T','I')
  }else{
    stop("Only methods 'Tm_Wallace' or 'Tm_GC' or 'Tm_NN' is allowed")
  }
  final_seq <- NULL
  #i='A'
  for(i in my_seq){
    if(i %in% base_set){
      final_seq <- append(final_seq,i)
    }
  }
  return(final_seq)
}
