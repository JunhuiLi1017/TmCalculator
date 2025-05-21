#' complement and reverse complement base of nucleotide sequences
#' 
#' get reverse complement and complement base of nucleotide sequences
#' 
#' @param input_seq Sequence (5' to 3') of one strand of the nucleic acid duplex
#'   as string or vector of characters
#' 
#' @param reverse Logical value, TRUE is reverse complement sequence, FALSE is not.
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 
#' @author Junhui Li
#' 
#' @examples
#' 
#' complement("ATCGYCGYsWwsaVv")
#' complement("ATCGYCGYsWwsaVv",reverse=TRUE)
#' 
#' @export complement

complement <- function(input_seq,reverse=FALSE){
  complement_table <- matrix(c('A','T','B','V','C','G','D','H','G','C','H','D','M','K','N','N','R','Y','S',
                               'S','T','A','U','A','V','B','W','W','X','X','Y','R','a','t','b','v','c','g',
                               'd','h','g','c','h','d','m','k','n','n','r','y','s','s','t','a','u','a','v',
                               'b','w','w','x','x','y','r'),ncol = 2, byrow = TRUE)
  new_seq <- NULL
  if(length(input_seq)==1){
    line_seq <- s2c(input_seq)
  } else {
    line_seq <- input_seq
  }
  for (i in line_seq){
    if(i %in% complement_table){
      n <- which(complement_table[,1] %in% i)
      new_seq <- append(new_seq,complement_table[n,2])
    }
  }
  if(reverse==TRUE){
    new_seq <- rev(new_seq)
  }
  rc_seq <- c2s(new_seq)
  return(rc_seq)
}
