#' convert a string into a vector of characters
#' 
#' Simply convert a single string such as "HelloWorld" into a vector of characters such as c("H","e","l","l","o","W","o","r","l","d")
#'
#' @param strings A single string such as "HelloWorld" 
#' 
#' @returns Retrun a vector of characters
#' 
#' @author Junhui Li
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 
#' @examples
#' 
#' s2c(c("HelloWorld"))
#' 
#' @export s2c
s2c <- function(strings){
  vecChar <- vector()
  for (i in nchar(strings):1){
    vecChar <- append(unlist(strsplit(strings,''))[i],vecChar)
  }
  return(vecChar)
}
