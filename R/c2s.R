#' convert a vector of characters into a string
#' 
#' Simply convert a vector of characters such as c("H","e","l","l","o","W","o","r","l","d") into a single string "HelloWorld".
#'
#' @param characters A vector of characters 
#' 
#' @returns Retrun a strings
#' 
#' @author Junhui Li
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 
#' @examples
#' 
#' c2s(c("H","e","l","l","o","W","o","r","l","d"))
#' 
#' @export c2s

c2s <- function(characters){
  strings <- paste0(characters,collapse = "")
  return(strings)
}
