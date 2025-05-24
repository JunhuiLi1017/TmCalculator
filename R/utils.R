#' Filter and replace invalid bases in nucleotide sequences
#' 
#' This function processes nucleotide sequences by converting characters to uppercase and replacing invalid bases with "." 
#' based on the specified method. The function preserves the sequence length and attributes (name and Tm) of each sequence.
#'
#' @param seq_list Input sequence in 5' to 3' direction. Must be provided as:
#'   - A list of sequences with attributes (name and Tm)
#'   
#' @param method Method to determine valid bases:
#' 
#' TM_Wallace: Valid bases are "A","B","C","D","G","H","I","K","M","N","R","S","T","V","W" and "Y"
#' 
#' TM_GC: Valid bases are "A","B","C","D","G","H","I","K","M","N","R","S","T","V","W", "X" and "Y"
#' 
#' TM_NN: Valid bases are "A","C","G","I" and "T" 
#' 
#' @returns Returns a list of sequences with the same structure as input, where invalid bases are replaced with "."
#' 
#' @author Junhui Li
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 

filter_seq <- function(seq_list, method) {
  if(is.null(seq_list) || length(seq_list) == 0) {
    stop("Input sequence list cannot be NULL or empty")
  }
  
  # Process each sequence in the list
  result <- lapply(seq_list, function(i) {
    if(is.null(i)) {
      return(NULL)
    }
    
    # Convert to uppercase
    i_s2c <- s2c(toupper(i))
    
    # Filter based on method
    if(method == "tm_wallace") {
      baseset <- c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","Y")
    } else if(method == "tm_nn") {
      baseset = c('A','C','G','T','I')
    } else if(method == "tm_gc") {
      baseset <- c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","X","Y")
    } else {
      stop("Invalid method specified")
    }
    
    filtered <- NULL
    for(idx in i_s2c){
      if(idx %in% baseset){
        filtered <- append(filtered,idx)
      }
    }
    
    # Create result with attributes
    attr(i, "filtered_seq") <- paste(filtered, collapse = "")
    
    return(i)
  })
  
  # Remove NULL results
  result <- result[!sapply(result, is.null)]
  
  if(length(result) == 0) {
    stop("No valid sequences after filtering")
  }
  
  return(result)
}

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
#' @keywords internal

c2s <- function(characters){
  strings <- paste0(characters,collapse = "")
  return(strings)
}

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
#' @keywords internal

s2c <- function(strings){
  vec_char <- vector()
  for (i in nchar(strings):1){
    vec_char <- append(unlist(strsplit(strings,''))[i],vec_char)
  }
  return(vec_char)
}

