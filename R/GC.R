#' Calculate G and C content of nucleotide sequences
#' 
#' Calculate G and C content of nucleotide sequences. The function calculates the percentage of G and C bases
#' relative to the total number of A, T, G, and C bases in the sequence.
#' 
#' @param input_seq Sequence (5' to 3') of one strand of the nucleic acid duplex. Can be provided as either:
#'   - A character string (e.g., "ATGCG")
#'   - A vector of characters (e.g., c("A","T","G","C","G"))
#' 
#' @param ambiguous Logical. If TRUE, ambiguous bases are taken into account when computing the G and C content.
#'   The function handles various ambiguous bases (S, W, M, K, R, Y, V, H, D, B) by proportionally
#'   distributing their contribution to GC content based on their possible nucleotide compositions.
#'   For example:
#'   - S (G or C) contributes fully to GC content
#'   - W (A or T) contributes fully to AT content
#'   - M (A or C) contributes proportionally based on the ratio of A to C in the sequence
#'   - And so on for other ambiguous bases
#' 
#' @returns Content of G and C as a percentage (range from 0 to 100%)
#' 
#' @examples 
#' 
#' # Calculate GC content of a DNA sequence
#' GC(c("a","t","c","t","g","g","g","c","c","a","g","t","a"))  # 53.85%
#' 
#' # Calculate GC content including ambiguous bases
#' GC("GCATSWSYK", ambiguous = TRUE)  # 55.56%
#' 
#' @author Junhui Li
#' 
#' @export GC
GC <- function(input_seq, ambiguous = FALSE) {
  if (length(input_seq) == 1 && is.na(input_seq)) {
    return(NA)
  }
  
  if (inherits(input_seq, "character")) {
    input_seq <- toupper(input_seq)
    if (length(input_seq) == 1 && nchar(input_seq) > 1) {
      vec_seq <- s2c(input_seq)
    } else if (length(input_seq) > 1) {
      vec_seq <- input_seq
    }
  } else {
    stop("sequence is not characters")
  }
  
  n_seq <- length(vec_seq)
  if (!all(vec_seq %in% c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","Y"))) {
    warning("None Nucleic Acid Base are found in input Sequence")
  }
  
  nc <- sum(vec_seq %in% "C")
  ng <- sum(vec_seq %in% "G")
  na <- sum(vec_seq %in% "A")
  nt <- sum(vec_seq %in% "T")
  
  if (ambiguous == FALSE) {
    ngc <- ng + nc
    nat <- na + nt
  } else {
    ngc <- ng + nc + sum(vec_seq %in% "S")
    nat <- na + nt + sum(vec_seq %in% "W")
    # for other ambiguous nucleatide acid base
    if (na + nc != 0) {     #M
      nm <- sum(vec_seq %in% "M")
      ngc <- ngc + nm * nc/(na + nc)
      nat <- nat + nm * na/(na + nc)
    }
    if (ng + nt != 0) {     #K
      nk <- sum(vec_seq %in% "K")
      ngc <- ngc + nk * ng/(ng + nt)
      nat <- nat + nk * nt/(ng + nt)
    }
    if (ng + na != 0) {    #R
      nr <- sum(vec_seq %in% "R")
      ngc <- ngc + nr * ng/(ng + na)
      nat <- nat + nr * na/(ng + na)
    }
    if (nc + nt != 0) {    #Y
      ny <- sum(vec_seq %in% "Y")
      ngc <- ngc + ny * nc/(nc + nt)
      nat <- nat + ny * nt/(nc + nt)
    }
    if (na + nc + ng != 0) {    #V
      nv <- sum(vec_seq %in% "V")
      ngc <- ngc + nv * (nc + ng)/(na + nc + ng)
      nat <- nat + nv * na/(na + nc + ng)
    }
    if (na + nc + nt != 0) {    #H
      nh <- sum(vec_seq %in% "H")
      ngc <- ngc + nh * nc/(na + nc + nt)
      nat <- nat + nh * (na + nt)/(na + nc + nt)
    }
    if (na + ng + nt != 0) {    #D
      nd <- sum(vec_seq %in% "D")
      ngc <- ngc + nd * ng/(na + ng + nt)
      nat <- nat + nd * (na + nt)/(na + ng + nt)
    }
    if (nc + ng + nt != 0) {    #B
      nb <- sum(vec_seq %in% "B")
      ngc <- ngc + nb * (nc + ng)/(nc + ng + nt)
      nat <- nat + nb * nt/(nc + ng + nt)
    }
  }
  
  if (ngc + nat == 0) {
    pt_gc <- NA
  } else {
    pt_gc <- 100 * ngc/(ngc + nat)
  }
  
  return(pt_gc)
}

