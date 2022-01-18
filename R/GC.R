#' Calculate G and C content of nucleotide sequences
#' 
#' Calculate G and C content of nucleotide sequences. The number of G and C in sequence is divided by length of sequence(when totalnt is TRUE) or the number of all A,T,C,G and ambiguous base.
#' 
#' @param ntseq Sequence (5' to 3') of one strand of the nucleic acid duplex as string or vector of characters.
#' 
#' @param ambiguous Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE.
#' 
#' @param totalnt Sum of 'G' and 'C' bases divided by the length of the sequence when totalnt is TRUE.
#' 
#' @returns Content of G and C(range from 0 to 100%)
#' 
#' @examples 
#' 
#' GC(c("a","t","c","t","g","g","g","c","c","a","g","t","a"))#53.84615
#' GC("GCATSWSYK",ambiguous = TRUE)#55.55556
#' 
#' @author Junhui Li
#' 
#' @export GC
GC <- function (ntseq, ambiguous = FALSE, totalnt = TRUE){
  if (length(ntseq) == 1 && is.na(ntseq)){
    return(NA)
  }
  if(class(ntseq)=="character"){
    ntseq <- toupper(ntseq)
    if(length(ntseq)==1 && nchar(ntseq)>1){
      vecSeq <- s2c(ntseq)
    }else if(length(ntseq) > 1){
      vecSeq <- ntseq
    }
  }else{
    stop("sequence is not characters")
  }
  nSeq <- length(vecSeq)
  if(!all(vecSeq %in% c("A","B","C","D","G","H","I","K","M","N","R","S","T","V","W","Y"))){
    warning("None Nucleic Acid Base are found in input Sequence")
  }
  nc <- sum(vecSeq %in% "C")
  ng <- sum(vecSeq %in% "G")
  na <- sum(vecSeq %in% "A")
  nt <- sum(vecSeq %in% "T")
  
  if(ambiguous == FALSE){
    ngc <- ng+nc
    nat <- na+nt
  }else{
    ngc <- ng+nc+sum(vecSeq %in% "S")
    nat <- na+nt+sum(vecSeq %in% "W")
    # for other ambiguous nucleatide acid base
    if (na+nc != 0) {     #M
      nm <- sum(vecSeq %in% "M")
      ngc <- ngc+nm*nc/(na+nc)
      nat <- nat+nm*na/(na+nc)
    }
    if (ng+nt != 0) {     #K
      nk <- sum(vecSeq %in% "K")
      ngc <- ngc+nk*ng/(ng+nt)
      nat <- nat+nk*nt/(ng+nt)
    }
    if (ng+na != 0) {    #R
      nr <- sum(vecSeq %in% "R")
      ngc <- ngc+nr*ng/(ng+na)
      nat <- nat+nr*na/(ng+na)
    }
    if (nc+nt != 0) {    #Y
      ny <- sum(vecSeq %in% "Y")
      ngc <- ngc+ny*nc/(nc+nt)
      nat <- nat+ny*nt/(nc+nt)
    }
    if (na+nc+ng != 0) {    #V
      nv <- sum(vecSeq %in% "V")
      ngc <- ngc+nv*(nc+ng)/(na+nc+ng)
      nat <- nat+nv*na/(na+nc+ng)
    }
    if (na+nc+nt != 0) {    #H
      nh <- sum(vecSeq %in% "H")
      ngc <- ngc+nh*nc/(na+nc+nt)
      nat <- nat+nh*(na+nt)/(na+nc+nt)
    }
    if (na+ng+nt != 0) {    #D
      nd <- sum(vecSeq %in% "D")
      ngc <- ngc+nd*ng/(na+ng+nt)
      nat <- nat+nd*(na+nt)/(na+ng+nt)
    }
    if (nc+ng+nt != 0) {    #B
      nb <- sum(vecSeq %in% "B")
      ngc <- ngc+nb*(nc+ng)/(nc+ng+nt)
      nat <- nat+nb*nt/(nc+ng+nt)
    }
  }
  if(totalnt){
    cat("argument totalnt is deprecated\n")
    ptGC <- 100*(ngc)/nSeq
    return(ptGC)
  }else{
    if (ngc+nat == 0) {
      ptGC <- NA
    }else {
      ptGC <- 100*ngc/(ngc+nat)
    }
  }
  return(ptGC)
}
