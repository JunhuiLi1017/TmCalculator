#' Corrections of melting temperature with chemical substances
#' 
#' Corrections coefficient of melting temperature with DMSO and formamide and these corrections are rough approximations.
#' 
#' @param DMSO Percent DMSO
#' 
#' @param fmd Formamide concentration in percentage (fmdmethod="concentration") or molar (fmdmethod="molar").
#' 
#' @param DMSOfactor Coefficient of Tm decreases per percent DMSO. Default=0.75 von Ahsen N (2001) <PMID:11673362>. Other published values are 0.5, 0.6 and 0.675.
#' 
#' @param fmdmethod "concentration" method for formamide concentration in percentage and "molar" for formamide concentration in molar
#' 
#' @param fmdfactor Coefficient of Tm decrease per percent formamide. Default=0.65. Several papers report factors between 0.6 and 0.72.
#' 
#' @param ptGC Percentage of GC(\%).
#' 
#' @details 
#' 
#' fmdmethod = "concentration"
#' 
#' Correction = - factor*percentage_of_formamide
#' 
#' fmdmethod = "molar"
#' 
#' Correction = (0.453*GC/100 - 2.88) x formamide
#' 
#' @references 
#' 
#' von Ahsen N, Wittwer CT, Schutz E , et al. Oligonucleotide melting temperatures under PCR conditions: deoxynucleotide Triphosphate and Dimethyl sulfoxide concentrations with comparison to alternative empirical formulas. Clin Chem 2001, 47:1956-C1961.
#' 
#' @author Junhui Li
#' 
#' @examples
#' chem_correction(DMSO=3)
#' chem_correction(fmd=1.25, fmdmethod="molar", ptGC=50)
#' 
#' @export chem_correction

chem_correction <-function(DMSO=0,
                           fmd=0, 
                           DMSOfactor=0.75,
                           fmdmethod=c("concentration","molar"),
                           fmdfactor=0.65,
                           ptGC){
  #DMSOfactor <- match.arg(DMSOfactor)
  #fmdfactor <- match.arg(fmdfactor)
  if(!any(DMSOfactor %in% c(0.75,0.5,0.6,0.65,0.675))){
    stop("'DMSOfactor' shoule be one of 0.5,0.6,0.65,0.675,0.75")
  }
  if(!any(fmdfactor %in% c(0.65,0.6,0.72))){
    stop("'fmdfactor' shoule be one of 0.6,0.65,0.72")
  }
  fmdmethod <- match.arg(fmdmethod)

  corr <- 0
  ## for DMSO correction
  corr <- corr - DMSOfactor*DMSO
  ## for fmd correction
  if(fmdmethod == "concentration"){
    corr <- corr - fmdfactor*fmd
  }else if(fmdmethod == "molar"){
    if(is.null(ptGC)){
      stop("'ptGC' should not be NULL when fmdmethod = molar")
    }
    corr <- corr + (0.453*(ptGC/100)-2.88)*fmd
  }
  return(corr)
}
