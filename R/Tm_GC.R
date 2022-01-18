#' Calculate the melting temperature using empirical formulas based on GC content
#' 
#' Calculate the melting temperature using empirical formulas based on GC content with different options
#' 
#' @param ntseq Sequence (5' to 3') of one strand of the nucleic acid duplex as string or vector of characters.
#' 
#' @param ambiguous Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE.
#'
#' @param userset A vector of four coefficient values. Usersets override value sets.
#' 
#' @param variant Empirical constants coefficient with 8 variant: Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001
#' 
#' @param Na Millimolar concentration of Na, default is 0
#' 
#' @param K Millimolar concentration of K, default is 0
#' 
#' @param Tris Millimolar concentration of Tris, default is 0
#' 
#' @param Mg Millimolar concentration of Mg, default is 0
#' 
#' @param dNTPs Millimolar concentration of dNTPs, default is 0
#' 
#' @param saltcorr Salt correction method should be chosen when provide 'userset'. Options are "Schildkraut2010", "Wetmur1991","SantaLucia1996","SantaLucia1998-1","Owczarzy2004","Owczarzy2008". Note that "SantaLucia1998-2" is not available for this function.
#'   
#' @param mismatch If 'True' (default) every 'X' in the sequence is counted as mismatch
#' 
#' @param DMSO Percent DMSO
#' 
#' @param fmd Formamide concentration in percentage (fmdmethod="concentration") or molar (fmdmethod="molar").
#' 
#' @param DMSOfactor Coeffecient of Tm decreases per percent DMSO. Default=0.75 von Ahsen N (2001) <PMID:11673362>. Other published values are 0.5, 0.6 and 0.675.
#' 
#' @param fmdfactor Coeffecient of Tm decrease per percent formamide. Default=0.65. Several papers report factors between 0.6 and 0.72.
#' 
#' @param fmdmethod "concentration" method for formamide concentration in percentage and "molar" for formamide concentration in molar
#' 
#' @details 
#' 
#' Empirical constants coefficient with 8 variant:
#' 
#' Chester1993: Tm = 69.3 + 0.41(Percentage_GC) - 650/N
#'   
#' QuikChange: Tm = 81.5 + 0.41(Percentage_GC) - 675/N - Percentage_mismatch
#'   
#' Schildkraut1965: Tm = 81.5 + 0.41(Percentage_GC) - 675/N + 16.6 x log[Na+]
#'   
#' Wetmur1991_MELTING: Tm = 81.5 + 0.41(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#'   
#' Wetmur1991_RNA: Tm = 78 + 0.7(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#'   
#' Wetmur1991_RNA/DNA: Tm = 67 + 0.8(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#'   
#' Primer3Plus: Tm = 81.5 + 0.41(Percentage_GC) - 600/N + 16.6 x log[Na+]
#'   
#' vonAhsen2001: Tm = 77.1 + 0.41(Percentage_GC) - 528/N + 11.7 x log[Na+]
#' 
#' @references 
#' 
#' Marmur J , Doty P . Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature.[J]. Journal of Molecular Biology, 1962, 5(1):109-118.
#' 
#' Schildkraut C . Dependence of the melting temperature of DNA on salt concentration[J]. Biopolymers, 2010, 3(2):195-208.
#' 
#' Wetmur J G . DNA Probes: Applications of the Principles of Nucleic Acid Hybridization[J]. CRC Critical Reviews in Biochemistry, 1991, 26(3-4):33.
#' 
#' Untergasser A , Cutcutache I , Koressaar T , et al. Primer3--new capabilities and interfaces[J]. Nucleic Acids Research, 2012, 40(15):e115-e115.
#' 
#' von Ahsen N, Wittwer CT, Schutz E , et al. Oligonucleotide melting temperatures under PCR conditions: deoxynucleotide Triphosphate and Dimethyl sulfoxide concentrations with comparison to alternative empirical formulas. Clin Chem 2001, 47:1956-1961.
#' 
#' @author Junhui Li
#' 
#' @examples
#' 
#' ntseq <- c("ATCGTGCGTAGCAGTACGATCAGTAG")
#' out <- Tm_GC(ntseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)
#' out
#' out$Tm
#' out$Options
#' 
#' @export Tm_GC
Tm_GC <- function(ntseq,
                  ambiguous=FALSE,
                  userset=NULL,
                  variant=c("Primer3Plus",
                            "Chester1993",
                            "QuikChange",
                            "Schildkraut1965",
                            "Wetmur1991_MELTING",
                            "Wetmur1991_RNA",
                            "Wetmur1991_RNA/DNA",
                            "vonAhsen2001"),
                  Na=0,
                  K=0,
                  Tris=0,
                  Mg=0, 
                  dNTPs=0,
                  saltcorr=c("Schildkraut2010",
                             "Wetmur1991",
                             "SantaLucia1996",
                             "SantaLucia1998-1",
                             "Owczarzy2004",
                             "Owczarzy2008"),
                  mismatch=TRUE,
                  DMSO=0,
                  fmd=0, 
                  DMSOfactor=0.75,
                  fmdfactor=0.65,
                  fmdmethod=c("concentration","molar")){
  variant <- match.arg(variant)
  saltcorr <- match.arg(saltcorr)
  mySeq <- check_filter(ntseq,method='Tm_GC')
  nSeq <- length(mySeq)
  ptGC <- GC(mySeq,ambiguous=ambiguous)
  varTab <- data.frame(A=c(69.3,81.5,81.5,81.5,78.0,67.0,81.5,77.1),
                       B=c(0.41,0.41,0.41,0.41,0.70,0.80,0.41,0.41),
                       C=c(650,675,675,500,500,500,600,528),
                       D=rep(1,8),
                       saltcorr=c(NA,NA,"Schildkraut2010",
                                  rep("Wetmur1991",3),"Schildkraut2010","SantaLucia1998-1"))
  rownames(varTab) <- c("Chester1993","QuikChange","Schildkraut1965","Wetmur1991_MELTING","Wetmur1991_RNA","Wetmur1991_RNA/DNA","Primer3Plus","vonAhsen2001")
  if(is.null(userset)){
    if(!variant %in% rownames(varTab)){
      stop("only Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001 are allowed in variant")
    }else{
      gcCoef <- varTab[variant,]
      saltcorr <- varTab[variant,"saltcorr"]
    }
  }else{
    gcCoef <- as.numeric(userset)
    saltcorr <- saltcorr
  }

  Tm = gcCoef[1]+gcCoef[2]*ptGC-gcCoef[3]/nSeq
  if(!is.na(saltcorr)){
    corrSalt <- salt_correction(Na=Na,K=K,Tris=Tris,Mg=Mg,dNTPs=dNTPs,method=saltcorr,ntseq=mySeq,ambiguous = ambiguous)
    Tm <- Tm+corrSalt
  }
  if(mismatch == TRUE){
    Tm <- Tm-gcCoef[4]*(sum(mySeq %in% 'X')*100/nSeq)
  }
  
  corrChem <- chem_correction(DMSO=DMSO,fmd=fmd,DMSOfactor=DMSOfactor,fmdmethod=fmdmethod,fmdfactor=fmdfactor,ptGC=ptGC)
  Tm <- Tm + corrChem
  
  resultList <- vector('list',2L)
  names(resultList) <- c("Tm","Options")
  resultList$Tm <- as.numeric(Tm)
  resultList$Options <- list("Sequence"=ntseq,"Check filter"=c2s(mySeq),"Variant"=variant,"Na"=Na,
                            "K"=K,"Tris"=Tris,"Mg"=Mg,"dNTPs"=dNTPs,
                            "Salt correlation"=saltcorr,"Ambiguous"=ambiguous,
                            "Mismatch"=mismatch)
  class(resultList) <- c("TmCalculator","list")
  attr(resultList, "nonhidden") <- "Tm"
  return(resultList)
}

