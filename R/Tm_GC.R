#' Calculate the melting temperature using empirical formulas based on GC content
#' 
#' Calculate the melting temperature using empirical formulas based on GC content with different options
#' 
#' @param input_seq Sequence (5' to 3') of one strand of the nucleic acid duplex as string or vector of characters.
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
#' @param salt_corr Salt correction method should be chosen when provide 'userset'. Options are "Schildkraut2010", "Wetmur1991","SantaLucia1996","SantaLucia1998-1","Owczarzy2004","Owczarzy2008". Note that "SantaLucia1998-2" is not available for this function.
#'   
#' @param mismatch If 'True' (default) every 'X' in the sequence is counted as mismatch
#' 
#' @param DMSO Percent DMSO concentration in the reaction mixture. Default is 0.
#' 
#' @param formamide_value_unit A list containing formamide concentration value and unit. The list should have two elements:
#'   - value: numeric value of formamide concentration
#'   - unit: character string specifying the unit, either "percent" or "molar"
#' 
#' @param dmso_factor Coefficient of Tm decreases per percent DMSO. Default=0.75 von Ahsen N (2001) <PMID:11673362>. Other published values are 0.5, 0.6 and 0.675.
#' 
#' @param formamide_factor Coefficient of Tm decrease per percent formamide. Default=0.65. Several papers report factors between 0.6 and 0.72.
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
#' input_seq <- c("ATCGTGCGTAGCAGTACGATCAGTAG")
#' out <- Tm_GC(input_seq, ambiguous = TRUE, variant = "Primer3Plus", Na = 50, mismatch = TRUE)
#' out
#' out$tm
#' out$options
#' 
#' @export Tm_GC
Tm_GC <- function(input_seq,
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
                  salt_corr=c("Schildkraut2010",
                             "Wetmur1991",
                             "SantaLucia1996",
                             "SantaLucia1998-1",
                             "Owczarzy2004",
                             "Owczarzy2008"),
                  mismatch=TRUE,
                  DMSO=0,
                  formamide_value_unit=list(value=0, unit="percent"),
                  dmso_factor=0.75,
                  formamide_factor=0.65){
  variant <- match.arg(variant)
  salt_corr <- match.arg(salt_corr)
  my_seq <- check_filter(input_seq,method='Tm_GC')
  n_seq <- length(my_seq)
  pt_gc <- GC(my_seq,ambiguous=ambiguous)
  var_table <- data.frame(A=c(69.3,81.5,81.5,81.5,78.0,67.0,81.5,77.1),
                       B=c(0.41,0.41,0.41,0.41,0.70,0.80,0.41,0.41),
                       C=c(650,675,675,500,500,500,600,528),
                       D=rep(1,8),
                       salt_corr=c(NA,NA,"Schildkraut2010",
                                  rep("Wetmur1991",3),"Schildkraut2010","SantaLucia1998-1"))
  rownames(var_table) <- c("Chester1993","QuikChange","Schildkraut1965","Wetmur1991_MELTING","Wetmur1991_RNA","Wetmur1991_RNA/DNA","Primer3Plus","vonAhsen2001")
  if(is.null(userset)){
    if(!variant %in% rownames(var_table)){
      stop("only Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001 are allowed in variant")
    }else{
      gc_coef <- var_table[variant,]
      salt_corr <- var_table[variant,"salt_corr"]
    }
  }else{
    gc_coef <- as.numeric(userset)
    salt_corr <- salt_corr
  }

  tm = gc_coef[1]+gc_coef[2]*pt_gc-gc_coef[3]/n_seq
  if(!is.na(salt_corr)){
    corr_salt <- salt_correction(Na=Na,K=K,Tris=Tris,Mg=Mg,dNTPs=dNTPs,method=salt_corr,input_seq=my_seq,ambiguous = ambiguous)
    tm <- tm+corr_salt
  }
  if(mismatch == TRUE){
    tm <- tm-gc_coef[4]*(sum(my_seq %in% 'X')*100/n_seq)
  }
  
  corr_chem <- chem_correction(DMSO=DMSO,formamide_value_unit=formamide_value_unit,dmso_factor=dmso_factor,formamide_factor=formamide_factor,pt_gc=pt_gc)
  tm <- tm + corr_chem
  
  result_list <- vector('list',2L)
  names(result_list) <- c("Tm","Options")
  result_list$Tm <- as.numeric(tm)
  result_list$Options <- list("Sequence"=input_seq,"Check filter"=c2s(my_seq),"Variant"=variant,"Na"=Na,
                            "K"=K,"Tris"=Tris,"Mg"=Mg,"dNTPs"=dNTPs,
                            "Salt correlation"=salt_corr,"Ambiguous"=ambiguous,
                            "Mismatch"=mismatch,"Percent of DMSO"=DMSO,
                            "Formamide concentration"=formamide_value_unit$value,
                            "Method for formamide concentration"=formamide_value_unit$unit,
                            "Coeffecient of Tm decreases per percent DMSO"=dmso_factor,
                            "Coefficient of Tm decrease per percent formamide"=formamide_factor)
  class(result_list) <- c("TmCalculator","list")
  attr(result_list, "nonhidden") <- "Tm"
  return(result_list)
}

