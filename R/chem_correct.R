#' Corrections of melting temperature with chemical substances
#' 
#' Apply corrections to melting temperature calculations based on the presence of DMSO and formamide.
#' These corrections are rough approximations and should be used with caution.
#' 
#' @param DMSO Percent DMSO concentration in the reaction mixture. Default: 0
#'   DMSO can lower the melting temperature of nucleic acid duplexes.
#' 
#' @param formamide_value_unit A list containing formamide concentration information:
#'   - value: numeric value of formamide concentration
#'   - unit: character string specifying the unit ("percent" or "molar")
#'   Default: list(value=0, unit="percent")
#' 
#' @param dmso_factor Coefficient of melting temperature (Tm) decrease per percent DMSO.
#'   Default: 0.75 (von Ahsen N, 2001, PMID:11673362)
#'   Other published values: 0.5, 0.6, 0.675
#' 
#' @param formamide_factor Coefficient of melting temperature (Tm) decrease per percent formamide.
#'   Default: 0.65
#'   Literature reports values ranging from 0.6 to 0.72
#' 
#' @param pt_gc Percentage of GC content in the sequence (0-100%)
#'   This is used in molar formamide corrections.
#' 
#' @details 
#' 
#' When formamide_value_unit$unit = "percent":
#' Correction = - factor * percentage_of_formamide
#' 
#' When formamide_value_unit$unit = "molar":
#' Correction = (0.453 * GC/100 - 2.88) * formamide
#' 
#' @references 
#' 
#' von Ahsen N, Wittwer CT, Schutz E, et al. Oligonucleotide melting temperatures under PCR conditions: 
#' deoxynucleotide Triphosphate and Dimethyl sulfoxide concentrations with comparison to alternative 
#' empirical formulas. Clin Chem 2001, 47:1956-1961.
#' 
#' @author Junhui Li
#' 
#' @examples
#' # DMSO correction
#' chem_correct(DMSO = 3)
#' 
#' # Formamide correction (percent)
#' chem_correct(formamide_value_unit = list(value = 1.25, unit = "percent"), pt_gc = 50)
#' 
#' # Formamide correction (molar)
#' chem_correct(formamide_value_unit = list(value = 1.25, unit = "molar"), pt_gc = 50)
#' 
#' @export chem_correct

chem_correct <- function(DMSO = 0,
                           formamide_value_unit = list(value = 0, unit = "percent"),
                           dmso_factor = 0.75,
                           formamide_factor = 0.65,
                           pt_gc = NULL){
  if(DMSO < 0 | formamide_value_unit$value < 0){
    stop("all parameters 'DMSO','formamide_value_unit$value' should not be less than 0")
  }
  
  if(!any(dmso_factor %in% c(0.75,0.5,0.6,0.65,0.675))){
    stop("'dmso_factor' shoule be one of 0.5,0.6,0.65,0.675,0.75")
  }
  if(!any(formamide_factor %in% c(0.65,0.6,0.72))){
    stop("'formamide_factor' shoule be one of 0.6,0.65,0.72")
  }
  if(!formamide_value_unit$unit %in% c("percent", "molar")){
    stop("formamide_value_unit$unit must be either 'percent' or 'molar'")
  }

  corr <- 0
  
  if(DMSO > 0){
    corr <- corr - dmso_factor*DMSO
  }
  
  if(formamide_value_unit$value > 0){
    if(formamide_value_unit$unit == "percent"){
      corr <- corr - formamide_factor*formamide_value_unit$value
    }else if(formamide_value_unit$unit == "molar"){
      if(is.null(pt_gc)){
        stop("'pt_gc' should not be NULL when formamide_value_unit$unit = 'molar'")
      }
      corr <- corr + (0.453*(pt_gc/100)-2.88)*formamide_value_unit$value
    }
  }
  
  return(corr)
}
