#' Calculate the melting temperature using empirical formulas based on GC content
#' 
#' Calculate the melting temperature using empirical formulas based on GC content with different options.
#' The function returns a list of sequences with updated Tm attributes and calculation options.
#' 
#' @param raw_seq Pre-processed sequence(s) in 5' to 3' direction. This should be the output from
#'   process_seq() function.
#' 
#' @param ambiguous Logical. If TRUE, ambiguous bases are taken into account when computing the G and C content.
#'   The function handles various ambiguous bases (S, W, M, K, R, Y, V, H, D, B) by proportionally
#'   distributing their contribution to GC content based on their possible nucleotide compositions.
#' 
#' @param userset A vector of four coefficient values. Usersets override value sets.
#' 
#' @param variant Empirical constants coefficient with 8 variants:
#'   - Chester1993: Tm = 69.3 + 0.41(Percentage_GC) - 650/N
#'   - QuikChange: Tm = 81.5 + 0.41(Percentage_GC) - 675/N - Percentage_mismatch
#'   - Schildkraut1965: Tm = 81.5 + 0.41(Percentage_GC) - 675/N + 16.6 x log[Na+]
#'   - Wetmur1991_MELTING: Tm = 81.5 + 0.41(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#'   - Wetmur1991_RNA: Tm = 78 + 0.7(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#'   - Wetmur1991_RNA/DNA: Tm = 67 + 0.8(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#'   - Primer3Plus: Tm = 81.5 + 0.41(Percentage_GC) - 600/N + 16.6 x log[Na+]
#'   - vonAhsen2001: Tm = 77.1 + 0.41(Percentage_GC) - 528/N + 11.7 x log[Na+]
#' 
#' @param Na Millimolar concentration of sodium ions. Default: 50
#' 
#' @param K Millimolar concentration of potassium ions. Default: 0
#' 
#' @param Tris Millimolar concentration of Tris buffer. Default: 0
#' 
#' @param Mg Millimolar concentration of magnesium ions. Default: 0
#' 
#' @param dNTPs Millimolar concentration of deoxynucleotide triphosphates. Default: 0
#' 
#' @param salt_corr_method Salt correction method. Options are:
#'   - "Schildkraut2010": Schildkraut & Lifson 1965
#'   - "Wetmur1991": Wetmur 1991
#'   - "SantaLucia1996": SantaLucia 1996
#'   - "SantaLucia1998-1": SantaLucia 1998 (Method 1)
#'   - "Owczarzy2004": Owczarzy 2004
#'   - "Owczarzy2008": Owczarzy 2008
#'   Note: "SantaLucia1998-2" is not available for this function.
#'   
#' @param mismatch Logical. If TRUE (default), every '.' in the sequence is counted as a mismatch
#' 
#' @param DMSO Percent DMSO concentration in the reaction mixture. Default: 0
#' 
#' @param formamide_value_unit List containing formamide concentration value and unit. Default: list(value = 0, unit = "percent")
#'   - value: Numeric value of formamide concentration
#'   - unit: Either "percent" or "molar"
#' 
#' @param dmso_factor Coefficient of Tm decreases per percent DMSO. Default: 0.75 (von Ahsen et al. 2001)
#'   Other published values are 0.5, 0.6 and 0.675.
#' 
#' @param formamide_factor Coefficient of Tm decrease per percent formamide. Default: 0.65
#'   Several papers report factors between 0.6 and 0.72.
#' 
#' @returns Returns a list with two components:
#'   - Tm: A list of sequences with updated Tm attributes
#'   - Options: A list containing calculation parameters and method information
#' 
#' @references 
#' 
#' Marmur J, Doty P. Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. Journal of Molecular Biology, 1962, 5(1):109-118.
#' 
#' Schildkraut C. Dependence of the melting temperature of DNA on salt concentration. Biopolymers, 2010, 3(2):195-208.
#' 
#' Wetmur JG. DNA Probes: Applications of the Principles of Nucleic Acid Hybridization. CRC Critical Reviews in Biochemistry, 1991, 26(3-4):33.
#' 
#' Untergasser A, Cutcutache I, Koressaar T, et al. Primer3--new capabilities and interfaces. Nucleic Acids Research, 2012, 40(15):e115-e115.
#' 
#' von Ahsen N, Wittwer CT, Schutz E, et al. Oligonucleotide melting temperatures under PCR conditions: deoxynucleotide Triphosphate and Dimethyl sulfoxide concentrations with comparison to alternative empirical formulas. Clin Chem 2001, 47:1956-1961.
#' 
#' @author Junhui Li
#' 
#' @examples
#' 
#' # Example with multiple sequences
#' input_seq <- c("ATCGTGCGTAGCAGTACGATCAGTAG", "ATCGTGCGTAGCAGTACGATCAGTAG")
#' raw_seq <- process_seq(input_seq)
#' out <- tm_gc(raw_seq, ambiguous = TRUE, variant = "Primer3Plus", Na = 50, mismatch = TRUE)
#' out
#' out$Options
#' 
#' @export tm_gc
tm_gc <- function(raw_seq,
                  ambiguous = FALSE,
                  userset = NULL,
                  variant = c("Primer3Plus",
                            "Chester1993",
                            "QuikChange",
                            "Schildkraut1965",
                            "Wetmur1991_MELTING",
                            "Wetmur1991_RNA",
                            "Wetmur1991_RNA/DNA",
                            "vonAhsen2001"),
                  Na = 50,
                  K = 0,
                  Tris = 0,
                  Mg = 0,
                  dNTPs = 0,
                  salt_corr_method = c("Schildkraut2010",
                                    "Wetmur1991",
                                    "SantaLucia1996",
                                    "SantaLucia1998-1",
                                    "Owczarzy2004",
                                    "Owczarzy2008"),
                  mismatch = TRUE,
                  DMSO = 0,
                  formamide_value_unit = list(value = 0, unit = "percent"),
                  dmso_factor = 0.75,
                  formamide_factor = 0.65) {
  variant <- match.arg(variant)
  salt_corr_method <- match.arg(salt_corr_method)
  
  if (is.null(userset)) {
    if (!variant %in% rownames(thermodynamic_gc_params)) {
      stop("only Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001 are allowed in variant")
    } else {
      gc_coef <- thermodynamic_gc_params[variant,]
      salt_corr_method <- thermodynamic_gc_params[variant,"salt_correction"]
    }
  } else {
    gc_coef <- as.numeric(userset)
    salt_corr_method <- salt_corr_method
  }

  # Filter sequence
  seq_checked <- filter_seq(raw_seq, method = 'tm_gc')

  # Calculate Tm for each sequence and update seq_checked
  seq_tm <- lapply(seq_checked, function(x) {
    filtered_seq <- attr(x, "filtered_seq")
    n_seq <- nchar(filtered_seq)
    pt_gc <- gc(filtered_seq, ambiguous = ambiguous)
    tm <- gc_coef[1] + gc_coef[2]*pt_gc - gc_coef[3]/n_seq
    if (mismatch == TRUE) {
      mismatch_count <- sum(filtered_seq %in% 'X')
      tm <- tm - gc_coef[4]*(mismatch_count*100/n_seq)
    }
    if (!is.null(salt_corr_method)) {
      corr_salt <- salt_correct(Na = Na, 
                                   K = K, 
                                   Tris = Tris, 
                                   Mg = Mg, 
                                   dNTPs = dNTPs, 
                                   method = salt_corr_method, 
                                   input_seq = filtered_seq, 
                                   ambiguous = ambiguous)
      tm <- tm + corr_salt
    }
    if (!is.na(DMSO)) {
      corr_chem <- chem_correct(DMSO = DMSO, 
                                   formamide_value_unit = formamide_value_unit, 
                                   dmso_factor = dmso_factor, 
                                   formamide_factor = formamide_factor, 
                                   pt_gc = pt_gc)
      tm <- tm + corr_chem
    }
    attr(x, "Tm") <- as.numeric(tm)
    x
  })
  
  # Create result list with proper structure
  result_list <- list(
    Tm = seq_tm,
    Options = list(
      Ambiguous = ambiguous,
      Method = paste0(variant, " (", 
                     if (variant == "Chester1993") "Chester & Marshak 1993" else
                     if (variant == "QuikChange") "QuikChange Site-Directed Mutagenesis" else
                     if (variant == "Schildkraut1965") "Schildkraut & Lifson 1965" else
                     if (variant == "Wetmur1991_MELTING") "Wetmur 1991 (MELTING)" else
                     if (variant == "Wetmur1991_RNA") "Wetmur 1991 (RNA)" else
                     if (variant == "Wetmur1991_RNA/DNA") "Wetmur 1991 (RNA/DNA)" else
                     if (variant == "Primer3Plus") "Primer3Plus" else
                     "von Ahsen et al. 2001", ")"),
      Na = Na,
      K = K,
      Tris = Tris,
      Mg = Mg,
      dNTPs = dNTPs,
      "Salt correction" = salt_corr_method,
      Mismatch = mismatch,
      "Percent of DMSO" = DMSO,
      "Formamide concentration" = formamide_value_unit$value,
      "Method for formamide concentration" = formamide_value_unit$unit,
      "Coeffecient of Tm decreases per percent DMSO" = dmso_factor,
      "Coefficient of Tm decrease per percent formamide" = formamide_factor
    )
  )
  
  # Set class and attributes
  class(result_list) <- c("TmCalculator", "list")
  attr(result_list, "nonhidden") <- "Tm"
  
  return(result_list)
}
