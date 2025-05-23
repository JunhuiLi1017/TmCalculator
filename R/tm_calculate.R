#' Calculate melting temperature using multiple methods
#' 
#' A wrapper function that calculates melting temperature using multiple methods:
#' - Nearest Neighbor thermodynamics (tm_nn)
#' - GC content-based method (tm_gc)
#' - Wallace rule (tm_wallace)
#' 
#' @section Available Options:
#' 
#' \strong{Method Selection:}
#' \itemize{
#'   \item \code{method}: c("tm_nn", "tm_gc", "tm_wallace")
#' }
#' 
#' \strong{Nearest Neighbor (NN) Method Options:}
#' \itemize{
#'   \item \code{nn_table}: 
#'     \itemize{
#'       \item "DNA_NN_Breslauer_1986"
#'       \item "DNA_NN_Sugimoto_1996"
#'       \item "DNA_NN_Allawi_1998"
#'       \item "DNA_NN_SantaLucia_2004" (default)
#'       \item "RNA_NN_Freier_1986"
#'       \item "RNA_NN_Xia_1998"
#'       \item "RNA_NN_Chen_2012"
#'       \item "RNA_DNA_NN_Sugimoto_1995"
#'     }
#'   \item \code{tmm_table} (Terminal Mismatches):
#'     \itemize{
#'       \item "DNA_TMM_Bommarito_2000" (default)
#'     }
#'   \item \code{imm_table} (Internal Mismatches):
#'     \itemize{
#'       \item "DNA_IMM_Peyret_1999" (default)
#'     }
#'   \item \code{de_table} (Dangling Ends):
#'     \itemize{
#'       \item "DNA_DE_Bommarito_2000" (default)
#'       \item "RNA_DE_Turner_2010"
#'     }
#' }
#' 
#' \strong{GC Method Options:}
#' \itemize{
#'   \item \code{variant}:
#'     \itemize{
#'       \item "Primer3Plus" (default)
#'       \item "Chester1993"
#'       \item "QuikChange"
#'       \item "Schildkraut1965"
#'       \item "Wetmur1991_MELTING"
#'       \item "Wetmur1991_RNA"
#'       \item "Wetmur1991_RNA/DNA"
#'       \item "vonAhsen2001"
#'     }
#' }
#' 
#' \strong{Salt Correction Options:}
#' \itemize{
#'   \item \code{salt_corr_method}:
#'     \itemize{
#'       \item "Schildkraut2010" (default)
#'       \item "Wetmur1991"
#'       \item "SantaLucia1996"
#'       \item "SantaLucia1998-1"
#'       \item "SantaLucia1998-2"
#'       \item "Owczarzy2004"
#'       \item "Owczarzy2008"
#'     }
#' }
#' 
#' \strong{Formamide Unit Options:}
#' \itemize{
#'   \item \code{formamide_value_unit$unit}:
#'     \itemize{
#'       \item "percent" (default)
#'       \item "molar"
#'     }
#' }
#' 
#' \strong{Other Parameters:}
#' \itemize{
#'   \item \code{ambiguous}: TRUE/FALSE (default: FALSE)
#'   \item \code{shift}: Integer value (default: 0)
#'   \item \code{dnac_high}: Numeric value in nM (default: 25)
#'   \item \code{dnac_low}: Numeric value in nM (default: 25)
#'   \item \code{self_comp}: TRUE/FALSE (default: FALSE)
#'   \item \code{Na}: Millimolar concentration (default: 50)
#'   \item \code{K}: Millimolar concentration (default: 0)
#'   \item \code{Tris}: Millimolar concentration (default: 0)
#'   \item \code{Mg}: Millimolar concentration (default: 0)
#'   \item \code{dNTPs}: Millimolar concentration (default: 0)
#'   \item \code{DMSO}: Percent concentration (default: 0)
#'   \item \code{dmso_factor}: Numeric value (default: 0.75)
#'   \item \code{formamide_factor}: Numeric value (default: 0.65)
#'   \item \code{mismatch}: TRUE/FALSE (default: TRUE)
#' }
#' 
#' @param input_seq Input sequence(s) in 5' to 3' direction. Can be provided as either:
#'   - A character string (e.g., "ATGCG")
#'   - A path to a FASTA file containing the sequence(s)
#' 
#' @param complement_seq Complementary sequence(s) in 3' to 5' direction. If not provided,
#'   the function will automatically generate it from input_seq. This is the template/target
#'   sequence that the input sequence will hybridize with.
#'   - A character string (e.g., "ATGCG")
#'   - A path to a FASTA file containing the sequence(s)
#'   - A NULL value (default)
#' 
#' @param method Method(s) to use for Tm calculation. Can be one or more of:
#'   - "tm_nn": Nearest Neighbor thermodynamics (default)
#'   - "tm_gc": GC content-based method
#'   - "tm_wallace": Wallace rule
#'   Default: c("tm_nn", "tm_gc", "tm_wallace")
#' 
#' @param ambiguous Logical. If TRUE, ambiguous bases are taken into account when computing
#'   the G and C content. The function handles various ambiguous bases (S, W, M, K, R, Y, V, H, D, B)
#'   by proportionally distributing their contribution to GC content based on their possible
#'   nucleotide compositions. Default: FALSE
#' 
#' @param shift Integer value controlling the alignment offset between primer and template sequences.
#'   Only applicable for the NN method. Default: 0
#' 
#' @param nn_table Thermodynamic nearest-neighbor parameters for different nucleic acid hybridizations.
#'   Only applicable for the NN method. Default: "DNA_NN_SantaLucia_2004"
#' 
#' @param tmm_table Thermodynamic parameters for terminal mismatches. Only applicable for the NN method.
#'   Default: "DNA_TMM_Bommarito_2000"
#' 
#' @param imm_table Thermodynamic parameters for internal mismatches. Only applicable for the NN method.
#'   Default: "DNA_IMM_Peyret_1999"
#' 
#' @param de_table Thermodynamic parameters for dangling ends. Only applicable for the NN method.
#'   Default: "DNA_DE_Bommarito_2000"
#' 
#' @param dnac_high Concentration of the higher concentrated strand in nM. Only applicable for the NN method.
#'   Default: 25
#' 
#' @param dnac_low Concentration of the lower concentrated strand in nM. Only applicable for the NN method.
#'   Default: 25
#' 
#' @param self_comp Logical value indicating if the sequence is self-complementary. Only applicable
#'   for the NN method. Default: FALSE
#' 
#' @param variant Empirical constants coefficient for GC method. Only applicable for the GC method.
#'   Default: "Primer3Plus"
#' 
#' @param userset A vector of four coefficient values for GC method. Only applicable for the GC method.
#'   Usersets override value sets. Default: NULL
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
#' @param salt_corr_method Method for calculating salt concentration corrections to the melting temperature.
#'   Available options:
#'   - "Schildkraut2010": Updated salt correction method
#'   - "Wetmur1991": Classic salt correction method
#'   - "SantaLucia1996": DNA-specific salt correction
#'   - "SantaLucia1998-1": Improved DNA salt correction
#'   - "SantaLucia1998-2": Alternative DNA salt correction
#'   - "Owczarzy2004": Comprehensive salt correction
#'   - "Owczarzy2008": Updated comprehensive salt correction
#'   Default: "Schildkraut2010"
#' 
#' @param DMSO Percent DMSO concentration in the reaction mixture. Default: 0
#' 
#' @param formamide_value_unit List containing formamide concentration value and unit. Default: list(value = 0, unit = "percent")
#'   - value: Numeric value of formamide concentration
#'   - unit: Either "percent" or "molar"
#' 
#' @param dmso_factor Coefficient of Tm decreases per percent DMSO. Default: 0.75
#'   Other published values are 0.5, 0.6 and 0.675.
#' 
#' @param formamide_factor Coefficient of Tm decrease per percent formamide. Default: 0.65
#'   Several papers report factors between 0.6 and 0.72.
#' 
#' @param mismatch Logical. If TRUE, every '.' in the sequence is counted as a mismatch.
#'   Only applicable for the GC method. Default: TRUE
#' 
#' @details
#' The function calculates melting temperature using the specified method(s). For each method:
#' - NN: Uses nearest neighbor thermodynamics with detailed sequence analysis
#' - GC: Uses GC content-based calculation with various empirical formulas
#' - Wallace: Uses the simple Wallace rule (2°C per A/T, 4°C per G/C)
#' 
#' The function processes the input sequence once and applies it to all selected methods,
#' making it more efficient than calling each method separately.
#' 
#' @return A list containing Tm values and options for each method used. The structure includes:
#'   - Tm: A list of sequences with updated Tm attributes
#'   - Options: A list containing calculation parameters and method information
#' 
#' @examples
#' # Calculate Tm using all methods
#' input_seq <- c("ATGCGATGCG", "ATGCGATGCGCCCGGAGATAG")
#' result <- tm_calculate(input_seq)
#' 
#' # Calculate Tm with specific method parameters
#' result <- tm_calculate(
#'   input_seq,
#'   method = "tm_nn",
#'   nn_table = "DNA_NN_SantaLucia_2004",
#'   salt_corr_method = "Owczarzy2008"
#' )
#' 
#' @export
tm_calculate <- function(input_seq,
                        complement_seq = NULL,
                        method = c("tm_nn", "tm_gc", "tm_wallace"),
                        ambiguous = FALSE,
                        shift = 0,
                        nn_table = c("DNA_NN_SantaLucia_2004",
                                    "DNA_NN_Breslauer_1986",
                                    "DNA_NN_Sugimoto_1996",
                                    "DNA_NN_Allawi_1998",
                                    "RNA_NN_Freier_1986",
                                    "RNA_NN_Xia_1998",
                                    "RNA_NN_Chen_2012",
                                    "RNA_DNA_NN_Sugimoto_1995"),
                        tmm_table = "DNA_TMM_Bommarito_2000",
                        imm_table = "DNA_IMM_Peyret_1999",
                        de_table = c("DNA_DE_Bommarito_2000",
                                    "RNA_DE_Turner_2010"),
                        dnac_high = 25,
                        dnac_low = 25,
                        self_comp = FALSE,
                        variant = c("Primer3Plus",
                                    "Chester1993",
                                    "QuikChange",
                                    "Schildkraut1965",
                                    "Wetmur1991_MELTING",
                                    "Wetmur1991_RNA",
                                    "Wetmur1991_RNA/DNA",
                                    "vonAhsen2001"),
                        userset = NULL,
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
                        DMSO = 0,
                        formamide_value_unit = list(value = 0, unit = "percent"),
                        dmso_factor = 0.75,
                        formamide_factor = 0.65,
                        mismatch = TRUE) {
  # Validate method argument
  method <- match.arg(method, several.ok = FALSE)
  
  # Process sequence once for all methods
  raw_seq <- process_seq(input_seq)
  # Initialize result list
  result <- list()
  
  # Calculate Tm using each selected method
  if ("tm_nn" %in% method) {
    result$tm_nn <- tm_nn(
      raw_seq = raw_seq,
      complement_seq = complement_seq,
      ambiguous = ambiguous,
      shift = shift,
      nn_table = nn_table,
      tmm_table = tmm_table,
      imm_table = imm_table,
      de_table = de_table,
      dnac_high = dnac_high,
      dnac_low = dnac_low,
      self_comp = self_comp,
      Na = Na,
      K = K,
      Tris = Tris,
      Mg = Mg,
      dNTPs = dNTPs,
      salt_corr_method = salt_corr_method,
      DMSO = DMSO,
      formamide_value_unit = formamide_value_unit,
      dmso_factor = dmso_factor,
      formamide_factor = formamide_factor
    )
  }
  
  if ("tm_gc" %in% method) {
    result$tm_gc <- tm_gc(
      raw_seq = raw_seq,
      ambiguous = ambiguous,
      userset = userset,
      variant = variant,
      Na = Na,
      K = K,
      Tris = Tris,
      Mg = Mg,
      dNTPs = dNTPs,
      salt_corr_method = salt_corr_method,
      mismatch = mismatch,
      DMSO = DMSO,
      formamide_value_unit = formamide_value_unit,
      dmso_factor = dmso_factor,
      formamide_factor = formamide_factor
    )
  }
  
  if ("tm_wallace" %in% method) {
    result$tm_wallace <- tm_wallace(
      raw_seq = raw_seq,
      ambiguous = ambiguous
    )
  }
  
  return(result)
} 
