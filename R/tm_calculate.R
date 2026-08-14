#' Calculate melting temperature using multiple methods
#'
#' Calculates nucleic acid melting temperature (Tm) by one of three methods, and
#' returns the result as a \code{GRanges} object so that Tm can be used directly
#' as a quantitative genomic feature alongside other assays:
#' \itemize{
#'   \item \strong{Nearest neighbor} (\code{tm_nn}) sums the stacking free
#'     energies of adjacent base-pair steps using experimentally derived
#'     enthalpy and entropy parameters, and so resolves sequences of identical
#'     base composition but different order. It supports salt and chemical
#'     corrections, mismatches and dangling ends, and is the default.
#'   \item \strong{GC content} (\code{tm_gc}) computes Tm as an empirical
#'     function of GC percentage with corrections for length and ionic strength.
#'     Cheaper than NN, but blind to sequence order.
#'   \item \strong{Wallace rule} (\code{tm_wallace}) assigns a fixed
#'     contribution per base. It is calibrated for short oligonucleotides,
#'     typically 14 to 20 bp, and ignores sequence context and reaction
#'     conditions, so it is not appropriate for long sequences or for
#'     genome-wide windows.
#' }
#'
#' The input sequence is processed once and passed to the selected method, which
#' is faster than calling the individual functions separately.
#'
#' @section Salt handling:
#' Most nearest-neighbor parameter sets were fitted at a single reference sodium
#' concentration, and other conditions are reached through the
#' \code{salt_method} correction formulas. The Weber/VarGibbs sets were instead
#' fitted directly at a stated sodium concentration and are meant to replace
#' salt correction. When such a set is selected and \code{Na} matches the
#' concentration it was fitted at, correction is skipped automatically; when it
#' does not, correction is applied with a warning. See \code{\link{tm_nn}} for
#' details.
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
#'       \item DNA/DNA: "DNA_NN_Breslauer_1986", "DNA_NN_Sugimoto_1996",
#'         "DNA_NN_Allawi_1998", "DNA_NN_SantaLucia_2004" (default)
#'       \item DNA/DNA, salt-optimized: "DNA_NN_Weber_2015" (1020 mM),
#'         "DNA_NN_Weber_OW04_69", "..._119", "..._220", "..._621",
#'         "..._1020" (fitted at 69 to 1020 mM sodium)
#'       \item RNA/RNA: "RNA_NN_Freier_1986", "RNA_NN_Xia_1998",
#'         "RNA_NN_Chen_2012"
#'       \item RNA/RNA, salt-optimized: "RNA_NN_Weber_VIF_71", "..._121",
#'         "..._221", "..._621", "..._1021" and the corresponding
#'         "RNA_NN_Weber_FIF_*" sets
#'       \item RNA/DNA: "RNA_DNA_NN_Sugimoto_1995",
#'         "RNA_DNA_NN_Weber_2019_FT", "RNA_DNA_NN_Weber_2019_VH" (1000 mM),
#'         "RNA_DNA_NN_Weber_2019_LS" (100 mM)
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
#'   \item \code{salt_method}:
#'     \itemize{
#'       \item "Schildkraut2010" (default)
#'       \item "Wetmur1991"
#'       \item "SantaLucia1996"
#'       \item "SantaLucia1998-1"
#'       \item "Owczarzy2004"
#'       \item "Owczarzy2008"
#'       \item "none" (also selected automatically when \code{nn_table} was
#'         fitted at the requested \code{Na})
#'     }
#' }
#' 
#' \strong{Formamide Unit Options:}
#' \itemize{
#'   \item \code{formamide_unit$unit}:
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
#'   - A GRanges object with sequence and complement metadata should be provided if mismatch is TRUE
#'   - A character vector where each element is a string in the format "chr:start-end:strand:species" (e.g., "chr1:100-200:+:BSgenome.Hsapiens.UCSC.hg38"). Strand is "+" for positive (default if not provided) or "-" for negative.
#'     - chr: Chromosome ID
#'     - start: Start position
#'     - end: End position
#'     - strand: positive or negtive strand
#'     - species:  Species name for reference genome (e.g., "BSgenome.Hsapiens.UCSC.hg38"), see \code{BSgenome::available.genomes()} for all available genomes. please make sure the genome package is installed, otherwise the function will stop.
#' 
#' @param complement_seq Complementary sequence(s) in 3' to 5' direction. If not provided,
#'   the function will automatically generate it from input_seq. This is the template/target
#'   sequence that the input sequence will hybridize with. Can be provided as input_seq format besides A NULL value(default)
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
#'   Only applicable for the NN method. Sets whose name encodes a sodium
#'   concentration were fitted at that condition and are not salt-corrected
#'   again. See \code{\link{tm_nn}} for the full list and guidance on choosing
#'   between them. Default: "DNA_NN_SantaLucia_2004"
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
#' @param salt_method Salt correction method for Tm. Default: "Schildkraut2010"
#'   Available options:
#'   - "none": Disables salt correction. Also selected automatically when the
#'     chosen \code{nn_table} was fitted at the requested \code{Na}.
#'   - "Schildkraut2010": Updated salt correction method
#'   - "Wetmur1991": Classic salt correction method
#'   - "SantaLucia1996": DNA-specific salt correction
#'   - "SantaLucia1998-1": Improved DNA salt correction
#'   - "Owczarzy2004": Comprehensive salt correction
#'   - "Owczarzy2008": Updated comprehensive salt correction
#'   Default: "Schildkraut2010"
#' 
#' @param DMSO Percent DMSO concentration in the reaction mixture. Default: 0
#' 
#' @param formamide_unit Formamide concentration as `list(value, unit)`. Default: list(value = 0, unit = "percent")
#'   - value: Numeric value of formamide concentration
#'   - unit: Either "percent" or "molar"
#' 
#' @param dmso_factor Coefficient of Tm decreases per percent DMSO. Default: 0.75
#'   Other published values are 0.5, 0.6 and 0.675.
#' 
#' @param formamide_factor Tm decrease per percent formamide. Default: 0.65
#'   Several papers report factors between 0.6 and 0.72.
#' 
#' @param mismatch Logical. If TRUE, every '.' in the sequence is counted as a mismatch.
#'   Only applicable for the GC method. Default: TRUE
#' 
#' @details
#' The three methods differ in resolution and in the range of sequence lengths
#' over which they are calibrated, so they are not interchangeable.
#'
#' \code{tm_nn} is the appropriate default. Because it sums sequence-dependent
#' stacking terms, it distinguishes sequences of identical GC content but
#' different base order, which the other two cannot. Its parameters were derived
#' from short duplexes under a two-state assumption; when applied to long
#' sequences or to fixed-width genomic windows the resulting value is best read
#' as a relative measure of local thermodynamic stability for comparison across
#' windows, rather than as an absolute experimental melting temperature.
#'
#' \code{tm_gc} computes Tm from GC percentage using one of several published
#' empirical formulas selected by \code{variant}, with corrections for length
#' and ionic strength. It extends to longer sequences at low computational cost
#' but cannot resolve base order.
#'
#' \code{tm_wallace} applies the 2 + 4 rule, adding 2 \eqn{^{\circ}}C per A or T
#' and 4 \eqn{^{\circ}}C per G or C. It is calibrated for short oligonucleotides,
#' typically 14 to 20 bp, and takes no account of sequence context, salt or
#' chemical additives. Accuracy degrades quickly with length, so it is retained
#' for compatibility rather than recommended for genome-scale work.
#'
#' Salt and chemical corrections apply to \code{tm_nn} and \code{tm_gc} only.
#' The input sequence is parsed and validated once and reused by the selected
#' method, which is faster than calling the individual functions directly.
#' 
#' @return A \code{TmCalculator} list with:
#'   \item{\code{gr}}{The input \code{GRanges} with metadata columns \code{Tm}
#'     and \code{GC} (melting temperature in \eqn{^{\circ}}C and GC percent).}
#'   \item{\code{options}}{Calculation parameters and method information. For
#'     the nearest-neighbor method this includes \code{Salt correction applied}
#'     (logical) and \code{Parameter set fitted at [Na+] (mM)}, which record
#'     whether a salt correction was actually performed.}
#' 
#' @encoding UTF-8
#' @author Junhui Li
#' 
#' @export
#' 
#' @importFrom BSgenome available.genomes
#' @importFrom GenomeInfoDb genome
#' 
#' @examples
#' \dontrun{
#' input_seq <- c("chr1:1000100-1000150:+:BSgenome.Hsapiens.UCSC.hg38")
#' result <- tm_calculate(
#'   input_seq,
#'   method = "tm_nn",
#'   nn_table = "DNA_NN_SantaLucia_2004",
#'   salt_method = "Owczarzy2008"
#' )
#'
#' # A hybrid parameter set fitted at 100 mM sodium. Salt correction is
#' # skipped automatically because Na matches the fitted condition.
#' result_ls <- tm_calculate(
#'   input_seq,
#'   method = "tm_nn",
#'   nn_table = "RNA_DNA_NN_Weber_2019_LS",
#'   Na = 100
#' )
#' }
#'
#' @seealso \code{\link{tm_nn}} for the nearest-neighbor method and the full
#'   list of thermodynamic parameter sets.
#'
#' @export tm_calculate
tm_calculate <- function(input_seq,
                        method = c("tm_nn", "tm_gc", "tm_wallace"),
                        complement_seq = NULL,
                        ambiguous = FALSE,
                        shift = 0,
                        nn_table = c("DNA_NN_SantaLucia_2004",
                                    "DNA_NN_Breslauer_1986",
                                    "DNA_NN_Sugimoto_1996",
                                    "DNA_NN_Allawi_1998",
                                    "RNA_NN_Freier_1986",
                                    "RNA_NN_Xia_1998",
                                    "RNA_NN_Chen_2012",
                                    "RNA_DNA_NN_Sugimoto_1995",
                                    "DNA_NN_Weber_2015",
                                    "DNA_NN_Weber_OW04_69",
                                    "DNA_NN_Weber_OW04_119",
                                    "DNA_NN_Weber_OW04_220",
                                    "DNA_NN_Weber_OW04_621",
                                    "DNA_NN_Weber_OW04_1020",
                                    "RNA_NN_Weber_VIF_71",
                                    "RNA_NN_Weber_VIF_121",
                                    "RNA_NN_Weber_VIF_221",
                                    "RNA_NN_Weber_VIF_621",
                                    "RNA_NN_Weber_VIF_1021",
                                    "RNA_NN_Weber_FIF_71",
                                    "RNA_NN_Weber_FIF_121",
                                    "RNA_NN_Weber_FIF_221",
                                    "RNA_NN_Weber_FIF_621",
                                    "RNA_NN_Weber_FIF_1021",
                                    "RNA_DNA_NN_Weber_2019_FT",
                                    "RNA_DNA_NN_Weber_2019_VH",
                                    "RNA_DNA_NN_Weber_2019_LS"),
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
                        salt_method = c("Schildkraut2010",
                                            "Wetmur1991",
                                            "SantaLucia1996",
                                            "SantaLucia1998-1",
                                            "Owczarzy2004",
                                            "Owczarzy2008",
                                            "none"),
                        DMSO = 0,
                        formamide_unit = list(value = 0, unit = "percent"),
                        dmso_factor = 0.75,
                        formamide_factor = 0.65,
                        mismatch = TRUE) {
  # Validate method argument
  method <- match.arg(method, several.ok = FALSE)
  
  # convert input_seq to genomic ranges
  if (inherits(input_seq, "GRanges")) {
    gr <- input_seq
  } else {  
    gr <- to_genomic_ranges(input_seq=input_seq, complement_seq = complement_seq)
  }

  # check and filter the sequence
  #gr$sequence <- check_filter_seq(gr$sequence, method)
  #gr$complement <- check_filter_seq(gr$complement, method)

  # Calculate Tm using each selected method
  if ("tm_nn" %in% method) {
    result <- tm_nn(
      gr_seq = gr,
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
      salt_method = salt_method,
      DMSO = DMSO,
      formamide_unit = formamide_unit,
      dmso_factor = dmso_factor,
      formamide_factor = formamide_factor
    )
  }
  
  if ("tm_gc" %in% method) {
    result <- tm_gc(
      gr_seq = gr,
      ambiguous = ambiguous,
      userset = userset,
      variant = variant,
      Na = Na,
      K = K,
      Tris = Tris,
      Mg = Mg,
      dNTPs = dNTPs,
      salt_method = salt_method,
      mismatch = mismatch,
      DMSO = DMSO,
      formamide_unit = formamide_unit,
      dmso_factor = dmso_factor,
      formamide_factor = formamide_factor
    )
  }
  
  if ("tm_wallace" %in% method) {
    result <- tm_wallace(
      gr_seq = gr,
      ambiguous = ambiguous
    )
  }

  # Ensure a data.frame representation is always available
  if (!is.null(result$gr) && is.null(result$df)) {
    result$df <- as.data.frame(result$gr)
  }

  result
}
