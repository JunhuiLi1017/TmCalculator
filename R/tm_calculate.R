#' Calculate melting temperature using multiple methods
#'
#' Calculates melting temperature using multiple methods:
#' \itemize{
#'   \item Nearest-Neighbor thermodynamics (\code{\link{tm_nn}})
#'   \item GC content-based method (\code{\link{tm_gc}})
#'   \item Wallace rule (\code{\link{tm_wallace}})
#' }
#'
#' For \code{method = "tm_nn"} this is a thin pass-through to \code{\link{tm_nn}};
#' all 18 model-selection arguments mirror that function's signature and cover
#' the full rmelting (Bioconductor) catalogue (Tutorial sections 4.2.1 – 4.4):
#' canonical NN, GU wobble, single / tandem / terminal mismatches, single /
#' double / long dangling ends, internal loops, single / long bulge loops, CNG
#' triplet repeats, inosine, hydroxyadenine, azobenzene and single /
#' consecutive locked nucleic acids (with or without a single mismatch).
#'
#' @section Available Options:
#'
#' \strong{Method Selection:}
#' \itemize{
#'   \item \code{method}: c("tm_nn", "tm_gc", "tm_wallace")
#' }
#'
#' \strong{Nearest-Neighbor (NN) model tables (all map 1:1 to rmelting):}
#' \tabular{lll}{
#'   rmelting argument                       \tab tm_calculate argument \tab rmelting section \cr
#'   method.nn                               \tab \code{nn_table}       \tab 4.2.1 \cr
#'   method.GU                               \tab \code{gu_table}       \tab 4.2.2 \cr
#'   method.singleMM                         \tab \code{imm_table}      \tab 4.2.3 \cr
#'   method.tandemMM                         \tab \code{tandem_table}   \tab 4.2.4 \cr
#'   method.single.dangle                    \tab \code{de_table}       \tab 4.2.5 \cr
#'   method.double.dangle                    \tab \code{dde_table}      \tab 4.2.6 \cr
#'   method.long.dangle                      \tab \code{lde_table}      \tab 4.2.7 \cr
#'   method.internal.loop                    \tab \code{il_table}       \tab 4.2.8 \cr
#'   method.single.bulge.loop                \tab \code{sbl_table}      \tab 4.2.9 \cr
#'   method.long.bulge.loop                  \tab \code{lbl_table}      \tab 4.2.10 \cr
#'   method.CNG                              \tab \code{cng_table}      \tab 4.2.11 \cr
#'   method.inosine                          \tab \code{ino_table}      \tab 4.2.12 \cr
#'   method.hydroxyadenine                   \tab \code{ha_table}       \tab 4.2.13 \cr
#'   method.azobenzenes                      \tab \code{azb_table}      \tab 4.2.14 \cr
#'   method.locked                           \tab \code{lna_table}      \tab 4.2.15 \cr
#'   method.consecutive.locked               \tab \code{clna_table}     \tab 4.3 \cr
#'   method.consecutive.locked.singleMM      \tab \code{clna_mm_table}  \tab 4.4 \cr
#' }
#'
#' Any *_table argument accepts the string \code{"none"} (the default for
#' optional features) to disable that contribution. Selecting a table that is
#' not yet present in \code{\link{thermodynamic_nn_params}} triggers a one-time
#' package-startup notice and the corresponding contribution is silently
#' skipped — the rest of the calculation proceeds normally.
#'
#' \strong{GC Method Options:}
#' \itemize{
#'   \item \code{variant}: "Primer3Plus" (default), "Chester1993", "QuikChange",
#'         "Schildkraut1965", "Wetmur1991_MELTING", "Wetmur1991_RNA",
#'         "Wetmur1991_RNA/DNA", "vonAhsen2001"
#' }
#'
#' \strong{Salt Correction Options:}
#' \itemize{
#'   \item \code{salt_method}: "Schildkraut2010" (default), "Wetmur1991",
#'         "SantaLucia1996", "SantaLucia1998-1", "SantaLucia1998-2",
#'         "Owczarzy2004", "Owczarzy2008"
#' }
#'
#' \strong{Formamide Unit Options:}
#' \itemize{
#'   \item \code{formamide_unit$unit}: "percent" (default), "molar"
#' }
#'
#' @param input_seq Input sequence(s) in 5' to 3' direction. Can be provided as either:
#'   - A character string (e.g., "ATGCG")
#'   - A path to a FASTA file containing the sequence(s)
#'   - A GRanges object with sequence and complement metadata should be provided if mismatch is TRUE
#'   - A character vector where each element is a string in the format
#'     "chr:start-end:strand:species" (e.g., "chr1:100-200:+:BSgenome.Hsapiens.UCSC.hg38").
#'     Strand is "+" for positive (default if not provided) or "-" for negative.
#'
#' @param complement_seq Complementary sequence(s) in 3' to 5' direction. If not provided,
#'   the function will automatically generate it from input_seq.
#'
#' @param method Method(s) to use for Tm calculation. One of "tm_nn" (default),
#'   "tm_gc", or "tm_wallace".
#'
#' @param ambiguous Logical. If \code{TRUE}, ambiguous IUPAC bases are kept and
#'   contribute fractionally to GC content. Default: \code{FALSE}.
#'
#' @param shift Integer alignment offset (NN only). Default: 0.
#'
#' @param nn_table NN parameter table for perfectly matched steps (rmelting 4.2.1).
#' @param tmm_table Terminal mismatch table.
#' @param imm_table Internal (single) mismatch table (rmelting 4.2.3).
#' @param tandem_table Tandem mismatch table (rmelting 4.2.4).
#' @param de_table Single dangling-end table (rmelting 4.2.5).
#' @param dde_table Double dangling-end table (rmelting 4.2.6).
#' @param lde_table Long dangling-end table (rmelting 4.2.7).
#' @param il_table Internal loop table (rmelting 4.2.8).
#' @param sbl_table Single bulge loop table (rmelting 4.2.9).
#' @param lbl_table Long bulge loop table (rmelting 4.2.10).
#' @param gu_table GU wobble table (rmelting 4.2.2).
#' @param cng_table CNG triplet repeat table (rmelting 4.2.11).
#' @param ino_table Inosine table (rmelting 4.2.12).
#' @param ha_table Hydroxyadenine table (rmelting 4.2.13).
#' @param azb_table Azobenzene table (rmelting 4.2.14).
#' @param lna_table Single locked-nucleic-acid table (rmelting 4.2.15).
#' @param clna_table Consecutive LNA table (rmelting 4.3).
#' @param clna_mm_table Consecutive LNA with single mismatch table (rmelting 4.4).
#'
#' @param dnac_high Concentration of the higher-concentrated strand (nM). Default 25.
#' @param dnac_low  Concentration of the lower-concentrated strand  (nM). Default 25.
#' @param self_comp Logical. \code{TRUE} if sequence is self-complementary. Default FALSE.
#'
#' @param variant Empirical-constants coefficient set for the GC method. Default "Primer3Plus".
#' @param userset A vector of four coefficient values for GC method. Overrides \code{variant} when not NULL.
#'
#' @param Na,K,Tris,Mg,dNTPs Ion concentrations in mM. Defaults: 50, 0, 0, 0, 0.
#' @param salt_method Salt correction method for Tm. Default "Schildkraut2010".
#' @param DMSO Percent DMSO in the reaction mixture. Default 0.
#' @param formamide_unit \code{list(value, unit)}. Default \code{list(value = 0, unit = "percent")}.
#' @param dmso_factor Tm drop per percent DMSO. Default 0.75.
#' @param formamide_factor Tm drop per percent formamide. Default 0.65.
#' @param mismatch Logical. If \code{TRUE}, every '.' in the sequence is counted as
#'   a mismatch (GC method only). Default \code{TRUE}.
#'
#' @details
#' The function calculates melting temperature using the specified method.
#' For NN it uses sequence-resolved nearest-neighbor thermodynamics with the
#' full rmelting model catalogue; GC uses empirical formulas based on GC
#' percentage; Wallace uses the simple 2°C/4°C rule.
#'
#' @return A \code{TmCalculator} list with:
#'   \item{\code{gr}}{The input \code{GRanges} with metadata columns \code{Tm}
#'     and \code{GC} (melting temperature in \eqn{^{\circ}}C and GC percent).}
#'   \item{\code{options}}{Calculation parameters and method information.}
#'   \item{\code{df}}{A data.frame representation of \code{gr}.}
#'
#' @encoding UTF-8
#' @author Junhui Li
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
#' }
#'
#' @export tm_calculate
tm_calculate <- function(input_seq,
                         method          = c("tm_nn", "tm_gc", "tm_wallace"),
                         complement_seq  = NULL,
                         ambiguous       = FALSE,
                         shift           = 0,

                         # ---- NN model tables (rmelting 4.2.1 – 4.4) -------
                         nn_table        = c("DNA_NN_SantaLucia_2004",
                                             "DNA_NN_Breslauer_1986",
                                             "DNA_NN_Sugimoto_1996",
                                             "DNA_NN_Allawi_1998",
                                             "DNA_NN_AllSan_1997",
                                             "DNA_NN_SantaLucia_1996",
                                             "DNA_NN_Tanaka_2004",
                                             "RNA_NN_Freier_1986",
                                             "RNA_NN_Xia_1998",
                                             "RNA_NN_Chen_2012",
                                             "RNA_DNA_NN_Sugimoto_1995",
                                             "MeRNA_RNA_NN_Kierzek_2006"),
                         tmm_table       = "DNA_TMM_Bommarito_2000",
                         imm_table       = c("DNA_IMM_Peyret_1999",
                                             "DNA_RNA_IMM_Watkins_2011",
                                             "RNA_IMM_Lu_2006",
                                             "RNA_IMM_Davis_Znosko_2007",
                                             "RNA_IMM_Davis_Znosko_2008"),
                         tandem_table    = c("none",
                                             "DNA_TandemMM_AllSanPey",
                                             "RNA_TandemMM_Mathews_1999"),
                         de_table        = c("DNA_DE_Bommarito_2000",
                                             "RNA_DE_Turner_2010",
                                             "DNA_DE_Ohmichi_2002",
                                             "RNA_DE_Ohmichi_2002",
                                             "RNA_DE_Serra_2008"),
                         dde_table       = c("none",
                                             "DNA_DDE_Ohmichi_2002",
                                             "RNA_DDE_Ohmichi_2002",
                                             "RNA_DDE_OToole_2005",
                                             "RNA_DDE_OToole_2006"),
                         lde_table       = c("none",
                                             "DNA_LDE_Ohmichi_2002",
                                             "RNA_LDE_Ohmichi_2002"),
                         il_table        = c("none",
                                             "DNA_IL_SantaLucia_2004",
                                             "RNA_IL_Lu_2006",
                                             "RNA_IL_Badhwar_2007"),
                         sbl_table       = c("none",
                                             "DNA_SBL_Tanaka_2007",
                                             "DNA_SBL_SantaLucia_2004",
                                             "RNA_SBL_Lu_2006",
                                             "RNA_SBL_Blose_2007"),
                         lbl_table       = c("none",
                                             "DNA_LBL_SantaLucia_2004",
                                             "RNA_LBL_Lu_2006"),
                         gu_table        = c("none",
                                             "RNA_GU_Mathews_1999",
                                             "RNA_GU_Chen_2012"),
                         cng_table       = c("none",
                                             "DNA_CNG_Broda_2005"),
                         ino_table       = c("none",
                                             "DNA_INO_SantaLucia_2005",
                                             "RNA_INO_Wright_2007"),
                         ha_table        = c("none",
                                             "DNA_HA_Kawakami_2001"),
                         azb_table       = c("none",
                                             "DNA_AZB_Asanuma_2005"),
                         lna_table       = c("none",
                                             "DNA_LNA_Owczarzy_2011",
                                             "DNA_LNA_McTigue_2004"),
                         clna_table      = c("none",
                                             "DNA_cLNA_Owczarzy_2011"),
                         clna_mm_table   = c("none",
                                             "DNA_cLNA_MM_Owczarzy_2011"),

                         # ---- Concentration / strand --------------------
                         dnac_high       = 25,
                         dnac_low        = 25,
                         self_comp       = FALSE,

                         # ---- GC-method options -------------------------
                         variant         = c("Primer3Plus",
                                             "Chester1993",
                                             "QuikChange",
                                             "Schildkraut1965",
                                             "Wetmur1991_MELTING",
                                             "Wetmur1991_RNA",
                                             "Wetmur1991_RNA/DNA",
                                             "vonAhsen2001"),
                         userset         = NULL,

                         # ---- Ions / salt correction --------------------
                         Na              = 50,
                         K               = 0,
                         Tris            = 0,
                         Mg              = 0,
                         dNTPs           = 0,
                         salt_method     = c("Schildkraut2010",
                                             "Wetmur1991",
                                             "SantaLucia1996",
                                             "SantaLucia1998-1",
                                             "SantaLucia1998-2",
                                             "Owczarzy2004",
                                             "Owczarzy2008"),

                         # ---- Chemical corrections ----------------------
                         DMSO            = 0,
                         formamide_unit  = list(value = 0, unit = "percent"),
                         dmso_factor     = 0.75,
                         formamide_factor = 0.65,
                         mismatch        = TRUE) {

  # Validate method argument
  method <- match.arg(method, several.ok = FALSE)

  # convert input_seq to genomic ranges
  if (inherits(input_seq, "GRanges")) {
    gr <- input_seq
  } else {
    gr <- to_genomic_ranges(input_seq = input_seq, complement_seq = complement_seq)
  }

  # Calculate Tm using the selected method
  if (method == "tm_nn") {
    result <- tm_nn(
      gr_seq          = gr,
      ambiguous       = ambiguous,
      shift           = shift,
      nn_table        = nn_table,
      tmm_table       = tmm_table,
      imm_table       = imm_table,
      tandem_table    = tandem_table,
      de_table        = de_table,
      dde_table       = dde_table,
      lde_table       = lde_table,
      il_table        = il_table,
      sbl_table       = sbl_table,
      lbl_table       = lbl_table,
      gu_table        = gu_table,
      cng_table       = cng_table,
      ino_table       = ino_table,
      ha_table        = ha_table,
      azb_table       = azb_table,
      lna_table       = lna_table,
      clna_table      = clna_table,
      clna_mm_table   = clna_mm_table,
      dnac_high       = dnac_high,
      dnac_low        = dnac_low,
      self_comp       = self_comp,
      Na              = Na,
      K               = K,
      Tris            = Tris,
      Mg              = Mg,
      dNTPs           = dNTPs,
      salt_method     = salt_method,
      DMSO            = DMSO,
      formamide_unit  = formamide_unit,
      dmso_factor     = dmso_factor,
      formamide_factor = formamide_factor
    )
  } else if (method == "tm_gc") {
    result <- tm_gc(
      gr_seq          = gr,
      ambiguous       = ambiguous,
      userset         = userset,
      variant         = variant,
      Na              = Na,
      K               = K,
      Tris            = Tris,
      Mg              = Mg,
      dNTPs           = dNTPs,
      salt_method     = salt_method,
      mismatch        = mismatch,
      DMSO            = DMSO,
      formamide_unit  = formamide_unit,
      dmso_factor     = dmso_factor,
      formamide_factor = formamide_factor
    )
  } else if (method == "tm_wallace") {
    result <- tm_wallace(
      gr_seq    = gr,
      ambiguous = ambiguous
    )
  }

  # Ensure a data.frame representation is always available
  if (!is.null(result$gr) && is.null(result$df)) {
    result$df <- as.data.frame(result$gr)
  }

  result
}
