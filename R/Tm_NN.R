#' Calculate melting temperature using nearest-neighbor thermodynamics
#'
#' @description
#' Calculate Tm of a nucleic acid duplex using the nearest-neighbor (NN)
#' model. The function supports DNA/DNA, RNA/RNA, RNA/DNA and 2'-O-MeRNA/RNA
#' duplexes plus a wide set of structural and chemical perturbations:
#' GU wobble pairs, single / tandem mismatches, single / double / long
#' dangling ends, internal loops, single / long bulge loops, CNG triplet
#' repeats, inosine (I), hydroxyadenine (H or A*), azobenzene (X_T / X_C),
#' single locked-nucleic-acid (LNA), and consecutive LNAs with or without a
#' single mismatch.
#'
#' Coverage mirrors the rmelting (Bioconductor) vignette sections 4.2.1 – 4.4.
#'
#' @details
#'
#' \strong{Argument-to-rmelting-section map}
#'
#' \tabular{lll}{
#'   rmelting argument                       \tab tm_nn argument        \tab rmelting section \cr
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
#' \strong{Modified-base / structural-feature encoding in the input sequence}
#'
#' \itemize{
#'   \item \code{I} — inosine
#'   \item \code{H} or \code{*} after the modified base — hydroxyadenine (A*)
#'   \item \code{X_T} (trans) or \code{X_C} (cis) — azobenzene
#'   \item \code{L} immediately before the base — locked nucleic acid (e.g. \code{AL}, \code{GL})
#'   \item \code{-} — bulge / internal-loop gap on the strand opposite the bulged base
#'   \item \code{.} (or leading/trailing \code{-}) — dangling-end placeholder
#' }
#'
#' \strong{Resolution priority for interior dinucleotide steps}
#'
#' \code{clna_mm > clna > lna > azb > ha > ino > gu > tandem > imm > nn}
#'
#' For each step the first table in this order that contains the key is used;
#' if no table contains the key the step is ignored and a warning is emitted.
#' The standard \code{nn_table} is always loaded and serves as the final
#' fallback for canonical Watson–Crick steps.
#'
#' \strong{Note on parameter tables}
#'
#' The NN tables (with rownames = dinucleotide keys, columns = \code{c("left","right")}
#' for ΔH / ΔS) live in \code{thermodynamic_nn_params} (see
#' \code{R/thermodynamic_nn_params.R}). This function calls
#' \code{get_table(name)} to fetch them. \strong{Any model identifier
#' exposed here but not yet present in \code{thermodynamic_nn_params}
#' will trigger a package-startup notice the first time the corresponding
#' feature is detected in a sequence}, and that contribution is skipped
#' (the rest of the calculation proceeds normally). This lets the API
#' expose the full rmelting catalogue while parameter values are filled
#' in incrementally from the primary literature.
#'
#' @param gr_seq Pre-processed sequence(s) in 5' to 3' direction. Output from
#'   \code{to_genomic_ranges()}.
#' @param ambiguous Logical. If \code{TRUE}, ambiguous IUPAC bases are kept.
#'   Default: \code{FALSE}.
#' @param shift Integer alignment offset (for dangling ends). Default: 0.
#' @param nn_table NN parameter table for perfectly matched steps (rmelting 4.2.1).
#'   See package data \code{thermodynamic_nn_params}.
#' @param tmm_table Terminal mismatch table.
#' @param imm_table Internal (single) mismatch table (rmelting 4.2.3).
#' @param tandem_table Tandem mismatch table (rmelting 4.2.4).
#' @param de_table Single dangling end table (rmelting 4.2.5).
#' @param dde_table Double dangling end table (rmelting 4.2.6).
#' @param lde_table Long dangling end table (rmelting 4.2.7).
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
#' @param dnac_high Concentration of the higher-concentrated strand, nM. Default 25.
#' @param dnac_low  Concentration of the lower-concentrated strand,  nM. Default 25.
#' @param self_comp Logical. \code{TRUE} if sequence is self-complementary.
#' @param Na,K,Tris,Mg,dNTPs Ion concentrations in mM. Defaults: 50,0,0,0,0.
#' @param salt_method Salt correction method. See \code{\link{salt_correction}}.
#' @param DMSO Percent DMSO. Default 0.
#' @param formamide_unit \code{list(value, unit)}. Default \code{list(0,"percent")}.
#' @param dmso_factor Tm drop per percent DMSO. Default 0.75.
#' @param formamide_factor Tm drop per percent formamide. Default 0.65.
#'
#' @return A \code{TmCalculator} list with elements \code{gr} (GRanges with
#'   \code{Tm} and \code{GC} metadata) and \code{options}.
#'
#' @references
#' Allawi H & SantaLucia J (1997) \emph{Biochemistry} 36:10581.
#' Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>.
#' SantaLucia J et al. (1996) \emph{Biochemistry} 35:3555.
#' SantaLucia J & Hicks D (2004) <doi:10.1146/annurev.biophys.32.110601.141800>.
#' Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>.
#' Tanaka F et al. (2004) \emph{Biochemistry} 43:7143.
#' Freier S (1986) <doi:10.1073/pnas.83.24.9373>.
#' Xia T (1998) <doi:10.1021/bi9809425>.
#' Sugimoto N (1995) <doi:10.1016/S0048-9697(98)00088-6>.
#' Kierzek E et al. (2006) \emph{Biochemistry} 45:581 (2'-O-MeRNA/RNA).
#' Mathews D H & Turner D H (1999) <doi:10.1006/jmbi.1999.2700>  (GU wobble).
#' Chen JL / Serra MJ (2012) <doi:10.1021/bi3002709>             (GU wobble + RNA NN).
#' Peyret N (1999) <doi:10.1021/bi9825091>                       (single MM, DNA).
#' Watkins N E et al. (2011) \emph{Nucleic Acids Res} 39:1894    (single MM, DNA/RNA).
#' Lu Z J et al. (2006) \emph{Nucleic Acids Res} 34:4912         (RNA single MM, IL, bulge).
#' Davis A R & Znosko B M (2007, 2008) \emph{Biochemistry}       (RNA single MM).
#' Bommarito S (2000) <doi:10.1093/nar/28.9.1929>                (DNA dangling ends).
#' Turner D H (2010) <doi:10.1093/nar/gkp892>                    (RNA dangling ends).
#' Ohmichi T et al. (2002) \emph{J Am Chem Soc} 124:10367        (DNA / RNA dangling ends).
#' O'Toole A S, Miller S, Serra M J (2005, 2006, 2008)            (RNA dangling ends).
#' Badhwar J et al. (2007) \emph{Biochemistry} 46:14715          (RNA internal loop).
#' Tan Z J & Chen S J (2007) \emph{Biophys J} 92:3615            (DNA bulge).
#' Blose J M et al. (2007) \emph{Biochemistry} 46:15123          (RNA bulge).
#' Broda M et al. (2005) <doi:10.1021/bi0501447>                 (CNG repeats).
#' SantaLucia J (2005) \emph{Nucleic Acids Res} 33:6258          (inosine, DNA).
#' Wright D J et al. (2007) \emph{Biochemistry} 46:4625          (inosine, RNA).
#' Kawakami J / Sugimoto N et al. (2001) <doi:10.1021/bi010918b> (hydroxyadenine).
#' Asanuma H et al. (2005) \emph{Angew Chem Int Ed} 44:2747      (azobenzene).
#' McTigue P M et al. (2004) \emph{Biochemistry} 43:5388         (single LNA).
#' Owczarzy R et al. (2011) \emph{Biochemistry} 50:9352          (single + consecutive LNA).
#'
#' @author Junhui Li
#' @encoding UTF-8
#' @export tm_nn
#'
#' @examples
#' input_seq <- c("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC",
#'                "AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
#' seqs <- to_genomic_ranges(input_seq)
#' out <- tm_nn(seqs, Na = 50)
#' out

tm_nn <- function(gr_seq,
                  ambiguous       = FALSE,
                  shift           = 0,
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
                  dnac_high       = 25,
                  dnac_low        = 25,
                  self_comp       = FALSE,
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
                  DMSO            = 0,
                  formamide_unit  = list(value = 0, unit = "percent"),
                  dmso_factor     = 0.75,
                  formamide_factor = 0.65) {

  # -- Validate args ----------------------------------------------------------
  nn_table       <- match.arg(nn_table)
  tmm_table      <- match.arg(tmm_table)
  imm_table      <- match.arg(imm_table)
  tandem_table   <- match.arg(tandem_table)
  de_table       <- match.arg(de_table)
  dde_table      <- match.arg(dde_table)
  lde_table      <- match.arg(lde_table)
  il_table       <- match.arg(il_table)
  sbl_table      <- match.arg(sbl_table)
  lbl_table      <- match.arg(lbl_table)
  gu_table       <- match.arg(gu_table)
  cng_table      <- match.arg(cng_table)
  ino_table      <- match.arg(ino_table)
  ha_table       <- match.arg(ha_table)
  azb_table      <- match.arg(azb_table)
  lna_table      <- match.arg(lna_table)
  clna_table     <- match.arg(clna_table)
  clna_mm_table  <- match.arg(clna_mm_table)
  salt_method    <- match.arg(salt_method)

  # -- Soft warning when GU table is layered onto Chen 2012 (already has GU)--
  if (nn_table == "RNA_NN_Chen_2012" && gu_table != "none") {
    warning(sprintf(
      "RNA_NN_Chen_2012 already includes GU wobble parameters; layering gu_table = '%s' may double-count.",
      gu_table))
  }

  # -- Load tables (NULL when 'none' or not yet populated) -------------------
  nn_tbl       <- .get_tm_table(nn_table,  required = TRUE)
  tmm_tbl      <- .get_tm_table(tmm_table, required = FALSE)
  imm_tbl      <- .get_tm_table(imm_table, required = FALSE)
  de_tbl       <- .get_tm_table(de_table,  required = FALSE)
  tandem_tbl   <- .get_tm_table(tandem_table)
  dde_tbl      <- .get_tm_table(dde_table)
  lde_tbl      <- .get_tm_table(lde_table)
  il_tbl       <- .get_tm_table(il_table)
  sbl_tbl      <- .get_tm_table(sbl_table)
  lbl_tbl      <- .get_tm_table(lbl_table)
  gu_tbl       <- .get_tm_table(gu_table)
  cng_tbl      <- .get_tm_table(cng_table)
  ino_tbl      <- .get_tm_table(ino_table)
  ha_tbl       <- .get_tm_table(ha_table)
  azb_tbl      <- .get_tm_table(azb_table)
  lna_tbl      <- .get_tm_table(lna_table)
  clna_tbl     <- .get_tm_table(clna_table)
  clna_mm_tbl  <- .get_tm_table(clna_mm_table)

  # -- N-filter + region-id setup --------------------------------------------
  region_ids <- names(gr_seq)
  if (is.null(region_ids)) region_ids <- rep("", length(gr_seq))
  empty_id <- is.na(region_ids) | region_ids == ""
  if (any(empty_id)) {
    region_ids[empty_id] <- paste0(
      as.character(GenomicRanges::seqnames(gr_seq))[empty_id], ":",
      GenomicRanges::start(gr_seq)[empty_id], "-",
      GenomicRanges::end(gr_seq)[empty_id]
    )
  }
  filtered_seq <- check_filter_seq(
    list(sequence   = gr_seq$sequence,
         complement = gr_seq$complement,
         region_ids = region_ids),
    method = "tm_nn"
  )
  gr_seq_dropoff      <- gr_seq[!filtered_seq$kept]
  gr_seq              <- gr_seq[filtered_seq$kept]
  gr_seq$sequence     <- filtered_seq$sequence
  gr_seq$complement   <- filtered_seq$complement

  # -- Per-sequence calculation ----------------------------------------------
  n   <- length(gr_seq)
  tm  <- numeric(n)
  gc  <- numeric(n)

  all_seqs  <- as.character(mcols(gr_seq)$sequence)
  all_cseqs <- as.character(mcols(gr_seq)$complement)

  for (i in seq_len(n)) {
    seq_str  <- all_seqs[i]
    cseq_str <- all_cseqs[i]
    if (nchar(seq_str) < 2L) { tm[i] <- NA_real_; gc[i] <- NA_real_; next }
    res <- tryCatch(
      .tm_nn_core(
        seq_str, cseq_str, ambiguous, shift,
        nn_tbl, tmm_tbl, imm_tbl, de_tbl,
        tandem_tbl, dde_tbl, lde_tbl,
        il_tbl, sbl_tbl, lbl_tbl,
        gu_tbl, cng_tbl, ino_tbl, ha_tbl, azb_tbl,
        lna_tbl, clna_tbl, clna_mm_tbl,
        cng_name = cng_table,
        dnac_high, dnac_low, self_comp,
        Na, K, Tris, Mg, dNTPs,
        salt_fn = salt_method, DMSO,
        dmso_factor, formamide_factor, formamide_unit
      ),
      error = function(e) list(Tm = NA_real_, GC = NA_real_)
    )
    tm[i] <- res$Tm
    gc[i] <- res$GC
  }
  if (!"GC" %in% names(GenomicRanges::mcols(gr_seq))) gr_seq$GC <- gc
  gr_seq$Tm <- tm

  # -- Reference list ---------------------------------------------------------
  ref_list <- .tm_nn_ref_list()
  fmt_ref <- function(tbl) {
    if (is.null(tbl) || identical(tbl, "none")) return("none")
    paste0(tbl, ": ", ifelse(is.null(ref_list[[tbl]]), "(reference TBD)", ref_list[[tbl]]))
  }

  gr_seq <- .normalize_tm_gc_metadata(gr_seq)

  result_list <- list(
    gr = gr_seq,
    options = list(
      "Ambiguous"                                  = ambiguous,
      "Shift"                                      = shift,
      "Thermodynamic NN values"                    = fmt_ref(nn_table),
      "Thermodynamic values for terminal mismatches" = fmt_ref(tmm_table),
      "Thermodynamic values for single (internal) mismatches" = fmt_ref(imm_table),
      "Thermodynamic values for tandem mismatches" = fmt_ref(tandem_table),
      "Thermodynamic values for single dangling ends" = fmt_ref(de_table),
      "Thermodynamic values for double dangling ends" = fmt_ref(dde_table),
      "Thermodynamic values for long dangling ends"   = fmt_ref(lde_table),
      "Thermodynamic values for internal loops"     = fmt_ref(il_table),
      "Thermodynamic values for single bulge loops" = fmt_ref(sbl_table),
      "Thermodynamic values for long bulge loops"   = fmt_ref(lbl_table),
      "GU wobble table"                             = fmt_ref(gu_table),
      "CNG repeat table"                            = fmt_ref(cng_table),
      "Inosine table"                               = fmt_ref(ino_table),
      "Hydroxyadenine table"                        = fmt_ref(ha_table),
      "Azobenzene table"                            = fmt_ref(azb_table),
      "Single LNA table"                            = fmt_ref(lna_table),
      "Consecutive LNA table"                       = fmt_ref(clna_table),
      "Consecutive LNA + single MM table"           = fmt_ref(clna_mm_table),
      "Concentration of higher strand (nM)"         = dnac_high,
      "Concentration of lower strand (nM)"          = dnac_low,
      "Sequence self-complementary"                 = self_comp,
      "Na (mM)"  = Na, "K (mM)" = K, "Tris (mM)" = Tris,
      "Mg (mM)"  = Mg, "dNTPs (mM)" = dNTPs,
      "Salt correction method"   = salt_method,
      "Percent DMSO"             = DMSO,
      "Formamide concentration"  = formamide_unit$value,
      "Formamide unit"           = formamide_unit$unit,
      "DMSO factor"              = dmso_factor,
      "Formamide factor"         = formamide_factor,
      "Skipped regions containing N" = gr_seq_dropoff
    )
  )
  class(result_list) <- c("TmCalculator", "list")
  attr(result_list, "nonhidden") <- "gr"
  result_list
}


# =============================================================================
#  Internal helpers
# =============================================================================

# -- Fetch a table from the package data, returning NULL when absent ----------
.get_tm_table <- function(name, required = FALSE) {
  if (is.null(name) || identical(name, "none")) return(NULL)
  tbl <- tryCatch(get_table(name), error = function(e) NULL)
  if (is.null(tbl) && required) {
    stop(sprintf("Required thermodynamic table '%s' not found in package data.", name))
  }
  if (is.null(tbl)) {
    packageStartupMessage(sprintf(
      "tm_nn: table '%s' is selected but not yet populated in thermodynamic_nn_params; contributions skipped.",
      name))
  }
  tbl
}

# -- Single source of truth for reference strings -----------------------------
.tm_nn_ref_list <- function() {
  list(
    # Perfect matching (4.2.1)
    "DNA_NN_Breslauer_1986"      = "Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>",
    "DNA_NN_Sugimoto_1996"       = "Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>",
    "DNA_NN_Allawi_1998"         = "Allawi H (1998) <doi:10.1093/nar/26.11.2694>",
    "DNA_NN_SantaLucia_2004"     = "SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>",
    "DNA_NN_AllSan_1997"         = "Allawi H & SantaLucia J (1997) Biochemistry 36:10581",
    "DNA_NN_SantaLucia_1996"     = "SantaLucia J et al. (1996) Biochemistry 35:3555",
    "DNA_NN_Tanaka_2004"         = "Tanaka F et al. (2004) Biochemistry 43:7143",
    "RNA_NN_Freier_1986"         = "Freier S (1986) <doi:10.1073/pnas.83.24.9373>",
    "RNA_NN_Xia_1998"            = "Xia T (1998) <doi:10.1021/bi9809425>",
    "RNA_NN_Chen_2012"           = "Chen JL / Serra MJ (2012) <doi:10.1021/bi3002709>",
    "RNA_DNA_NN_Sugimoto_1995"   = "Sugimoto N (1995) <doi:10.1016/S0048-9697(98)00088-6>",
    "MeRNA_RNA_NN_Kierzek_2006"  = "Kierzek E et al. (2006) Biochemistry 45:581",
    # Terminal MM
    "DNA_TMM_Bommarito_2000"     = "Bommarito S (2000) <doi:10.1093/nar/28.9.1929>",
    # Single MM (4.2.3)
    "DNA_IMM_Peyret_1999"        = "Peyret N (1999) <doi:10.1021/bi9825091>; Allawi & SantaLucia (1997/1998); SantaLucia (2005) <doi:10.1093/nar/gki918>",
    "DNA_RNA_IMM_Watkins_2011"   = "Watkins N E et al. (2011) Nucleic Acids Res 39:1894",
    "RNA_IMM_Lu_2006"            = "Lu Z J et al. (2006) Nucleic Acids Res 34:4912",
    "RNA_IMM_Davis_Znosko_2007"  = "Davis A R & Znosko B M (2007) Biochemistry 46:13425",
    "RNA_IMM_Davis_Znosko_2008"  = "Davis A R & Znosko B M (2008) Biochemistry 47:10178",
    # Tandem MM (4.2.4)
    "DNA_TandemMM_AllSanPey"     = "Allawi & SantaLucia (1997/1998); Peyret N (1999) <doi:10.1021/bi9825091>",
    "RNA_TandemMM_Mathews_1999"  = "Mathews D H et al. (1999) <doi:10.1006/jmbi.1999.2700>; Lu Z J et al. (2006)",
    # Single dangling end (4.2.5)
    "DNA_DE_Bommarito_2000"      = "Bommarito S (2000) <doi:10.1093/nar/28.9.1929>",
    "RNA_DE_Turner_2010"         = "Turner D H (2010) <doi:10.1093/nar/gkp892>",
    "DNA_DE_Ohmichi_2002"        = "Ohmichi T et al. (2002) J Am Chem Soc 124:10367",
    "RNA_DE_Ohmichi_2002"        = "Ohmichi T et al. (2002) J Am Chem Soc 124:10367",
    "RNA_DE_Serra_2008"          = "O'Toole A S, Miller S, Serra M J (2006/2008)",
    # Double dangling end (4.2.6)
    "DNA_DDE_Ohmichi_2002"       = "Ohmichi T et al. (2002) J Am Chem Soc 124:10367",
    "RNA_DDE_Ohmichi_2002"       = "Ohmichi T et al. (2002) J Am Chem Soc 124:10367",
    "RNA_DDE_OToole_2005"        = "O'Toole A S et al. (2005) Biochemistry 44:14914",
    "RNA_DDE_OToole_2006"        = "O'Toole A S et al. (2006) Nucleic Acids Res 34:3338",
    # Long dangling end (4.2.7)
    "DNA_LDE_Ohmichi_2002"       = "Ohmichi T et al. (2002) J Am Chem Soc 124:10367",
    "RNA_LDE_Ohmichi_2002"       = "Ohmichi T et al. (2002) J Am Chem Soc 124:10367",
    # Internal loop (4.2.8)
    "DNA_IL_SantaLucia_2004"     = "SantaLucia J & Hicks D (2004) <doi:10.1146/annurev.biophys.32.110601.141800>",
    "RNA_IL_Lu_2006"             = "Lu Z J et al. (2006) Nucleic Acids Res 34:4912",
    "RNA_IL_Badhwar_2007"        = "Badhwar J et al. (2007) Biochemistry 46:14715",
    # Single bulge loop (4.2.9)
    "DNA_SBL_Tanaka_2007"        = "Tan Z J & Chen S J (2007) Biophys J 92:3615",
    "DNA_SBL_SantaLucia_2004"    = "SantaLucia J & Hicks D (2004) <doi:10.1146/annurev.biophys.32.110601.141800>",
    "RNA_SBL_Lu_2006"            = "Lu Z J et al. (2006) Nucleic Acids Res 34:4912",
    "RNA_SBL_Blose_2007"         = "Blose J M et al. (2007) Biochemistry 46:15123",
    # Long bulge loop (4.2.10)
    "DNA_LBL_SantaLucia_2004"    = "SantaLucia J & Hicks D (2004) <doi:10.1146/annurev.biophys.32.110601.141800>",
    "RNA_LBL_Lu_2006"            = "Mathews D H et al. (1999); Lu Z J et al. (2006)",
    # GU wobble (4.2.2)
    "RNA_GU_Mathews_1999"        = "Mathews D H & Turner D H (1999) <doi:10.1006/jmbi.1999.2700>",
    "RNA_GU_Chen_2012"           = "Chen JL / Serra MJ (2012) <doi:10.1021/bi3002709>",
    # CNG (4.2.11)
    "DNA_CNG_Broda_2005"         = "Broda M et al. (2005) <doi:10.1021/bi0501447>",
    # Inosine (4.2.12)
    "DNA_INO_SantaLucia_2005"    = "Watkins & SantaLucia (2005) Nucleic Acids Res 33:6258",
    "RNA_INO_Wright_2007"        = "Wright D J et al. (2007) Biochemistry 46:4625",
    # Hydroxyadenine (4.2.13)
    "DNA_HA_Kawakami_2001"       = "Kawakami J / Sugimoto N et al. (2001) <doi:10.1021/bi010918b>",
    # Azobenzene (4.2.14)
    "DNA_AZB_Asanuma_2005"       = "Asanuma H et al. (2005) Angew Chem Int Ed 44:2747",
    # Single LNA (4.2.15)
    "DNA_LNA_Owczarzy_2011"      = "Owczarzy R et al. (2011) Biochemistry 50:9352",
    "DNA_LNA_McTigue_2004"       = "McTigue P M et al. (2004) Biochemistry 43:5388",
    # Consecutive LNA (4.3)
    "DNA_cLNA_Owczarzy_2011"     = "Owczarzy R et al. (2011) Biochemistry 50:9352",
    # Consecutive LNA + MM (4.4)
    "DNA_cLNA_MM_Owczarzy_2011"  = "Owczarzy R et al. (2011) Biochemistry 50:9352",
    "none"                       = "none"
  )
}

# -- Detect pure (CNG)n repeat (Broda 2005) -----------------------------------
.detect_cng <- function(seq_str) {
  if (grepl("^(C[ACGT]G){2,7}$", seq_str, perl = TRUE)) {
    list(is_cng = TRUE,
         N      = substr(seq_str, 2L, 2L),
         n      = nchar(seq_str) / 3L)
  } else {
    list(is_cng = FALSE)
  }
}

# -- Count leading / trailing dangling-end placeholders ('.' or '-') ----------
.count_dangle_ends <- function(s) {
  left  <- attr(regexpr("^[.-]*", s), "match.length")
  right <- attr(regexpr("[.-]*$", s), "match.length")
  c(left = max(0L, left), right = max(0L, right))
}

# -- Build a right-terminal dinucleotide key (reverse-orientation pair) -------
.right_key <- function(seq, cseq, len) {
  paste0(
    substr(cseq, len,     len),
    substr(cseq, len - 1L, len - 1L),
    "/",
    substr(seq,  len,     len),
    substr(seq,  len - 1L, len - 1L)
  )
}

# -- Add ΔH/ΔS for a key from the first table containing it -------------------
.lookup_step <- function(key, tables) {
  for (tbl in tables) {
    if (!is.null(tbl) && key %in% rownames(tbl)) {
      return(c(dH = tbl[key, 1L], dS = tbl[key, 2L], hit = 1))
    }
  }
  c(dH = 0, dS = 0, hit = 0)
}

# =============================================================================
#  Core single-sequence NN computation
# =============================================================================
.tm_nn_core <- function(seq_str, cseq_str, ambiguous, shift,
                         nn_tbl, tmm_tbl, imm_tbl, de_tbl,
                         tandem_tbl, dde_tbl, lde_tbl,
                         il_tbl, sbl_tbl, lbl_tbl,
                         gu_tbl, cng_tbl, ino_tbl, ha_tbl, azb_tbl,
                         lna_tbl, clna_tbl, clna_mm_tbl,
                         cng_name,
                         dnac_high, dnac_low, self_comp,
                         Na, K, Tris, Mg, dNTPs,
                         salt_fn, DMSO, dmso_factor, formamide_factor,
                         formamide_unit) {

  delta_h <- 0
  delta_s <- 0
  use_cng_path <- FALSE

  # -- 4.2.11 CNG fast-path (Broda 2005) ------------------------------------
  if (!is.null(cng_tbl)) {
    cng_info <- .detect_cng(seq_str)
    if (cng_info$is_cng) {
      cng_key <- paste0("C", cng_info$N, "G_", cng_info$n)
      if (cng_key %in% rownames(cng_tbl)) {
        delta_h <- cng_tbl[cng_key, 1L]
        delta_s <- cng_tbl[cng_key, 2L]
        use_cng_path <- TRUE
      } else {
        warning(sprintf("CNG key '%s' not found in table '%s'; falling back to NN.",
                        cng_key, cng_name))
      }
    }
  }

  if (!use_cng_path) {

    tmp_seq  <- seq_str
    tmp_cseq <- cseq_str

    # -- Pad sequences for shift / dangling-end alignment -------------------
    if (shift != 0 || nchar(tmp_seq) != nchar(tmp_cseq)) {
      if (shift > 0) {
        tmp_seq  <- paste0(strrep(".", shift), tmp_seq)
      } else if (shift < 0) {
        tmp_cseq <- paste0(strrep(".", -shift), tmp_cseq)
      }
      if (nchar(tmp_cseq) > nchar(tmp_seq))
        tmp_seq  <- paste0(tmp_seq,  strrep(".", nchar(tmp_cseq) - nchar(tmp_seq)))
      if (nchar(tmp_cseq) < nchar(tmp_seq))
        tmp_cseq <- paste0(tmp_cseq, strrep(".", nchar(tmp_seq)  - nchar(tmp_cseq)))

      while (substring(tmp_seq,  1L, 2L) == ".." ||
             substring(tmp_cseq, 1L, 2L) == "..") {
        tmp_seq  <- substring(tmp_seq,  2L, nchar(tmp_seq))
        tmp_cseq <- substring(tmp_cseq, 2L, nchar(tmp_cseq))
      }
      while ({ns <- nchar(tmp_seq);  substring(tmp_seq,  ns - 1L, ns)  == ".."} ||
             {nc <- nchar(tmp_cseq); substring(tmp_cseq, nc - 1L, nc) == ".."}) {
        tmp_seq  <- substring(tmp_seq,  1L, nchar(tmp_seq)  - 1L)
        tmp_cseq <- substring(tmp_cseq, 1L, nchar(tmp_cseq) - 1L)
      }
    }

    # -- 4.2.5 / 4.2.6 / 4.2.7  Dangling ends — single vs double vs long ----
    # Choose dangling-end table by the count of '.'/'-' at each terminus.
    ends_left  <- max(.count_dangle_ends(tmp_seq)["left"],
                      .count_dangle_ends(tmp_cseq)["left"])
    ends_right <- max(.count_dangle_ends(tmp_seq)["right"],
                      .count_dangle_ends(tmp_cseq)["right"])

    pick_de <- function(k) {
      if (k <= 0L) return(NULL)
      if (k == 1L) return(list(de_tbl))
      if (k == 2L) return(list(dde_tbl, de_tbl))    # fall through to single-DE
      return(list(lde_tbl, dde_tbl, de_tbl))        # ≥3 → long-DE
    }

    # -- Build interior dinucleotide keys ----------------------------------
    n     <- nchar(tmp_seq)
    n_int <- n - 1L
    if (n_int < 1L) {
      return(list(Tm = NA_real_, GC = .GC_fast(seq_str, ambiguous = ambiguous)))
    }
    fwd     <- substring(tmp_seq,  1L:n_int, 2L:(n_int + 1L))
    cfwd    <- substring(tmp_cseq, 1L:n_int, 2L:(n_int + 1L))
    keys_fr <- paste0(fwd, "/", cfwd)
    keys_t_left  <- keys_fr[1L]
    keys_t_right <- .right_key(tmp_seq, tmp_cseq, n)

    # -- Apply dangling-end contributions (terminal step only) -------------
    de_left_tables  <- pick_de(ends_left)
    de_right_tables <- pick_de(ends_right)

    if (!is.null(de_left_tables)) {
      hit <- .lookup_step(keys_t_left, de_left_tables)
      if (hit["hit"] == 1) {
        delta_h <- delta_h + hit["dH"]; delta_s <- delta_s + hit["dS"]
        keys_fr <- keys_fr[-1L]
        tmp_seq  <- substring(tmp_seq,  2L, nchar(tmp_seq))
        tmp_cseq <- substring(tmp_cseq, 2L, nchar(tmp_cseq))
        if (length(keys_fr) > 0L) keys_t_left <- keys_fr[1L]
      }
    }
    if (!is.null(de_right_tables) && length(keys_fr) > 0L) {
      hit <- .lookup_step(keys_t_right, de_right_tables)
      if (hit["hit"] == 1) {
        delta_h <- delta_h + hit["dH"]; delta_s <- delta_s + hit["dS"]
        keys_fr <- keys_fr[-length(keys_fr)]
        n <- nchar(tmp_seq) - 1L
        tmp_seq  <- substring(tmp_seq,  1L, n)
        tmp_cseq <- substring(tmp_cseq, 1L, n)
        keys_t_right <- .right_key(tmp_seq, tmp_cseq, n)
      }
    }

    # -- Terminal mismatches ----------------------------------------------
    if (!is.null(tmm_tbl) && length(keys_fr) > 0L &&
        keys_t_left %in% rownames(tmm_tbl)) {
      delta_h <- delta_h + tmm_tbl[keys_t_left, 1L]
      delta_s <- delta_s + tmm_tbl[keys_t_left, 2L]
      keys_fr <- keys_fr[-1L]
      tmp_seq  <- substring(tmp_seq,  2L, nchar(tmp_seq))
      tmp_cseq <- substring(tmp_cseq, 2L, nchar(tmp_cseq))
    }
    if (!is.null(tmm_tbl) && length(keys_fr) > 0L &&
        keys_t_right %in% rownames(tmm_tbl)) {
      delta_h <- delta_h + tmm_tbl[keys_t_right, 1L]
      delta_s <- delta_s + tmm_tbl[keys_t_right, 2L]
      keys_fr <- keys_fr[-length(keys_fr)]
      n <- nchar(tmp_seq) - 1L
      tmp_seq  <- substring(tmp_seq,  1L, n)
      tmp_cseq <- substring(tmp_cseq, 1L, n)
    }

    # -- Initiation parameters (canonical NN) -----------------------------
    delta_h <- delta_h + nn_tbl["init", 1L]
    delta_s <- delta_s + nn_tbl["init", 2L]
    if (nchar(tmp_seq) > 0L && substr(tmp_seq, 1L, 1L) == "T") {
      delta_h <- delta_h + nn_tbl["init_5T/A", 1L]
      delta_s <- delta_s + nn_tbl["init_5T/A", 2L]
    }
    if (nchar(tmp_seq) > 0L) {
      first_base <- substr(tmp_seq, 1L, 1L)
      last_base  <- substr(tmp_seq, nchar(tmp_seq), nchar(tmp_seq))
      gc_ends    <- sum(c(first_base, last_base) %in% c("G","C"))
      at_ends    <- 2L - gc_ends
      if (gc_ends == 0L) {
        delta_h <- delta_h + nn_tbl["init_allA/T", 1L]
        delta_s <- delta_s + nn_tbl["init_allA/T", 2L]
      } else {
        delta_h <- delta_h + nn_tbl["init_oneG/C", 1L]
        delta_s <- delta_s + nn_tbl["init_oneG/C", 2L]
      }
      delta_h <- delta_h + nn_tbl["init_A/T", 1L] * at_ends
      delta_s <- delta_s + nn_tbl["init_A/T", 2L] * at_ends
      delta_h <- delta_h + nn_tbl["init_G/C", 1L] * gc_ends
      delta_s <- delta_s + nn_tbl["init_G/C", 2L] * gc_ends
    }

    # -- Interior-step dispatch -------------------------------------------
    # Priority list — first table containing the key wins:
    # clna_mm > clna > lna > azb > ha > ino > gu > tandem > imm > nn
    priority_tables <- list(clna_mm_tbl, clna_tbl, lna_tbl,
                            azb_tbl, ha_tbl, ino_tbl,
                            gu_tbl, tandem_tbl, imm_tbl, nn_tbl)

    for (k in keys_fr) {
      hit <- .lookup_step(k, priority_tables)
      if (hit["hit"] == 1) {
        delta_h <- delta_h + hit["dH"]
        delta_s <- delta_s + hit["dS"]
      } else {
        # Unmatched step: could be inside an internal loop or a bulge —
        # consult IL / SBL / LBL tables. The exact loop / bulge sizing
        # logic is encoded as a key prefix in the table rownames, e.g.
        #   "IL_<size>"   for an internal loop of that size,
        #   "SBL_<base>"  for a single bulge of the given base,
        #   "LBL_<size>"  for a long bulge of that size.
        # Users add such rows when populating the structural tables.
        loop_hit <- .lookup_step(k, list(il_tbl, sbl_tbl, lbl_tbl))
        if (loop_hit["hit"] == 1) {
          delta_h <- delta_h + loop_hit["dH"]
          delta_s <- delta_s + loop_hit["dS"]
        }
        # else silently skip; missing-table notices are emitted at load time
      }
    }
  }  # end !use_cng_path

  # -- Self-complementary symmetry correction --------------------------------
  if (self_comp && "sym" %in% rownames(nn_tbl)) {
    delta_h <- delta_h + nn_tbl["sym", 1L]
    delta_s <- delta_s + nn_tbl["sym", 2L]
    k <- dnac_high * 1e-9
  } else {
    k <- (dnac_high - dnac_low / 2.0) * 1e-9
  }

  # -- Tm from thermodynamics + salt + chemical corrections ------------------
  R <- 1.987
  tm <- (1000 * delta_h) / (delta_s + R * log(k)) - 273.15

  if (!is.null(salt_fn) && !is.na(salt_fn)) {
    corr_salt <- salt_correction(
      Na = Na, K = K, Tris = Tris, Mg = Mg, dNTPs = dNTPs,
      method = salt_fn, input_seq = seq_str, ambiguous = ambiguous
    )
    if (salt_fn == "SantaLucia1998-2") {
      delta_s <- delta_s + corr_salt
      tm <- (1000 * delta_h) / (delta_s + R * log(k)) - 273.15
    } else if (salt_fn %in% c("Schildkraut2010", "Wetmur1991",
                              "SantaLucia1996",  "SantaLucia1998-1")) {
      tm <- tm + corr_salt
    } else if (salt_fn %in% c("Owczarzy2004", "Owczarzy2008")) {
      tm <- 1 / (1 / (tm + 273.15) + corr_salt) - 273.15
    }
  }

  pt_gc <- .GC_fast(seq_str, ambiguous = ambiguous)
  tm    <- tm + chem_correct(
    DMSO = DMSO,
    formamide_unit = formamide_unit,
    dmso_factor = dmso_factor,
    formamide_factor = formamide_factor,
    pt_gc = pt_gc
  )

  list(Tm = tm, GC = pt_gc)
}


# =============================================================================
#  Auxiliary helpers (kept for downstream callers)
# =============================================================================

# Package-private cache for thermodynamic tables (not .GlobalEnv)
.TM_TABLE_CACHE <- new.env(parent = emptyenv())

# Get a table from .TM_CONSTANTS / thermodynamic_nn_params with caching
get_table <- function(table_name) {
  if (!exists(table_name, envir = .TM_TABLE_CACHE, inherits = FALSE)) {
    val <- NULL
    if (exists(".TM_CONSTANTS", inherits = TRUE) &&
        is.list(.TM_CONSTANTS) && !is.null(.TM_CONSTANTS[[table_name]])) {
      val <- .TM_CONSTANTS[[table_name]]
    } else if (exists("thermodynamic_nn_params", inherits = TRUE) &&
               is.list(thermodynamic_nn_params) &&
               !is.null(thermodynamic_nn_params[[table_name]])) {
      val <- thermodynamic_nn_params[[table_name]]
    } else {
      stop(sprintf("Thermodynamic table '%s' not found.", table_name))
    }
    assign(table_name, val, envir = .TM_TABLE_CACHE)
  }
  get(table_name, envir = .TM_TABLE_CACHE, inherits = FALSE)
}

#' @keywords internal
.GC_fast <- function(seq_upper, ambiguous = FALSE) {
  n <- nchar(seq_upper)
  if (n == 0L) return(NA_real_)
  if (!ambiguous) {
    nGC <- n - nchar(gsub("[GC]", "", seq_upper, perl = TRUE))
    return(100 * nGC / n)
  }
  nG <- n - nchar(gsub("G", "", seq_upper, fixed = TRUE))
  nC <- n - nchar(gsub("C", "", seq_upper, fixed = TRUE))
  nS <- n - nchar(gsub("S", "", seq_upper, fixed = TRUE))
  nY <- n - nchar(gsub("Y", "", seq_upper, fixed = TRUE))
  nR <- n - nchar(gsub("R", "", seq_upper, fixed = TRUE))
  nK <- n - nchar(gsub("K", "", seq_upper, fixed = TRUE))
  nM <- n - nchar(gsub("M", "", seq_upper, fixed = TRUE))
  nB <- n - nchar(gsub("B", "", seq_upper, fixed = TRUE))
  nV <- n - nchar(gsub("V", "", seq_upper, fixed = TRUE))
  gc_count <- nG + nC + nS + 0.5 * (nY + nR + nK + nM) + (2/3) * (nB + nV)
  100 * gc_count / n
}

.rev_str <- function(s) paste(rev(strsplit(s, "", fixed = TRUE)[[1]]), collapse = "")

#' Fast complement and reverse complement
#'
#' Uses \code{chartr()} instead of character-by-character translation.
#'
#' @param seq_str Character string or vector.
#' @param rev Logical. Return reverse complement? Default \code{FALSE}.
#' @return Character string(s) with complement.
#' @export
complement_fast <- function(seq_str, rev = FALSE) {
  comp <- chartr("ACGTacgtWSRYKMBVDH", "TGCAtgcaWSYRMKVBHD", seq_str)
  if (rev) {
    vapply(comp, function(s) paste(rev(strsplit(s, "", fixed = TRUE)[[1]]), collapse = ""),
           character(1L), USE.NAMES = FALSE)
  } else {
    comp
  }
}

#' @keywords internal
benchmark_tm_nn <- function(n_seqs = 1000L, seq_len = 200L) {
  set.seed(42)
  bases <- c("A", "C", "G", "T")
  seqs <- vapply(seq_len(n_seqs), function(i) {
    paste(sample(bases, seq_len, replace = TRUE), collapse = "")
  }, character(1))
  gr_all <- to_genomic_ranges(seqs)
  t_batch <- system.time({
    out_b <- tm_nn(gr_all, Na = 50, salt_method = "Owczarzy2004")
    tm_batch <- as.numeric(S4Vectors::mcols(out_b$gr)$Tm)
  })
  cat(sprintf("Batch tm_nn() on %d seq x %d bp: %.2f sec\n",
              n_seqs, seq_len, t_batch["elapsed"]))
  invisible(tm_batch)
}
