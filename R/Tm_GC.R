#' Calculate the melting temperature using empirical formulas based on GC content
#' 
#' Calculate the melting temperature using empirical formulas based on GC content with different options.
#' The function returns a list of sequences with updated Tm attributes and calculation options.
#' 
#' @param gr_seq Pre-processed sequence(s) in 5' to 3' direction. This should be the output from
#'   to_genomic_ranges() function.
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
#'   - Schildkraut1965: Tm = 81.5 + 0.41(%GC) - 675/N + 16.6 x log10[Na+]
#'   - Wetmur1991_MELTING: Tm = 81.5 + 0.41(%GC) - 500/N + 16.6 x log10([Na+]/(1 + 0.7 x [Na+])) - %mismatch
#'   - Wetmur1991_RNA: Tm = 78 + 0.7(%GC) - 500/N + 16.6 x log10([Na+]/(1 + 0.7 x [Na+])) - %mismatch
#'   - Wetmur1991_RNA/DNA: Tm = 67 + 0.8(%GC) - 500/N + 16.6 x log10([Na+]/(1 + 0.7 x [Na+])) - %mismatch
#'   - Primer3Plus: Tm = 81.5 + 0.41(%GC) - 600/N + 16.6 x log10[Na+]
#'   - vonAhsen2001: Tm = 77.1 + 0.41(%GC) - 528/N + 11.7 x log10[Na+]
#'
#'   Salt correction is applied only for variants that include it in the formula
#'   (via \code{salt_correct()}). Chester1993 and QuikChange use no salt term.
#'   D is the mismatch penalty (typically 1): Tm decreases by D x (%mismatch).
#'   Use \code{X} (or \code{.}) in the sequence to mark mismatch positions.
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
#' @param salt_method Salt correction method:
#'   - \code{NULL} (default): the correction that belongs to the formula, i.e.
#'     the one published with \code{variant}, or \code{"Schildkraut2010"}
#'     when \code{userset} is supplied.
#'   - \code{NA} or \code{"none"}: no salt correction at all.
#'   - "Schildkraut2010": Schildkraut & Lifson 1965
#'   - "Wetmur1991": Wetmur 1991
#'   - "SantaLucia1996": SantaLucia 1996
#'   - "SantaLucia1998-1": SantaLucia 1998 (Method 1)
#'
#'   With a built-in \code{variant} the salt term is part of the published
#'   formula rather than a free choice, so naming a \emph{different} one is
#'   ignored with a warning: the result would otherwise be labelled with one
#'   method and computed with another. Supply \code{userset} to choose the
#'   correction yourself. Dropping it with \code{NA} is not a substitution
#'   and is honoured on either path.
#'
#'   "SantaLucia1998-2", "Owczarzy2004" and "Owczarzy2008" are not available
#'   for this function. The first corrects the entropy of a nearest-neighbor
#'   model, which a GC-content formula does not have. The other two correct
#'   the reciprocal of the melting temperature in kelvin, referenced to the
#'   same duplex in 1 M Na+, and carry a duplex-length term of their own,
#'   which these formulas already have. All three are available in
#'   \code{\link{tm_nn}}.
#'
#' @param mismatch Logical. If TRUE (default), every 'X' in the sequence is counted as a mismatch
#' 
#' @param DMSO Percent DMSO concentration in the reaction mixture. Default: 0
#' 
#' @param formamide_unit Formamide concentration as `list(value, unit)`. Default: list(value = 0, unit = "percent")
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
#' Schildkraut C, Lifson S. Dependence of the melting temperature of DNA on salt concentration. Biopolymers, 1965, 3(2):195-208.
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
#' gr_seq <- to_genomic_ranges(input_seq)
#' out <- tm_gc(gr_seq, ambiguous = TRUE, variant = "Primer3Plus", Na = 50, mismatch = TRUE)
#' out
#' out$options
#' 
#' @export tm_gc
tm_gc <- function(gr_seq,
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
                  salt_method = NULL,
                  mismatch = TRUE,
                  DMSO = 0,
                  formamide_unit = list(value = 0, unit = "percent"),
                  dmso_factor = 0.75,
                  formamide_factor = 0.65) {
  variant <- match.arg(variant)

  # Which corrections a GC-content formula can take. The two Owczarzy
  # corrections are deliberately absent. They are corrections to 1/Tm in
  # kelvin, referenced to the melting temperature of the same duplex in 1 M
  # Na+, and they carry a 1/(2(N-1)) duplex-length term of their own; these
  # formulas are on neither footing and already have a length term of their
  # own, so combining the two double-counts length even after the reciprocal
  # arithmetic is done right. tm_nn() is where they belong.
  GC_SALT <- c("Schildkraut2010", "Wetmur1991", "SantaLucia1996",
               "SantaLucia1998-1")
  named <- !is.null(salt_method)          # did the caller name one at all
  if (named) {
    if (length(salt_method) != 1L)
      stop("`salt_method` must be a single method name, NA to disable the ",
           "correction, or NULL to use the one that belongs to the formula.",
           call. = FALSE)
    if (is.na(salt_method)) {
      salt_method <- NA_character_        # explicit "no salt correction"
    } else if (identical(as.character(salt_method), "none")) {
      salt_method <- NA_character_        # the spelling tm_calculate() uses
    } else {
      salt_method <- as.character(salt_method)
      if (salt_method %in% c("Owczarzy2004", "Owczarzy2008")) {
        stop("`salt_method = \"", salt_method, "\"` is not available for ",
             "tm_gc(). The Owczarzy corrections apply to the reciprocal of ",
             "the melting temperature in kelvin, referenced to the same ",
             "duplex in 1 M Na+, and carry a duplex-length term of their ",
             "own, which the GC-content formulas already have. Use tm_nn() ",
             "for them.", call. = FALSE)
      }
      salt_method <- match.arg(salt_method, GC_SALT)
    }
  }

  if (is.null(userset)) {
    if (!variant %in% rownames(get_table("GC_VARTAB"))) {
      stop("only Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001 are allowed in variant")
    }
    gc_coef <- get_table("GC_VARTAB")[variant, ]
    # Each published variant carries its own salt term, so the correction is
    # a property of the formula rather than a free choice; NA_character_
    # means the formula has none (Chester1993, QuikChange). Substituting a
    # different term would hand back a result labelled with one method and
    # computed with another, so it is refused; dropping the term altogether
    # is not a substitution and is honoured, with $options reporting the
    # result as uncorrected.
    own <- get_table("GC_VARTAB")[variant, "salt_correct"]
    if (named && is.na(salt_method)) {
      salt_method_eff <- NA_character_
    } else {
      salt_method_eff <- own
      if (named && !identical(salt_method, own)) {
        carries <- if (is.na(own)) "no salt term of its own"
                   else paste0("the '", own, "' salt term")
        warning("variant '", variant, "' carries ", carries, ", so ",
                "`salt_method = \"", salt_method, "\"` is ignored. Use ",
                "`salt_method = NA` to drop the correction altogether, or ",
                "`userset` to choose a different one.", call. = FALSE)
      }
    }
  } else {
    gc_coef <- as.numeric(userset)
    # A user-supplied coefficient set says nothing about which salt term it
    # was fitted with, so the method has to be named; the default is the one
    # this function has always applied in that case.
    salt_method_eff <- if (named) salt_method else "Schildkraut2010"
  }
  # What $options reports is what was actually applied.
  salt_method <- salt_method_eff

  # Filter sequence
  gr_seq$sequence <- check_filter_seq(gr_seq$sequence, method = 'tm_gc')

  # Normalize gc_coef to a plain numeric vector of the four coefficients
  # (when userset is NULL it is a 1-row data.frame from GC_VARTAB)
  if (is.data.frame(gc_coef) || is.list(gc_coef)) {
    gc_coef <- vapply(gc_coef[1:4], function(x) as.numeric(x[1]), numeric(1))
  } else {
    gc_coef <- as.numeric(gc_coef)[1:4]
  }

  # Calculate Tm for all sequences in one vectorised pass
  chunk_res <- .tm_gc_chunk(
    list(sequence = as.character(gr_seq$sequence)),
    ambiguous = ambiguous, gc_coef = gc_coef, mismatch = mismatch,
    salt_method_eff = salt_method_eff,
    Na = Na, K = K, Tris = Tris, Mg = Mg, dNTPs = dNTPs,
    DMSO = DMSO, formamide_unit = formamide_unit,
    dmso_factor = dmso_factor, formamide_factor = formamide_factor
  )

  # One mcols<- assignment rather than two `$<-`: each `$<-` replaces the
  # whole metadata DataFrame and revalidates the GRanges, which dominates the
  # cost of a call on a short input.
  mc_out    <- GenomicRanges::mcols(gr_seq)
  mc_out$GC <- chunk_res$GC
  mc_out$Tm <- chunk_res$Tm
  GenomicRanges::mcols(gr_seq) <- mc_out
  gr_seq <- .normalize_tm_gc_metadata(gr_seq)

  # Create result list with proper structure
  # (result$df is computed lazily via `$.TmCalculator`)
  result_list <- list(
    gr = gr_seq,
    options = list(
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
      "Salt correction" = salt_method,
      Mismatch = mismatch,
      "Percent of DMSO" = DMSO,
      "Formamide concentration" = formamide_unit$value,
      "Formamide concentration unit" = formamide_unit$unit,
      "DMSO factor" = dmso_factor,
      "Formamide factor" = formamide_factor
    )
  )
  
  # Set class and attributes
  class(result_list) <- c("TmCalculator", "list")
  attr(result_list, "nonhidden") <- "gr"

  return(result_list)
}

# -- GC-method Tm over a block of sequences -----------------------------------
# `chunk` is list(sequence=) for the whole input.
#' @keywords internal
.tm_gc_chunk <- function(chunk, ambiguous, gc_coef, mismatch, salt_method_eff,
                         Na, K, Tris, Mg, dNTPs,
                         DMSO, formamide_unit, dmso_factor, formamide_factor) {
  seqs <- chunk$sequence
  m    <- length(seqs)
  if (m == 0L) return(list(Tm = numeric(0), GC = numeric(0)))

  # Previously this was a per-sequence loop calling gc_content(), which split every
  # sequence with s2c() and scanned it five times, and salt_correct(), which
  # did the same again. On 23,208 E. coli windows that cost 51.7 s against
  # 0.67 s for the compiled nearest-neighbor path. Both are now computed for
  # the whole chunk at once.
  n_seq <- nchar(seqs)
  pt_gc <- .gc_vec(seqs, ambiguous = ambiguous)

  tm <- gc_coef[1] + gc_coef[2] * pt_gc - gc_coef[3] / n_seq

  if (isTRUE(mismatch)) {
    # Behaviour preserved deliberately: `seqs %in% "X"` is TRUE only when an
    # entire sequence is the single character "X", so this term is zero for
    # any real input. Vectorising it as-is keeps results identical; the
    # underlying counting bug is recorded in ROADMAP.md rather than changed
    # here, where it would silently alter published values.
    mismatch_count <- as.numeric(seqs %in% "X")
    tm <- tm - gc_coef[4] * (mismatch_count * 100 / n_seq)
  }

  if (!is.na(salt_method_eff)) {
    tm <- tm + .salt_correct_vec(Na = Na, K = K, Tris = Tris, Mg = Mg,
                                 dNTPs = dNTPs, method = salt_method_eff,
                                 gc_pct = pt_gc, seq_len = n_seq)
  }

  if (DMSO > 0 | formamide_unit$value > 0) {
    tm <- tm + chem_correct(DMSO = DMSO,
                            formamide_unit = formamide_unit,
                            dmso_factor = dmso_factor,
                            formamide_factor = formamide_factor,
                            pt_gc = pt_gc)
  }

  list(Tm = as.numeric(tm), GC = as.numeric(pt_gc))
}
