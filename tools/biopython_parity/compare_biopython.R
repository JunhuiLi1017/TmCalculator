#!/usr/bin/env Rscript
# ===========================================================================
# compare_biopython.R -- cross-check every tm_nn() parameter against Biopython
#
#   Rscript tools/biopython_parity/compare_biopython.R [--quick] [--out DIR]
#
# Needs Biopython on the same machine (pip install biopython) and python3 on
# the PATH. Writes cases.csv, biopython.csv and parity_report.csv to --out
# (default tools/biopython_parity/out) and prints a summary.
#
# ---------------------------------------------------------------------------
# What the two implementations should and should not agree on
#
# Suite A, perfect duplexes: must agree EXACTLY, for every parameter set, every
#   salt method and every salt/strand condition. Nothing in either model's
#   handling of a fully paired duplex is supposed to differ, so any row here
#   that is not "exact" is a real defect on one side or the other. This is the
#   assertion the script exists for; the exit status is non-zero if it breaks.
#
# Suite B, terminal mismatches and dangling ends: differences are EXPECTED and
#   deliberate. Biopython consumes the terminal mismatch or dangling end and
#   then indexes all four initiation terms on the original sequence anyway;
#   TmCalculator indexes them on the duplex that is left, because SantaLucia &
#   Hicks (2004) define the terminal penalty as a property of the closing base
#   pair and a mismatch is not a base pair. See block 15 of
#   tests/testthat/test_regressions_1_1_2.R. Rows here are reported but do not
#   fail the run.
#
# Suite C, internal mismatches only: initiation is untouched by an internal
#   mismatch, so these must agree exactly as well.
#
# Suite D, uncovered stacks: Biopython raises, TmCalculator returns NA. Both
#   refusing is agreement; one answering while the other refuses is not.
#
# A note on DNA_NN_Breslauer_1986. It is the only shipped set where
# init_allA/T and init_oneG/C differ, and the two implementations ask slightly
# different questions to choose between them: Biopython tests the whole input
# sequence for G or C, TmCalculator tests the trimmed duplex for a G.C pair.
# Those agree on every perfect duplex, which is why suite A still demands
# exactness for it; they can differ in suites B and D.
# ===========================================================================

suppressPackageStartupMessages({
  library(TmCalculator)
  library(GenomicRanges)
})

# -- is the package we just loaded actually the one in this working tree? ----
# library() loads the INSTALLED copy. Rebuilding R/sysdata.rda in the source
# tree does nothing until the package is reinstalled, and a stale install shows
# up here as a wall of suite-A differences that look like a parameter-table
# bug. Check it up front instead of letting it masquerade as one.
lib <- find.package("TmCalculator")
message("TmCalculator ", as.character(utils::packageVersion("TmCalculator")),
        " loaded from ", lib)
inst_sysdata <- file.path(lib, "R", "sysdata.rdb")
src_newer <- Filter(function(f) file.exists(f) && file.exists(inst_sysdata) &&
                      file.mtime(f) > file.mtime(inst_sysdata),
                    c("R/zzz.R", "R/sysdata.rda", "src/tm_nn_core.cpp",
                      "R/tm_nn.R"))
if (length(src_newer)) {
  warning("These source files are newer than the installed package:\n  ",
          paste(src_newer, collapse = "\n  "),
          "\nReinstall before trusting this report:  R CMD INSTALL .",
          call. = FALSE, immediate. = TRUE)
}

args   <- commandArgs(trailingOnly = TRUE)
quick  <- "--quick" %in% args
outdir <- if ("--out" %in% args) args[match("--out", args) + 1L] else
  file.path("tools", "biopython_parity", "out")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

TOL <- 1e-6          # degrees C; both sides are double precision on the same formula

# -- parameter axes ---------------------------------------------------------

NN_DNA <- c("DNA_NN_Breslauer_1986", "DNA_NN_Sugimoto_1996",
            "DNA_NN_Allawi_1998", "DNA_NN_SantaLucia_2004")
NN_RNA <- c("RNA_NN_Freier_1986", "RNA_NN_Xia_1998", "RNA_NN_Chen_2012")
NN_HYB <- "RNA_DNA_NN_Sugimoto_1995"
NN_ALL <- c(NN_DNA, NN_RNA, NN_HYB)

SALT <- c("none", "Schildkraut2010", "Wetmur1991", "SantaLucia1996",
          "SantaLucia1998-1", "SantaLucia1998-2", "Owczarzy2004",
          "Owczarzy2008")

# Each row is one ionic condition. Owczarzy2004/2008 are only defined with
# magnesium present or absent in specific regimes, so the grid deliberately
# spans monovalent-only, mixed, and magnesium-dominated.
SALTCOND <- data.frame(
  Na    = c( 50, 1000,  50,  50,   0,  50),
  K     = c(  0,    0,  50,   0,  50,   0),
  Tris  = c(  0,    0,  10,   0,  10,   0),
  Mg    = c(  0,    0,   0,   3,   3,  10),
  dNTPs = c(  0,    0,   0,   0, 0.8, 0.8)
)

DNAC <- data.frame(dnac1 = c(25, 250,  50),
                   dnac2 = c(25,   0,  10))

SEQS <- c("GCATCGTAGGCTAGCT",
          "ATGCGCGCAT",
          "ATATATATATAT",
          "GCGCGCGCGCGC",
          "TTTTTTTTTTTTAAAA",
          "ACGTACGTACGTACGTACGT")
PALINDROME <- "GCGCGCGC"          # for self_comp = TRUE

if (quick) {
  NN_ALL   <- c("DNA_NN_SantaLucia_2004", "DNA_NN_Breslauer_1986", NN_HYB)
  SALT     <- c("none", "SantaLucia1998-2", "Owczarzy2008")
  SALTCOND <- SALTCOND[c(1, 4, 5), , drop = FALSE]
  DNAC     <- DNAC[1:2, , drop = FALSE]
  SEQS     <- SEQS[1:3]
}

de_for <- function(nn) if (startsWith(nn, "RNA_NN")) "RNA_DE_Turner_2010" else
  "DNA_DE_Bommarito_2000"

rc <- function(s) chartr("ACGT", "TGCA", s)
flip <- function(s) vapply(strsplit(s, "", fixed = TRUE),
                           function(x) paste(rev(x), collapse = ""), character(1))

mutate_at <- function(s, i, to) {
  ch <- strsplit(s, "", fixed = TRUE)[[1L]]
  ch[i] <- to
  paste(ch, collapse = "")
}

# -- case construction ------------------------------------------------------

new_cases <- function(suite, seq, c_seq, shift = 0L, self_comp = FALSE) {
  grid <- expand.grid(nn_table = NN_ALL, salt_method = SALT,
                      cond = seq_len(nrow(SALTCOND)), dnac = seq_len(nrow(DNAC)),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  data.frame(
    suite       = suite,
    seq         = seq,
    c_seq       = c_seq,
    shift       = as.integer(shift),
    nn_table    = grid$nn_table,
    tmm_table   = "DNA_TMM_Bommarito_2000",
    imm_table   = "DNA_IMM_Peyret_1999",
    de_table    = unname(vapply(grid$nn_table, de_for, character(1))),
    dnac1       = DNAC$dnac1[grid$dnac],
    dnac2       = DNAC$dnac2[grid$dnac],
    selfcomp    = self_comp,
    SALTCOND[grid$cond, , drop = FALSE],
    salt_method = grid$salt_method,
    stringsAsFactors = FALSE, row.names = NULL
  )
}

cases <- do.call(rbind, c(
  # A: perfect duplexes, both strands supplied explicitly
  lapply(SEQS, function(s) new_cases("A_perfect", s, rc(s))),
  # A: and the same duplex read from the other strand, which must also agree
  lapply(SEQS, function(s) new_cases("A_perfect_swapped", flip(rc(s)), flip(s))),
  # A: self-complementary
  list(new_cases("A_selfcomp", PALINDROME, rc(PALINDROME), self_comp = TRUE)),
  # B: one terminal mismatch, at the 5' end, the 3' end, and both
  list(new_cases("B_tmm_5p",   SEQS[1], mutate_at(rc(SEQS[1]), 1L, "G"))),
  list(new_cases("B_tmm_3p",   SEQS[1],
                 mutate_at(rc(SEQS[1]), nchar(SEQS[1]), "G"))),
  list(new_cases("B_tmm_both", SEQS[1],
                 mutate_at(mutate_at(rc(SEQS[1]), 1L, "G"),
                           nchar(SEQS[1]), "G"))),
  # B: dangling ends. The complement has to be the one that still pairs under
  # the shift, or every position becomes a mismatch and the case stops being
  # about dangling ends at all. Under both implementations seq[j] pairs with
  # c_seq[j + shift], so a 5' overhang needs a complement short by one at the
  # front and a 3' overhang needs one extra base at the front.
  list(new_cases("B_dangle_5p", SEQS[1], substring(rc(SEQS[1]), 2L),
                 shift = -1L)),
  list(new_cases("B_dangle_3p", SEQS[1], paste0("A", rc(SEQS[1])),
                 shift =  1L)),
  # and the same overhangs reached through unequal lengths at shift 0, where
  # the shorter strand is padded at its 3' end instead
  list(new_cases("B_short_cmp", SEQS[1],
                 substring(rc(SEQS[1]), 1L, nchar(SEQS[1]) - 1L))),
  list(new_cases("B_long_cmp",  SEQS[1], paste0(rc(SEQS[1]), "A"))),
  # C: internal mismatches only, which leave initiation alone
  list(new_cases("C_imm_mid",  SEQS[1], mutate_at(rc(SEQS[1]), 8L, "G"))),
  list(new_cases("C_imm_two",  SEQS[1],
                 mutate_at(mutate_at(rc(SEQS[1]), 5L, "G"), 11L, "C"))),
  # D: a stack neither table covers, where both sides should refuse
  list(new_cases("D_uncovered", "AGGTCA", "TGAGTA"))
))
cases$id <- seq_len(nrow(cases))
cases <- cases[c("id", "suite", "seq", "c_seq", "shift", "nn_table",
                 "tmm_table", "imm_table", "de_table", "dnac1", "dnac2",
                 "selfcomp", "Na", "K", "Tris", "Mg", "dNTPs", "salt_method")]

f_cases <- file.path(outdir, "cases.csv")
f_bio   <- file.path(outdir, "biopython.csv")
write.csv(cases, f_cases, row.names = FALSE)
message(sprintf("built %d cases across %d suites", nrow(cases),
                length(unique(cases$suite))))

# -- Biopython half ---------------------------------------------------------

py <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])),
                "biopython_tm.py")
if (!file.exists(py)) py <- "tools/biopython_parity/biopython_tm.py"
status <- system2("python3", c(shQuote(py), shQuote(f_cases), shQuote(f_bio)))
if (status != 0L) stop("biopython_tm.py failed with status ", status)
bio <- read.csv(f_bio, colClasses = c("integer", "numeric", "character"))

# -- TmCalculator half ------------------------------------------------------

tmc_one <- function(r) {
  gr <- TmCalculator::to_genomic_ranges(r$seq, complement_seq = r$c_seq)
  out <- withCallingHandlers(
    tryCatch(
      TmCalculator::tm_nn(gr, shift = r$shift, nn_table = r$nn_table,
                          tmm_table = r$tmm_table, imm_table = r$imm_table,
                          de_table = r$de_table,
                          dnac_high = r$dnac1, dnac_low = r$dnac2,
                          self_comp = r$selfcomp,
                          Na = r$Na, K = r$K, Tris = r$Tris, Mg = r$Mg,
                          dNTPs = r$dNTPs, salt_method = r$salt_method),
      error = function(e) structure(conditionMessage(e), class = "tmc_error")),
    warning = function(w) invokeRestart("muffleWarning"))
  if (inherits(out, "tmc_error")) return(list(tm = NA_real_, err = as.character(out)))
  list(tm = as.numeric(GenomicRanges::mcols(out$gr)$Tm)[1L], err = "")
}

res <- lapply(seq_len(nrow(cases)), function(i) tmc_one(cases[i, ]))
cases$tmc_tm    <- vapply(res, `[[`, numeric(1),   "tm")
cases$tmc_error <- vapply(res, `[[`, character(1), "err")

# -- compare ----------------------------------------------------------------

rep <- merge(cases, bio, by = "id", all.x = TRUE)
rep <- rep[order(rep$id), ]
rep$diff <- rep$tmc_tm - rep$bio_tm

tmc_na <- is.na(rep$tmc_tm)
bio_na <- is.na(rep$bio_tm)
rep$verdict <- ifelse(
  tmc_na & bio_na, "both refused",
  ifelse(tmc_na | bio_na, "ONE REFUSED",
         ifelse(abs(rep$diff) <= TOL, "exact", "differs")))

# Suites B and D are allowed to differ; A and C are not.
rep$ok <- rep$verdict %in% c("exact", "both refused") |
  startsWith(rep$suite, "B_") | startsWith(rep$suite, "D_")

f_rep <- file.path(outdir, "parity_report.csv")
write.csv(rep, f_rep, row.names = FALSE)

# -- summary ----------------------------------------------------------------

cat("\n", strrep("=", 72), "\n", sep = "")
tab <- table(rep$suite, rep$verdict)
print(tab)

cat("\nlargest absolute difference per suite (degrees C)\n")
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_      # a suite where both sides refused
}                                          # every case has no difference to report
agg <- aggregate(abs(rep$diff), list(suite = rep$suite), safe_max)
names(agg)[2] <- "max_abs_diff"
print(agg, row.names = FALSE)

bad <- rep[!rep$ok, ]
if (nrow(bad)) {
  # Group rather than dump rows: a systematic fault repeats across every salt
  # and ionic condition, and 200 near-identical lines hide which axis it is on.
  cat("\n", nrow(bad), " case(s) that must have matched and did not.\n",
      "By suite, parameter set and verdict:\n\n", sep = "")
  grp <- aggregate(
    list(n = bad$id, worst = abs(bad$diff)),
    list(suite = bad$suite, nn_table = bad$nn_table, verdict = bad$verdict),
    function(x) x)
  grp$n     <- lengths(grp$n)
  grp$worst <- vapply(grp$worst, safe_max, numeric(1))
  print(grp[order(grp$suite, grp$nn_table), ], row.names = FALSE)

  cat("\nsalt methods involved: ",
      paste(sort(unique(bad$salt_method)), collapse = ", "), "\n", sep = "")
  cat("one example per group:\n\n")
  ex <- bad[!duplicated(bad[c("suite", "nn_table", "verdict")]), ]
  print(head(ex[c("suite", "nn_table", "salt_method", "seq", "c_seq",
                  "tmc_tm", "bio_tm", "diff")], 20), row.names = FALSE)

  if (length(src_newer)) {
    cat("\nNOTE: the working tree is newer than the installed package (see the",
        "\nwarning above). Reinstall and re-run before reading anything into",
        "\nthese differences.\n")
  }
  cat("\nfull report: ", f_rep, "\n", sep = "")
  quit(status = 1L)
}

cat("\nEvery perfect-duplex and internal-mismatch case matches Biopython to ",
    format(TOL), " degrees C.\n", sep = "")
cat("Suite B differences are the documented initiation convention; see the",
    "\nheader of this script and block 15 of test_regressions_1_1_2.R.\n")
cat("full report: ", f_rep, "\n", sep = "")
