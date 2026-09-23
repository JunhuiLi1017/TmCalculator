#!/usr/bin/env Rscript
# ===========================================================================
# bench_crosstool.R -- output consistency and computational cost of
#                      TmCalculator relative to existing Tm calculators
#
# WHAT THIS DOES AND DOES NOT MEASURE
#
# The comparison is deliberately restricted to a plain list of
# oligonucleotides, which is the one task all three tools were designed to
# perform. Benchmarking MELTING 5 or Biopython on millions of genomic
# windows would compare tools built for different purposes and prove
# nothing; the coordinate model that only TmCalculator provides is a
# capability difference and is reported as a feature table, not as a timing.
#
# The timings are therefore a comparison of THREE TOOLS AS DISTRIBUTED, not
# of three languages. TmCalculator's nearest-neighbor loop is compiled C++,
# Biopython's Tm_NN is pure Python, and MELTING 5 is Java. That difference
# is an implementation choice by each project and is exactly what a user
# experiences, but the results must not be described as "R is faster than
# Python".
#
# Three quantities are produced:
#   1. Output consistency -- mean and maximum |dTm| between tools on an
#      identical input under an identical model. This is the part that is
#      independent of language, hardware and intended use, and it is the
#      part the reviewer asked for by name.
#   2. Throughput -- sequences per second at increasing input size, with
#      interpreter/JVM start-up measured separately so that it does not
#      dominate the small input sizes.
#   3. Peak resident set size of the whole process. R, Java and Python
#      reserve memory differently, so this indicates the practical
#      footprint of each workflow rather than the memory the calculation
#      itself requires.
#
# BEFORE TRUSTING ANY TIMING, RUN THE CALIBRATION STEP
#
#   Rscript inst/scripts/bench_crosstool.R --calibrate
#
# It pushes a handful of sequences through all three tools and prints the
# per-sequence Tm side by side. The three projects use different argument
# names, different units for strand concentration, and different defaults
# for the salt correction. If the calibration output does not agree to
# within rounding, the benchmark is measuring default-value differences
# rather than implementation differences, and its conclusions are void.
# Adjust the TOOL PARAMETERS block below until it agrees, then run:
#
#   Rscript inst/scripts/bench_crosstool.R --outdir bench_crosstool_out
#
# Two diagnostics help when the calibration does not agree:
#
#   --grid     sweep MELTING 5's nucleic.acid.conc conventions against its
#              ion corrections and report which pair matches
#   --nosalt   remove the salt term from all three tools, which separates a
#              disagreement in the salt correction from one in the dH/dS
#              accumulation. Combine with --calibrate.
#
# REQUIREMENTS
#   R:      TmCalculator, rmelting (which pulls in melting5jars), a Java runtime
#   Python: python3 with biopython
#   Shell:  GNU time (`gtime` on macOS via `brew install gnu-time`) or the
#           BSD /usr/bin/time -l that ships with macOS
# ===========================================================================

## -- Command line ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
has_flag  <- function(flag) flag %in% args

outdir    <- argval("--outdir", "bench_crosstool_out")
n_rep     <- as.integer(argval("--reps", "3"))
seq_len   <- as.integer(argval("--seqlen", "200"))
timeout_s <- as.integer(argval("--timeout", "7200"))
sizes     <- as.integer(strsplit(argval("--sizes", "1000,10000,100000"), ",")[[1]])
python    <- argval("--python", "python3")
calibrate <- has_flag("--calibrate")
grid      <- has_flag("--grid")
nosalt    <- has_flag("--nosalt")
seed      <- as.integer(argval("--seed", "20260903"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# TOOL PARAMETERS -- all three tools must be given the SAME model.
#
# Nearest-neighbor set: SantaLucia & Hicks 2004 (the "unified" parameters).
#   TmCalculator  nn_table = "DNA_NN_SantaLucia_2004"
#   Biopython     nn_table = mt.DNA_NN4
#   MELTING 5     see TODO below
#
# MELTING 5 defaults to method.nn = "all97" (Allawi & SantaLucia 1997) and
# correction.ion = "ahs01" (von Ahsen 2001). Neither is what we want, and
# leaving them unset was the reason the first calibration run showed a
# systematic 2.4 C offset. Both are now set explicitly.
#
# Salt correction: 12.5 x log10[Na+] (SantaLucia et al. 1996). This is the
# only plain additive log term present in all three tools. The 11.7 x log10
# form (SantaLucia 1998 method 1) has no entry in MELTING 5's correction.ion
# list, so it cannot serve as common ground however natural it is elsewhere
# in this package.
#   TmCalculator  salt_method    = "SantaLucia1996"
#   Biopython     saltcorr       = 3
#   MELTING 5     correction.ion = "san96"
#
# Strand concentration is the remaining unknown. Biopython forms
# k = (dnac1 - dnac2/2) * 1e-9 and uses R*ln(k); MELTING 5 takes a single
# `nucleic.acid.conc` whose convention the R wrapper does not document.
# Rather than guess, `--grid` sweeps the plausible conventions against the
# candidate ion corrections and reports which combination agrees.
# ===========================================================================
NA_MM      <- 50      # mM sodium
DNAC_HIGH  <- 25      # nM, the strand in excess
DNAC_LOW   <- 25      # nM, the limiting strand
SELF_COMP  <- FALSE

SALT_TMCALC    <- "SantaLucia1996"
SALT_BIOPYTHON <- 3L
MELTING_ION    <- "san96"
MELTING_NN     <- "san04"
MELTING_NA_M   <- NA_MM / 1000        # mol/L

# Which model MELTING 5 uses is decided by `method.nn`, which the runner sets
# explicitly to "san04". A nearest-neighbor method given explicitly is applied
# whatever the sequence length, so the 200 bp windows used by the case study
# are evaluated with the same model as in the other two tools. The
# `size.threshold` default of 60 governs the automatic choice between the
# nearest-neighbor and approximative formulas, not an override of an explicit
# one.
#
# It is passed anyway, raised above the sequence length, so that the run does
# not depend on that default at all. --melting-threshold 60 restores MELTING's
# own value, which is a cheap way to confirm the two agree.
MELTING_SIZE_THRESHOLD <- local({
  v <- argval("--melting-threshold", "auto")
  if (identical(v, "auto")) max(60L, seq_len + 1L) else as.integer(v)
})
message("MELTING 5 size.threshold = ", MELTING_SIZE_THRESHOLD,
        " (sequence length ", seq_len, ")")
# Effective concentration Biopython uses, in mol/L. The grid varies this.
MELTING_CONC   <- (DNAC_HIGH - DNAC_LOW / 2) * 1e-9

## -- --nosalt: bisect the discrepancy ---------------------------------------
# Removing the salt term splits the comparison in two. If the tools then
# agree, the disagreement lives in the salt correction (the composition of
# [Mon], or where the correction is applied); if it survives, it lives in
# the dH/dS accumulation and the two-state formula.
#
# TmCalculator and Biopython can switch the correction off outright.
# MELTING 5 cannot: `correction.ion` has no "none" entry. It is instead
# evaluated at 1 M [Na+], where every log10-type correction -- including
# san96, 12.5 x log10(1) -- is exactly zero. That is also the reference
# condition at which the parameters were fitted, so it is the right place
# to compare thermodynamics rather than a trick.
if (nosalt) {
  message("--nosalt: salt correction disabled; MELTING 5 evaluated at 1 M [Na+]")
  SALT_TMCALC    <- "none"
  SALT_BIOPYTHON <- 0L
  MELTING_NA_M   <- 1.0
}

## -- Environment -----------------------------------------------------------
suppressPackageStartupMessages(library(TmCalculator))

have_rmelting <- requireNamespace("rmelting", quietly = TRUE)
if (!have_rmelting)
  message("NOTE: rmelting is not installed; it will be skipped.")

have_python <- nzchar(Sys.which(python)) &&
  system2(python, c("-c", shQuote("import Bio; print(Bio.__version__)")),
          stdout = NULL, stderr = NULL) == 0L
if (!have_python)
  message("NOTE: ", python, " with biopython not found; it will be skipped.")

# GNU time reports peak RSS in kB after "Maximum resident set size";
# the BSD time that ships with macOS reports it in bytes on a line ending
# "maximum resident set size". Detect which one is available.
time_cmd <- local({
  if (nzchar(Sys.which("gtime")))            list(bin = "gtime", flag = "-v", unit = 1024)
  else if (identical(Sys.info()[["sysname"]], "Darwin"))
                                             list(bin = "/usr/bin/time", flag = "-l", unit = 1)
  else if (nzchar(Sys.which("/usr/bin/time"))) list(bin = "/usr/bin/time", flag = "-v", unit = 1024)
  else NULL
})
if (is.null(time_cmd))
  stop("No usable `time` binary found; install GNU time (brew install gnu-time).")

## -- 1. Input sets ---------------------------------------------------------
# Sequences are sampled from the E. coli chromosome rather than generated
# uniformly at random, so that the GC distribution -- which drives both the
# Tm range and the Owczarzy-style corrections -- is realistic. Any window
# containing N is discarded rather than repaired.
make_sequences <- function(n, len, seed) {
  set.seed(seed)
  pkg <- "BSgenome.Ecoli.NCBI.ASM584v2"
  if (requireNamespace(pkg, quietly = TRUE)) {
    genome <- BSgenome::getBSgenome(pkg)
    chr    <- GenomeInfoDb::seqnames(genome)[1]
    L      <- GenomeInfoDb::seqlengths(genome)[[chr]]
    starts <- sample.int(L - len, n * 1.1)     # oversample, then drop N
    s <- as.character(Biostrings::Views(
      genome[[chr]], start = starts, width = len))
  } else {
    message("NOTE: ", pkg, " not installed; falling back to uniform random ",
            "sequences. The GC distribution will not be realistic.")
    s <- vapply(seq_len(n * 1.1), function(i)
      paste(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = ""),
      character(1))
  }
  s <- s[!grepl("[^ACGTacgt]", s)]
  toupper(head(s, n))
}

seq_file <- function(n) file.path(outdir, sprintf("seqs_%d.txt", n))
for (n in sizes) {
  if (!file.exists(seq_file(n)))
    writeLines(make_sequences(n, seq_len, seed), seq_file(n))
}
message("input sets written: ", paste(sizes, collapse = ", "), " sequences")

## -- 2. Runner scripts -----------------------------------------------------
# Each tool runs in its OWN process. This is the only way to attribute peak
# RSS to a tool rather than to the accumulated state of one long session,
# and it also means a tool that crashes or hangs does not take the others
# with it. Each runner reads a sequence file and writes one Tm per line.

make_runner_tmcalculator <- function(salt = SALT_TMCALC,
                                     hi = DNAC_HIGH, lo = DNAC_LOW) sprintf('
## Four timestamps, not two. Reading the input and writing the results are
## artefacts of running each tool as a separate process, and they scale with
## the input, so folding them into "start-up" makes that quantity grow with n
## and describe nothing.
tA   <- proc.time()[["elapsed"]]
suppressPackageStartupMessages(library(TmCalculator))
tB   <- proc.time()[["elapsed"]]
a    <- commandArgs(trailingOnly = TRUE)
seqs <- readLines(a[1])
tC   <- proc.time()[["elapsed"]]
res  <- tm_calculate(seqs, method = "tm_nn",
                     nn_table    = "DNA_NN_SantaLucia_2004",
                     salt_method = "%s",
                     Na = %s, dnac_high = %s, dnac_low = %s,
                     self_comp = %s)
tD   <- proc.time()[["elapsed"]]
writeLines(format(res$gr$Tm, digits = 10), a[2])
tE   <- proc.time()[["elapsed"]]
## timing goes to a file, not stdout: rmelting drives the Java engine through
## rJava, which redirects stdout, so a printed value can be swallowed
tim <- sprintf(c("LOAD_SECONDS %%.4f", "READ_SECONDS %%.4f",
                 "COMPUTE_SECONDS %%.4f", "WRITE_SECONDS %%.4f"),
               c(tB - tA, tC - tB, tD - tC, tE - tD))
if (length(a) >= 3L) writeLines(tim, a[3])
cat(tim, sep = "\\n")
', salt, NA_MM, hi, lo, ifelse(SELF_COMP, "TRUE", "FALSE"))

runner_tmcalculator <- make_runner_tmcalculator()

# method.nn and correction.ion MUST be given. Their defaults ("all97" and
# "ahs01") are a different model from the one the other two tools are
# running, which is what the first calibration run detected.
make_runner_rmelting <- function(conc, ion, nn, na_m = MELTING_NA_M) sprintf('
tA   <- proc.time()[["elapsed"]]
suppressPackageStartupMessages(library(rmelting))
tB   <- proc.time()[["elapsed"]]
a    <- commandArgs(trailingOnly = TRUE)
seqs <- readLines(a[1])
tC   <- proc.time()[["elapsed"]]
tm <- vapply(seqs, function(s) {
  r <- rmelting::melting(sequence           = s,
                         nucleic.acid.conc  = %s,
                         hybridisation.type = "dnadna",
                         Na.conc            = %s,
                         method.nn          = "%s",
                         correction.ion     = "%s",
                         size.threshold     = %s)
  as.numeric(r$Results[["Melting temperature (C)"]])
}, numeric(1), USE.NAMES = FALSE)
tD   <- proc.time()[["elapsed"]]
writeLines(format(tm, digits = 10), a[2])
tE   <- proc.time()[["elapsed"]]
tim <- sprintf(c("LOAD_SECONDS %%.4f", "READ_SECONDS %%.4f",
                 "COMPUTE_SECONDS %%.4f", "WRITE_SECONDS %%.4f"),
               c(tB - tA, tC - tB, tD - tC, tE - tD))
if (length(a) >= 3L) writeLines(tim, a[3])
cat(tim, sep = "\\n")
', format(conc, scientific = TRUE), format(na_m, scientific = FALSE),
   nn, ion, format(MELTING_SIZE_THRESHOLD, scientific = FALSE))

runner_rmelting <- make_runner_rmelting(MELTING_CONC, MELTING_ION, MELTING_NN)

make_runner_biopython <- function(saltcorr = SALT_BIOPYTHON,
                                  hi = DNAC_HIGH, lo = DNAC_LOW) sprintf('
import sys, time
tA = time.perf_counter()
from Bio.SeqUtils import MeltingTemp as mt
tB = time.perf_counter()
seqs = [l.strip() for l in open(sys.argv[1]) if l.strip()]
t0 = time.perf_counter()
tm = [mt.Tm_NN(s, nn_table=mt.DNA_NN4, saltcorr=%s,
               Na=%s, K=0, Tris=0, Mg=0, dNTPs=0,
               dnac1=%s, dnac2=%s, selfcomp=%s) for s in seqs]
t1 = time.perf_counter()
open(sys.argv[2], "w").write("\\n".join("%%.10f" %% v for v in tm) + "\\n")
t2 = time.perf_counter()
tim = ("LOAD_SECONDS %%.4f\\nREAD_SECONDS %%.4f\\n"
       "COMPUTE_SECONDS %%.4f\\nWRITE_SECONDS %%.4f\\n"
       %% (tB - tA, t0 - tB, t1 - t0, t2 - t1))
if len(sys.argv) > 3:
    open(sys.argv[3], "w").write(tim)
print(tim)
', saltcorr, NA_MM, hi, lo, ifelse(SELF_COMP, "True", "False"))

runner_biopython <- make_runner_biopython()

writeLines(runner_tmcalculator, file.path(outdir, "run_tmcalculator.R"))
writeLines(runner_rmelting,     file.path(outdir, "run_rmelting.R"))
writeLines(runner_biopython,    file.path(outdir, "run_biopython.py"))

tools <- list(
  TmCalculator = list(cmd = "Rscript", script = "run_tmcalculator.R",
                      available = TRUE,
                      engine = "R with a compiled C++ core"),
  rmelting     = list(cmd = "Rscript", script = "run_rmelting.R",
                      available = have_rmelting,
                      engine = "R calling the MELTING 5 Java engine"),
  Biopython    = list(cmd = python,    script = "run_biopython.py",
                      available = have_python,
                      engine = "pure Python")
)

## -- 3. One measured run ---------------------------------------------------
# Wall clock comes from the wrapper (it includes start-up); COMPUTE_SECONDS
# is printed by the runner itself and excludes it. Reporting both is what
# keeps the small input sizes interpretable: a JVM that takes two seconds to
# start would otherwise dominate the 1,000-sequence result.
parse_rss_bytes <- function(txt, unit) {
  ln <- grep("[Mm]aximum resident set size", txt, value = TRUE)
  if (!length(ln)) return(NA_real_)
  v <- suppressWarnings(as.numeric(gsub("[^0-9]", "", ln[1])))
  v * unit
}

run_once <- function(tool, n) {
  sf  <- seq_file(n)
  of  <- file.path(outdir, sprintf("tm_%s_%d.txt", tool, n))
  ef  <- tempfile()
  spec <- tools[[tool]]

  tf <- file.path(outdir, sprintf("time_%s_%d.txt", tool, n))
  unlink(c(of, tf))

  argv <- c(time_cmd$flag, spec$cmd, file.path(outdir, spec$script), sf, of, tf)
  if (nzchar(Sys.which("timeout")))
    argv <- c(time_cmd$flag, "timeout", timeout_s,
              spec$cmd, file.path(outdir, spec$script), sf, of, tf)

  t0  <- proc.time()[["elapsed"]]
  out <- suppressWarnings(
    system2(time_cmd$bin, argv, stdout = TRUE, stderr = ef))
  wall <- proc.time()[["elapsed"]] - t0
  err  <- readLines(ef, warn = FALSE)

  # Prefer the file the runner wrote; fall back to whatever reached stdout.
  src <- if (file.exists(tf)) readLines(tf, warn = FALSE) else c(out, err)
  getv <- function(key) {
    v <- suppressWarnings(as.numeric(
      sub(paste0(key, " "), "", grep(paste0("^", key), src, value = TRUE)[1])))
    if (length(v)) v else NA_real_
  }
  load_s  <- getv("LOAD_SECONDS")
  read_s  <- getv("READ_SECONDS")
  compute <- getv("COMPUTE_SECONDS")
  write_s <- getv("WRITE_SECONDS")

  # Four quantities, kept apart because they answer different questions.
  #   launch   process spawn and interpreter initialisation, before the runner
  #            takes its first timestamp
  #   startup  launch + attaching the package: what a user waits through before
  #            anything happens, and constant in the input
  #   io       reading the input file and writing the results. This exists only
  #            because each tool is run as a separate process so that peak
  #            memory can be attributed to it; a user passes data in memory.
  #            It scales with n, so folding it into startup (as this script did
  #            previously) made "start-up" grow with the input.
  #   compute  the calculation itself
  io_s    <- sum(c(read_s, write_s), na.rm = TRUE)
  known   <- sum(c(load_s, read_s, compute, write_s), na.rm = TRUE)
  launch  <- max(wall - known, 0)

  list(wall_s    = wall,
       compute_s = compute,
       load_s    = load_s,
       io_s      = if (all(is.na(c(read_s, write_s)))) NA_real_ else io_s,
       launch_s  = launch,
       startup_s = launch + (if (is.na(load_s)) 0 else load_s),
       rss_gb    = parse_rss_bytes(err, time_cmd$unit) / 1e9,
       ok        = file.exists(of) && length(readLines(of, warn = FALSE)) == n,
       out_file  = of)
}

## -- 4z. Recover dH and dS from each tool ----------------------------------
# Once the salt term is excluded, any remaining disagreement lives in the
# thermodynamics. None of the three tools reports dH and dS through its
# public interface, but both can be recovered from Tm measured at two
# strand concentrations, because
#
#   1/T1 - 1/T2 = R * ln(k1/k2) / (1000 * dH)
#
# so  dH = R * ln(k1/k2) / (1000 * (1/T1 - 1/T2))  and  dS follows.
#
# dH depends only on the RATIO of the two concentrations. It is therefore
# recoverable from MELTING 5 as well, even though we still do not know what
# its `nucleic.acid.conc` means: whatever the convention, a tenfold change
# of the input is a tenfold change of the internal k. Any offset in the
# convention is absorbed entirely by dS. That splits the remaining question
# in two: a dH disagreement is a table or summation problem, whereas dH
# agreeing while dS does not points at the initiation terms or at the
# concentration convention.
if (has_flag("--thermo")) {
  Rgas <- 1.987
  cal  <- file.path(outdir, "seqs_calibrate.txt")
  writeLines(head(readLines(seq_file(sizes[1])), 8), cal)
  seqs_cal <- readLines(cal)

  # Factor of ten between the two runs, salt off in every tool.
  runs <- list(
    TmCalculator = function(f) list(
      make_runner_tmcalculator("none", DNAC_HIGH,      DNAC_LOW),
      make_runner_tmcalculator("none", DNAC_HIGH * 10, DNAC_LOW * 10)),
    Biopython    = function(f) list(
      make_runner_biopython(0L, DNAC_HIGH,      DNAC_LOW),
      make_runner_biopython(0L, DNAC_HIGH * 10, DNAC_LOW * 10)),
    rmelting     = function(f) list(
      make_runner_rmelting(MELTING_CONC,      MELTING_ION, MELTING_NN, 1.0),
      make_runner_rmelting(MELTING_CONC * 10, MELTING_ION, MELTING_NN, 1.0))
  )
  ext <- c(TmCalculator = ".R", Biopython = ".py", rmelting = ".R")

  thermo <- list()
  for (tool in names(runs)) {
    if (!tools[[tool]]$available) next
    tm <- vector("list", 2L)
    for (j in 1:2) {
      rf <- file.path(outdir, paste0("run_thermo_", tool, "_", j, ext[[tool]]))
      writeLines(runs[[tool]](NULL)[[j]], rf)
      of <- file.path(outdir, sprintf("tm_thermo_%s_%d.txt", tool, j))
      unlink(of)
      system2(tools[[tool]]$cmd, c(rf, cal, of), stdout = NULL, stderr = NULL)
      tm[[j]] <- if (file.exists(of)) as.numeric(readLines(of, warn = FALSE))
                 else rep(NA_real_, length(seqs_cal))
    }
    T1 <- tm[[1]] + 273.15
    T2 <- tm[[2]] + 273.15
    k1 <- (DNAC_HIGH - DNAC_LOW / 2) * 1e-9         # ratio is all that matters
    dH <- Rgas * log(k1 / (k1 * 10)) / (1000 * (1 / T1 - 1 / T2))
    dS <- 1000 * dH / T1 - Rgas * log(k1)
    thermo[[tool]] <- data.frame(dH = dH, dS = dS)
  }

  cat("\n=== dH (kcal/mol) recovered from two concentrations ===\n")
  print(data.frame(sequence = seqs_cal,
                   lapply(thermo, `[[`, "dH"), check.names = FALSE),
        digits = 6, row.names = FALSE)
  cat("\n=== dS (cal/mol/K); an offset here can also be the concentration",
      "convention ===\n")
  print(data.frame(sequence = seqs_cal,
                   lapply(thermo, `[[`, "dS"), check.names = FALSE),
        digits = 6, row.names = FALSE)

  if (length(thermo) > 1L) {
    ref <- names(thermo)[1]
    cat("\nmax deviation vs ", ref, ":\n", sep = "")
    for (k in names(thermo)[-1])
      cat(sprintf("  %-14s dH %.4f kcal/mol   dS %.4f cal/mol/K\n", k,
                  max(abs(thermo[[k]]$dH - thermo[[ref]]$dH)),
                  max(abs(thermo[[k]]$dS - thermo[[ref]]$dS))))
  }
  quit(save = "no")
}

## -- 4a. Grid search over the MELTING 5 conventions -------------------------
# The R wrapper does not document what `nucleic.acid.conc` means, and the
# ion corrections are named by citation rather than by formula. Instead of
# guessing, run the combinations and let the data say which one matches.
# The reference is TmCalculator on the same eight sequences.
if (grid) {
  if (!have_rmelting) stop("rmelting is not installed.")

  cal <- file.path(outdir, "seqs_calibrate.txt")
  writeLines(head(readLines(seq_file(sizes[1])), 8), cal)

  ref_f <- file.path(outdir, "tm_TmCalculator_cal.txt")
  system2("Rscript", c(file.path(outdir, "run_tmcalculator.R"), cal, ref_f),
          stdout = NULL, stderr = NULL)
  ref <- as.numeric(readLines(ref_f, warn = FALSE))

  # Ct, Ct/2 and Ct/4 conventions around the Biopython effective value.
  concs <- c("Ct/4 = 12.5 nM" = 12.5e-9,
             "Ct/2 = 25 nM"   = 25e-9,
             "Ct   = 50 nM"   = 50e-9)
  # Only the plain additive log10 corrections; the Owczarzy and Tan entries
  # are GC- or Mg-dependent and are not what the other two tools are using.
  ions  <- c("san96", "san04", "schlif", "wet91", "ahs01")

  res <- expand.grid(conc = names(concs), ion = ions,
                     stringsAsFactors = FALSE)
  res$max_abs_dTm <- NA_real_
  res$mean_dTm    <- NA_real_

  tmp_runner <- file.path(outdir, "run_rmelting_grid.R")
  for (i in seq_len(nrow(res))) {
    writeLines(make_runner_rmelting(concs[[res$conc[i]]], res$ion[i],
                                    MELTING_NN), tmp_runner)
    of <- file.path(outdir, "tm_rmelting_grid.txt")
    unlink(of)
    ok <- tryCatch(
      system2("Rscript", c(tmp_runner, cal, of),
              stdout = NULL, stderr = NULL) == 0L,
      error = function(e) FALSE)
    if (ok && file.exists(of)) {
      d <- as.numeric(readLines(of, warn = FALSE)) - ref
      res$max_abs_dTm[i] <- max(abs(d))
      res$mean_dTm[i]    <- mean(d)
    }
    message(sprintf("  %-16s %-8s max|dTm| = %s", res$conc[i], res$ion[i],
                    ifelse(is.na(res$max_abs_dTm[i]), "failed",
                           sprintf("%.4f", res$max_abs_dTm[i]))))
  }

  res <- res[order(res$max_abs_dTm), ]
  cat("\n=== MELTING 5 convention grid, reference = TmCalculator ===\n")
  print(format(res, digits = 4), row.names = FALSE)
  cat("\nSet MELTING_CONC and MELTING_ION in the TOOL PARAMETERS block to the\n",
      "top row, then re-run --calibrate to confirm.\n", sep = "")
  quit(save = "no")
}

## -- 4b. Calibration -------------------------------------------------------
if (calibrate) {
  cat("\n=== rmelting::melting() signature, for the TODO above ===\n")
  if (have_rmelting) print(args(rmelting::melting)) else cat("(not installed)\n")

  cal <- file.path(outdir, "seqs_calibrate.txt")
  writeLines(head(readLines(seq_file(sizes[1])), 8), cal)

  cat("\n=== per-sequence Tm, all tools, identical input",
      if (nosalt) "(SALT CORRECTION OFF)" else
        sprintf("(salt: %s / %s / %s at %g mM Na)",
                SALT_TMCALC, SALT_BIOPYTHON, MELTING_ION, NA_MM),
      "===\n")
  vals <- list()
  for (tool in names(tools)) {
    if (!tools[[tool]]$available) next
    of <- file.path(outdir, sprintf("tm_%s_cal.txt", tool))
    system2(tools[[tool]]$cmd,
            c(file.path(outdir, tools[[tool]]$script), cal, of),
            stdout = NULL, stderr = NULL)
    if (file.exists(of))
      vals[[tool]] <- as.numeric(readLines(of, warn = FALSE))
  }
  cmp <- data.frame(sequence = readLines(cal), vals, check.names = FALSE)
  print(cmp, digits = 6, row.names = FALSE)

  if (length(vals) > 1L) {
    ref <- vals[[1]]
    cat("\nmax |dTm| vs ", names(vals)[1], ":\n", sep = "")
    for (k in names(vals)[-1])
      cat(sprintf("  %-14s %.4f C\n", k, max(abs(vals[[k]] - ref))))
    cat("\nIf any of these is not close to zero, the tools are NOT running\n",
        "the same model. Fix the TOOL PARAMETERS block before benchmarking.\n",
        sep = "")
  }
  quit(save = "no")
}

## -- 5. Sweep --------------------------------------------------------------
rows <- list()
skip <- character(0)

# --tools lets a slow tool be given its own, smaller set of sizes in a
# separate invocation writing into the same output directory. rmelting issues
# one Java call per sequence and has no batch interface, so pairing it with
# the sizes needed to characterise the other two would cost hours for a
# number that is already determined at n = 1000.
want <- argval("--tools", paste(names(tools), collapse = ","))
want <- trimws(strsplit(want, ",")[[1]])
if (!all(want %in% names(tools)))
  stop("unknown tool in --tools: ",
       paste(setdiff(want, names(tools)), collapse = ", "))

for (tool in want) {
  if (!tools[[tool]]$available) next
  for (n in sizes) {
    if (tool %in% skip) {
      message(sprintf("%-13s n=%-7d SKIPPED (exceeded %d s at a smaller n)",
                      tool, n, timeout_s))
      next
    }
    for (rep in seq_len(n_rep)) {
      message(sprintf("[rep %d] %-13s n = %d ...", rep, tool, n))
      r <- run_once(tool, n)
      rows[[length(rows) + 1L]] <- data.frame(
        tool = tool, engine = tools[[tool]]$engine, n = n, rep = rep,
        wall_s = r$wall_s, compute_s = r$compute_s,
        launch_s = r$launch_s, load_s = r$load_s,
        startup_s = r$startup_s, io_s = r$io_s,
        seqs_per_s = n / r$compute_s,
        rss_gb = r$rss_gb, ok = r$ok, stringsAsFactors = FALSE)
      if (!r$ok || r$wall_s > timeout_s) {
        message("  -> incomplete or over the time limit; ",
                "no larger input will be attempted for this tool")
        skip <- c(skip, tool)
        break
      }
    }
  }
}

S <- do.call(rbind, rows)

# Accumulate across invocations so that separate --tools runs build one table.
# Rows are keyed by tool/n/rep/seq_len; a repeat of the same configuration
# replaces the earlier measurement rather than being appended twice.
S$seq_len <- seq_len
# Stamp the package version. Merging is by configuration, not by run, so a
# CSV can otherwise end up holding measurements taken before and after a
# change to the package: the sizes no longer swept stay behind, and their
# start-up times belong to a build that no longer exists. Recording the
# version makes such a mixture visible instead of silently averaged.
S$pkg_version <- as.character(utils::packageVersion("TmCalculator"))
csv <- file.path(outdir, "crosstool_bench.csv")
n_new <- nrow(S)

if (has_flag("--fresh") && file.exists(csv)) {
  message("--fresh: discarding ", csv)
  file.remove(csv)
}

# Merge on the UNION of the columns, filling absences with NA. rbind() would
# demand identical column sets, so a CSV written by an earlier version of this
# script aborted the merge -- after the measurements had been taken and before
# anything was saved, which is the worst possible moment to fail.
if (file.exists(csv)) {
  prev <- tryCatch(utils::read.csv(csv, stringsAsFactors = FALSE),
                   error = function(e) NULL)
  if (!is.null(prev) && nrow(prev)) {
    if (!"seq_len" %in% names(prev)) prev$seq_len <- NA_integer_
    key  <- function(d) paste(d$tool, d$n, d$rep, d$seq_len, sep = "\r")
    prev <- prev[!key(prev) %in% key(S), , drop = FALSE]
    if (nrow(prev)) {
      all_cols <- union(names(prev), names(S))
      pad <- function(d) {
        for (k in setdiff(all_cols, names(d))) d[[k]] <- NA
        d[, all_cols, drop = FALSE]
      }
      S <- rbind(pad(prev), pad(S))
    }
  }
}
# Write to a temporary file first: a failure here must not destroy hours of
# measurements that are already in the old file.
tmpf <- paste0(csv, ".tmp")
utils::write.csv(S, tmpf, row.names = FALSE)
if (!file.rename(tmpf, csv)) stop("could not replace ", csv)
message(sprintf("wrote %d rows to %s (%d from this run)", nrow(S), csv, n_new))

## -- 6. Output consistency -------------------------------------------------
# Measured on the smallest input set, which every tool completes. The
# reference is TmCalculator, so a positive deviation means the other tool
# reports the higher Tm.
n_cons <- sizes[1]
cons <- NULL
ref_f <- file.path(outdir, sprintf("tm_TmCalculator_%d.txt", n_cons))
if (file.exists(ref_f)) {
  ref <- as.numeric(readLines(ref_f, warn = FALSE))
  cons <- do.call(rbind, lapply(setdiff(names(tools), "TmCalculator"), function(k) {
    f <- file.path(outdir, sprintf("tm_%s_%d.txt", k, n_cons))
    if (!file.exists(f)) return(NULL)
    v <- as.numeric(readLines(f, warn = FALSE))
    d <- v - ref
    data.frame(tool = k, n = n_cons,
               mean_dTm = mean(d), max_abs_dTm = max(abs(d)),
               pearson_r = stats::cor(v, ref), stringsAsFactors = FALSE)
  }))
  if (!is.null(cons))
    utils::write.csv(cons, file.path(outdir, "crosstool_consistency.csv"),
                     row.names = FALSE)
}

## -- 7. Report -------------------------------------------------------------
cat("\n== Output consistency (reference: TmCalculator, n =", n_cons, ") ==\n")
if (is.null(cons)) cat("(no comparable tool completed)\n") else
  print(format(cons, digits = 4), row.names = FALSE)

cat("\n== Throughput and memory, median of", n_rep, "repetitions ==\n")
usable <- S$ok & !is.na(S$compute_s)
if (any(usable)) {
  agg <- aggregate(cbind(compute_s, startup_s, io_s, seqs_per_s, rss_gb) ~
                     tool + n, data = S[usable, ], FUN = stats::median)
  print(format(agg[order(agg$tool, agg$n), ], digits = 4), row.names = FALSE)
} else cat("(no usable runs)\n")

# Runs are never dropped silently. A tool that failed, timed out, or did not
# report its own compute time still has to appear, because its absence would
# otherwise read as "not benchmarked" rather than "did not complete".
if (any(!usable)) {
  cat("\n== Runs excluded from the table above ==\n")
  bad <- S[!usable, c("tool", "n", "rep", "wall_s", "compute_s", "rss_gb", "ok")]
  bad$reason <- ifelse(!bad$ok, "wrong or missing output file",
                       "runner did not print COMPUTE_SECONDS")
  print(format(bad, digits = 4), row.names = FALSE)
  cat("\nFor these, wall_s still includes start-up and is a valid upper bound\n",
      "on the total cost. Inspect crosstool_bench.csv and the tm_*.txt files\n",
      "in the output directory to see how far the run got.\n", sep = "")
}

# Fixed per-call overhead dominates at small n: constructing the GRanges and
# converting five parameter tables for the compiled core cost the same
# whether one sequence is submitted or a million.
# Throughput measured below the size at which compute_s exceeds that overhead
# is a measurement of start-up, not of the calculation, and will understate a
# batch-oriented tool relative to a per-sequence one.
if (any(usable)) {
  flat <- do.call(rbind, lapply(split(S[usable, ], S$tool[usable]), function(d) {
    d <- d[order(d$n), ]
    data.frame(tool = d$tool[1],
               ratio_compute = max(d$compute_s) / min(d$compute_s),
               ratio_n       = max(d$n) / min(d$n),
               stringsAsFactors = FALSE)
  }))
  flat$overhead_bound <- flat$ratio_compute < 0.5 * flat$ratio_n
  if (any(flat$overhead_bound)) {
    cat("\n== Warning: these tools are start-up bound at the sizes tested ==\n")
    print(format(flat[flat$overhead_bound, ], digits = 3), row.names = FALSE)
    cat("Compute time barely grew with n, so seqs_per_s above reflects fixed\n",
        "overhead. Re-run with larger --sizes before quoting any throughput.\n",
        sep = "")
  }
}

cat("\n== Feature comparison (fill in from each tool's documentation) ==\n")
cat("Coordinate model / batch input / parallel interface / GRanges output\n")
cat("are capability differences and are NOT measured by this script.\n")

cat("\nWritten to ", normalizePath(outdir), "\n", sep = "")
cat("Environment: ", R.version.string, " on ",
    Sys.info()[["sysname"]], " ", Sys.info()[["release"]], "\n", sep = "")
