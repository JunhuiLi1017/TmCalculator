#!/usr/bin/env Rscript
# ===========================================================================
# benchmark_tools.R -- TmCalculator vs rmelting (MELTING 5 / melting5jars)
#                      vs Biopython Bio.SeqUtils.MeltingTemp
#
# Answers reviewer request for quantitative comparison of execution time,
# memory footprint and output consistency on an identical large-scale
# dataset: the E. coli K-12 MG1655 genome tiled into 23,208 non-overlapping
# 200 bp windows (the case study of the manuscript).
#
# DESIGN NOTES
#
# * Every tool runs in its OWN subprocess wrapped in /usr/bin/time, so peak
#   resident memory is measured the same way for R and for Python. Timing
#   the tools inside one R session would not give comparable memory numbers
#   and would let one tool's allocator state affect the next.
# * All tools receive the SAME input file (one sequence per line) and the
#   same thermodynamic conditions: SantaLucia 2004 nearest-neighbor
#   parameters, 50 mM Na+, 25 nM of each strand, Schildkraut/Lifson salt
#   correction.
# * Because rmelting calls the MELTING 5 Java program once per sequence
#   through rJava, it is orders of magnitude slower than the vectorised
#   tools. The script therefore runs a scaling series and, if the full
#   dataset would exceed `rmelting_max_n`, reports its full-dataset cost as
#   an extrapolation (clearly labelled as such) rather than leaving a blank.
# * The Tm values of every run are written to disk so that cross-tool
#   agreement can be quantified afterwards, not just speed.
#
# HOW TO RUN
#   1. Reboot or otherwise make sure the machine is idle and swap is clear.
#   2. Rscript inst/scripts/benchmark_tools.R
#      (from the package root; TmCalculator must be INSTALLED, not load_all'ed)
#
# OUTPUT
#   benchmark_tools_timing.csv       one row per tool x size
#   benchmark_tools_agreement.csv    pairwise Tm differences
#   plus a markdown summary printed to the console
# ===========================================================================
setwd("/Users/lij11/UMass Medical School Dropbox/Junhui Li/Project/UMMS/Github/JunhuiLi1017/TmCalculator")
## -- Configuration ---------------------------------------------------------
sizes         <- c(100L, 500L, 2000L, 5000L, 23208L)  # 23208 = full genome
rmelting_max_n <- 2000L      # raise if you are willing to wait
ecoli_pkg     <- "BSgenome.Ecoli.NCBI.ASM584v2"
chr_name      <- "U00096.3"
window        <- 200L
outdir        <- "benchmark_tools_out"

script_dir <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub("^--file=", "", a[grep("^--file=", a)][1]))
}, error = function(e) "inst/scripts")
if (is.na(script_dir) || !nzchar(script_dir)) script_dir <- "inst/scripts"

worker_R  <- file.path(script_dir, "bench_worker.R")
worker_py <- file.path(script_dir, "bench_worker.py")
stopifnot(file.exists(worker_R), file.exists(worker_py))
dir.create(outdir, showWarnings = FALSE)

## -- Availability of the comparison tools ----------------------------------
have_rmelting <- requireNamespace("rmelting", quietly = TRUE)
have_biopython <- system2("python3", c("-c", shQuote("import Bio")),
                          stdout = FALSE, stderr = FALSE) == 0L
message("rmelting available : ", have_rmelting)
message("Biopython available: ", have_biopython)
if (!have_rmelting)
  message("  install with: BiocManager::install('rmelting')  ",
          "(pulls in melting5jars and needs Java + rJava)")
if (!have_biopython)
  message("  install with: python3 -m pip install biopython")

## -- 1. Build the shared input --------------------------------------------
seqfile <- file.path(outdir, "ecoli_windows.txt")
if (!file.exists(seqfile)) {
  suppressPackageStartupMessages({
    library(TmCalculator)
    library(ecoli_pkg, character.only = TRUE)
  })
  genome  <- get("Ecoli", envir = asNamespace(ecoli_pkg))
  chr_len <- GenomeInfoDb::seqlengths(genome)[[chr_name]]
  bins <- make_genomiccoord(bsgenome = ecoli_pkg, chromosomes = chr_name,
                            window = window, slide = window,
                            start = 1, end = chr_len, strand = "+",
                            verbose = FALSE)
  gr <- to_genomic_ranges_fast(list(pkg_name = ecoli_pkg, seq = bins),
                               method = "preload_chr")
  writeLines(as.character(GenomicRanges::mcols(gr)$sequence), seqfile)
  message("wrote ", length(bins), " sequences to ", seqfile)
}
n_avail <- length(readLines(seqfile))
sizes <- sizes[sizes <= n_avail]
message("input sequences available: ", n_avail)

## -- 2. Timing helper: run a worker under /usr/bin/time --------------------
# macOS `time -l` reports "maximum resident set size" in bytes;
# GNU `time -v` reports "Maximum resident set size (kbytes)".
run_timed <- function(cmd, cmd_args, label) {
  is_mac <- Sys.info()[["sysname"]] == "Darwin"
  tf <- tempfile()
  full <- if (is_mac) c("-l", cmd, cmd_args) else c("-v", cmd, cmd_args)
  t_wall <- system.time(
    out <- suppressWarnings(
      system2("/usr/bin/time", full, stdout = TRUE, stderr = tf))
  )[["elapsed"]]
  err <- readLines(tf, warn = FALSE)
  unlink(tf)

  grab <- function(txt, pat) {
    hit <- grep(pat, txt, value = TRUE)
    if (!length(hit)) return(NA_real_)
    as.numeric(sub("^\\D*?([0-9]+).*$", "\\1", hit[1]))
  }
  rss <- if (is_mac) grab(err, "maximum resident set size") / 1e9 else
                     grab(err, "Maximum resident set size") / 1e6
  inner <- suppressWarnings(as.numeric(
    sub("BENCH_ELAPSED ", "", grep("^BENCH_ELAPSED", out, value = TRUE)[1])))
  n_na <- suppressWarnings(as.numeric(
    sub("BENCH_NA ", "", grep("^BENCH_NA", out, value = TRUE)[1])))
  prep <- suppressWarnings(as.numeric(
    sub("BENCH_PREP ", "", grep("^BENCH_PREP", out, value = TRUE)[1])))

  if (is.na(inner)) {
    message("  !! ", label, " failed; worker output:")
    message(paste("     ", utils::tail(c(out, err), 15), collapse = "\n"))
  }
  list(wall = t_wall, compute = inner, prep = prep, rss_gb = rss, n_na = n_na)
}

## -- 3. Scaling series -----------------------------------------------------
rows <- list()
tm_files <- list()

for (n in sizes) {
  # TmCalculator
  f <- file.path(outdir, sprintf("tm_tmcalculator_%d.csv", n))
  message(sprintf("[%6d] TmCalculator ...", n))
  r <- run_timed("Rscript", c(worker_R, "tmcalculator", seqfile, n, f),
                 "TmCalculator")
  rows[[length(rows) + 1]] <- data.frame(
    tool = "TmCalculator", n = n, compute_s = r$compute,
    prep_s = r$prep, total_s = r$wall, peak_rss_gb = r$rss_gb,
    n_failed = r$n_na, extrapolated = FALSE)
  tm_files[[paste0("TmCalculator_", n)]] <- f

  # Biopython
  if (have_biopython) {
    f <- file.path(outdir, sprintf("tm_biopython_%d.csv", n))
    message(sprintf("[%6d] Biopython ...", n))
    r <- run_timed("python3", c(worker_py, seqfile, n, f), "Biopython")
    rows[[length(rows) + 1]] <- data.frame(
      tool = "Bio.SeqUtils.MeltingTemp", n = n, compute_s = r$compute,
      prep_s = r$prep, total_s = r$wall, peak_rss_gb = r$rss_gb,
      n_failed = r$n_na, extrapolated = FALSE)
    tm_files[[paste0("Biopython_", n)]] <- f
  }

  # rmelting, two configurations (see bench_worker.R for why):
  #   defaults      -> approximative formula, because 200 bp > size.threshold
  #   forced NN     -> like-for-like nearest-neighbour comparison
  if (have_rmelting && n <= rmelting_max_n) {
    for (variant in c("rmelting", "rmelting_nn")) {
      lab <- if (variant == "rmelting") "rmelting (MELTING 5, defaults)" else
                                        "rmelting (MELTING 5, NN forced)"
      f <- file.path(outdir, sprintf("tm_%s_%d.csv", variant, n))
      message(sprintf("[%6d] %s ...", n, lab))
      r <- run_timed("Rscript", c(worker_R, variant, seqfile, n, f), lab)
      rows[[length(rows) + 1]] <- data.frame(
        tool = lab, n = n, compute_s = r$compute,
        prep_s = r$prep, total_s = r$wall, peak_rss_gb = r$rss_gb,
        n_failed = r$n_na, extrapolated = FALSE)
      tm_files[[paste0(variant, "_", n)]] <- f
    }
  }
}

timing <- do.call(rbind, rows)

## -- 4. Extrapolate rmelting to the full dataset if it was capped ----------
full_n <- max(sizes)
rm_rows <- timing[timing$tool == "rmelting (MELTING 5, NN forced)" &
                    !is.na(timing$compute_s), ]
if (nrow(rm_rows) >= 2 && !full_n %in% rm_rows$n) {
  # per-sequence cost is flat for a per-call Java interface; use the
  # largest measured size to extrapolate and label the row clearly
  per_seq <- max(rm_rows$compute_s) / max(rm_rows$n)
  timing <- rbind(timing, data.frame(
    tool = "rmelting (MELTING 5, NN forced)", n = full_n,
    compute_s = per_seq * full_n, prep_s = NA_real_, total_s = NA_real_,
    peak_rss_gb = max(rm_rows$peak_rss_gb, na.rm = TRUE),
    n_failed = NA_real_, extrapolated = TRUE))
}

utils::write.csv(timing, "benchmark_tools_timing.csv", row.names = FALSE)

## -- 5. Output consistency -------------------------------------------------
# Compare Tm on the largest size that every tool actually ran.
common_n <- min(vapply(split(timing$n[!timing$extrapolated],
                             timing$tool[!timing$extrapolated]),
                       max, numeric(1)))
read_tm <- function(key) {
  p <- tm_files[[key]]
  if (is.null(p) || !file.exists(p)) return(NULL)
  utils::read.csv(p)$Tm
}
agree <- NULL
a <- read_tm(paste0("TmCalculator_", common_n))
for (other in c("Biopython", "rmelting", "rmelting_nn")) {
  b <- read_tm(paste0(other, "_", common_n))
  if (is.null(a) || is.null(b)) next
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) next
  d <- a[ok] - b[ok]
  agree <- rbind(agree, data.frame(
    comparison = paste("TmCalculator vs", other), n = sum(ok),
    mean_diff_C = mean(d), median_diff_C = stats::median(d),
    sd_diff_C = stats::sd(d), max_abs_diff_C = max(abs(d)),
    pearson_r = stats::cor(a[ok], b[ok])))
}
if (!is.null(agree)) {
  utils::write.csv(agree, "benchmark_tools_agreement.csv", row.names = FALSE)
}

## -- 6. Report -------------------------------------------------------------
cat("\n\n### Runtime and memory\n\n")
cat("| tool | n sequences | Tm compute (s) | input prep (s) | peak RSS (GB) | note |\n")
cat("|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(timing))) {
  r <- timing[i, ]
  cat(sprintf("| %s | %s | %s | %s | %s | %s |\n", r$tool,
              format(r$n, big.mark = ","),
              ifelse(is.na(r$compute_s), "failed", sprintf("%.2f", r$compute_s)),
              ifelse(is.na(r$prep_s) || r$prep_s == 0, "-", sprintf("%.2f", r$prep_s)),
              ifelse(is.na(r$peak_rss_gb), "-", sprintf("%.2f", r$peak_rss_gb)),
              ifelse(r$extrapolated, "extrapolated from smaller sizes", "")))
}
if (!is.null(agree)) {
  cat(sprintf("\n### Output consistency (n = %d, identical conditions)\n\n",
              common_n))
  cat("| comparison | n | mean diff (C) | sd (C) | max |diff| (C) | r |\n")
  cat("|---|---|---|---|---|---|\n")
  for (i in seq_len(nrow(agree))) {
    r <- agree[i, ]
    cat(sprintf("| %s | %d | %.3f | %.3f | %.3f | %.5f |\n",
                r$comparison, r$n, r$mean_diff_C, r$sd_diff_C,
                r$max_abs_diff_C, r$pearson_r))
  }
  cat("\nResidual offsets are expected: the tools differ in documented",
      "\nconventions (e.g. how the strand concentration term is formed and",
      "\nwhich salt-correction variant is the default), not only in",
      "\nimplementation. A high correlation with a small, near-constant",
      "\noffset indicates agreement of the underlying model.\n")
}
cat("\nSaved: benchmark_tools_timing.csv, benchmark_tools_agreement.csv\n")
cat("Per-sequence Tm values are in ", outdir, "/\n", sep = "")
