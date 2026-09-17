#!/usr/bin/env Rscript
# ===========================================================================
# Worker-count sweep for tm_profile() on a whole genome.
#
#   Rscript bench_tm_profile.R --outdir results --workers 1,2,3,4,5,6 --reps 3
#
# This measures the shipped function rather than a hand-rolled dispatch
# loop, which is the difference from bench_parallel_strategy.R: that script
# implemented three partitioning strategies to decide which one tm_profile()
# should use, and this one times the result.
#
# WHAT IS MEASURED. Wall-clock time of one tm_profile() call, including
# worker start-up, since that is what a user waits through. Peak memory per
# worker is not measured inside R: tm_profile() offers no hook inside its
# tasks, and instrumenting it for a benchmark would mean shipping a
# benchmark's needs in a user-facing function. It is sampled from outside
# instead, by the sampler in bench_tm_profile.lsf, and joined onto each run
# by timestamp at the end of this script. On a machine without that sampler
# the memory columns are NA and the timings are unaffected.
#
# START-UP IS SEPARATED. Every call starts its workers and stops them again,
# and a PSOCK worker attaches TmCalculator and the BSgenome for itself, which
# is on the order of ten seconds each. That is a real cost a user pays and it
# stays in wall_s, but it is also a constant, so on a short run it swamps the
# calculation and the table reads as though parallelism made things worse. It
# is therefore measured on its own, per worker count, by timing a call over a
# 200 kb region that does almost no work, and reported as startup_s with
# compute_s = wall_s - startup_s beside it.
#
# WHAT IS CHECKED. Every configuration must return the same number of
# windows and the same Tm sum. A worker count that changed the answer would
# otherwise show up as a speedup.
# ===========================================================================
suppressPackageStartupMessages({
  library(TmCalculator)
  library(BiocParallel)
})

## -- arguments --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("missing value for ", flag)
  args[[i + 1L]]
}
OUTDIR   <- getarg("--outdir", "results")
WORKERS  <- as.integer(strsplit(getarg("--workers", "1,2,3,4,5,6"), ",")[[1]])
REPS     <- as.integer(getarg("--reps", "3"))
WINDOW   <- as.integer(getarg("--window", "200"))
SLIDE    <- as.integer(getarg("--slide", as.character(WINDOW)))
SEGSIZE  <- as.numeric(getarg("--segsize", "50e6"))
UNIT     <- getarg("--unit", "segment")
PKG      <- getarg("--genome", "BSgenome.Hsapiens.UCSC.hg38")
NN       <- getarg("--nn-table", "DNA_NN_Breslauer_1986")
NA_MM    <- as.numeric(getarg("--na", "50"))
REGIONS  <- getarg("--regions", NA_character_)   # e.g. "chr21,chr22"; NA = all
MEMLOG   <- getarg("--memlog", file.path(OUTDIR, "mem_sampling.tsv"))
TMPDIR   <- getarg("--tmpdir", tempdir())

CALIBRATE <- !identical(getarg("--no-calibrate", "no"), "yes")
regions <- if (is.na(REGIONS)) NULL else strsplit(REGIONS, ",")[[1]]
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

for (p in c(PKG, "BiocParallel"))
  if (!requireNamespace(p, quietly = TRUE)) stop("missing package: ", p)

cat(sprintf(paste0("tm_profile sweep\n",
                   "  genome    : %s\n  regions   : %s\n",
                   "  window    : %d   slide: %d   unit: %s   segment: %.0f\n",
                   "  workers   : %s\n  reps      : %d\n",
                   "  TmCalculator %s, R %s\n\n"),
            PKG, if (is.null(regions)) "standard chromosomes" else paste(regions, collapse = ","),
            WINDOW, SLIDE, UNIT, SEGSIZE, paste(WORKERS, collapse = ","), REPS,
            as.character(packageVersion("TmCalculator")),
            paste0(R.version$major, ".", R.version$minor)))

## -- one configuration -------------------------------------------------------
# A fresh SnowParam per run. Reusing one cluster across configurations would
# leave the workers warm for every run after the first, which is exactly the
# start-up cost this benchmark is meant to include.
run_one <- function(n_workers) {
  gc(reset = TRUE, full = TRUE)
  bp <- if (n_workers <= 1L) NULL else SnowParam(workers = n_workers)
  t0 <- as.numeric(Sys.time())
  el <- system.time({
    prof <- tm_profile(PKG, regions = regions,
                       window = WINDOW, slide = SLIDE,
                       unit = UNIT, segment_size = SEGSIZE,
                       BPPARAM = bp, tmpdir = TMPDIR, verbose = FALSE,
                       method = "tm_nn", nn_table = NN, Na = NA_MM)
  })[["elapsed"]]
  t1 <- as.numeric(Sys.time())
  out <- list(
    wall_s     = el,
    t_start    = t0,
    t_end      = t1,
    n_windows  = length(prof),
    # A cheap fingerprint of the answer. Rounded because the reduction order
    # inside a worker is fixed but the order in which task results are
    # concatenated is not, and floating-point addition is not associative.
    tm_sum     = round(sum(prof$Tm), 3),
    # gc() reports the high-water mark in Mb in a column whose name
    # carries the unit; matched rather than indexed by position, since
    # the column order has changed between R versions.
    gc_peak_gb = sum(gc()[, grep("max used", colnames(gc()))[1]]) / 1024)
  rm(prof); gc(full = TRUE)
  out
}

## -- worker start-up, measured on its own -----------------------------------
# A 200 kb slice of the shortest requested sequence: enough to exercise the
# whole path, little enough that what is timed is almost entirely the cost of
# bringing the workers up. Taken from the middle, since chromosome ends are
# assembly gaps and would tile to nothing.
startup_of <- stats::setNames(rep(NA_real_, length(WORKERS)), as.character(WORKERS))
if (CALIBRATE) {
  sl  <- GenomeInfoDb::seqlengths(get(PKG, envir = asNamespace(PKG)))
  if (!is.null(regions)) sl <- sl[intersect(regions, names(sl))]
  if (!length(sl)) sl <- GenomeInfoDb::seqlengths(get(PKG, envir = asNamespace(PKG)))
  ch  <- names(sl)[which.min(sl)]
  mid <- floor(sl[[ch]] / 2)
  tiny <- sprintf("%s:%d-%d", ch, mid, mid + 2e5)
  cat(sprintf("calibrating start-up on %s\n", tiny))
  for (w in WORKERS) {
    bp <- if (w <= 1L) NULL else SnowParam(workers = w)
    el <- system.time(tm_profile(PKG, regions = tiny, window = WINDOW,
                                 slide = SLIDE, unit = UNIT,
                                 segment_size = SEGSIZE, BPPARAM = bp,
                                 tmpdir = TMPDIR, verbose = FALSE,
                                 method = "tm_nn", nn_table = NN,
                                 Na = NA_MM))[["elapsed"]]
    startup_of[as.character(w)] <- el
    cat(sprintf("  %d workers: %5.1f s\n", w, el))
  }
  cat("\n")
}

## -- sweep ------------------------------------------------------------------
# Repetitions outermost: a machine that drifts over the hours of a sweep
# then affects every worker count in the same way, and the drift shows up as
# a wide range rather than as a trend across the x axis.
rows <- list()
for (rep in seq_len(REPS)) {
  for (w in WORKERS) {
    cat(sprintf("rep %d/%d  workers %d ... ", rep, REPS, w)); flush.console()
    r <- run_one(w)
    cat(sprintf("%7.1f s  %s windows\n", r$wall_s,
                format(r$n_windows, big.mark = ",")))
    rows[[length(rows) + 1L]] <- data.frame(
      genome = PKG, unit = UNIT, segment_size = SEGSIZE,
      window = WINDOW, slide = SLIDE, nn_table = NN, na_mm = NA_MM,
      n_workers = w, rep = rep,
      wall_s = r$wall_s, t_start = r$t_start, t_end = r$t_end,
      n_windows = r$n_windows, tm_sum = r$tm_sum,
      gc_peak_gb = r$gc_peak_gb,
      startup_s = startup_of[[as.character(w)]],
      host = Sys.info()[["nodename"]],
      stringsAsFactors = FALSE)
  }
}
bench <- do.call(rbind, rows)
bench$compute_s <- bench$wall_s - bench$startup_s

## -- the answer must not depend on the worker count -------------------------
if (length(unique(bench$n_windows)) != 1L)
  stop("window count differs between configurations: ",
       paste(unique(bench$n_windows), collapse = ", "))
if (length(unique(bench$tm_sum)) != 1L)
  stop("Tm sum differs between configurations: ",
       paste(unique(bench$tm_sum), collapse = ", "))
cat(sprintf("\nall %d runs agree: %s windows, Tm sum %.3f\n",
            nrow(bench), format(bench$n_windows[1], big.mark = ","),
            bench$tm_sum[1]))

## -- memory, joined from the external sampler -------------------------------
# The sampler writes one line per sample for the whole job; each run claims
# the samples that fall inside its own interval. Sampling is periodic, so
# these are floors on the peak, and a run shorter than the sampling interval
# may catch no sample at all and is left NA rather than guessed at.
bench$peak_worker_gb <- NA_real_
bench$peak_job_gb    <- NA_real_
bench$n_mem_samples  <- 0L
if (file.exists(MEMLOG)) {
  mem <- utils::read.delim(MEMLOG, stringsAsFactors = FALSE)
  if (nrow(mem) && all(c("unix_time", "rss_mb") %in% names(mem))) {
    has_max1 <- "max1_mb" %in% names(mem)
    for (i in seq_len(nrow(bench))) {
      s <- mem[mem$unix_time >= bench$t_start[i] & mem$unix_time <= bench$t_end[i], ]
      bench$n_mem_samples[i] <- nrow(s)
      if (nrow(s)) {
        bench$peak_job_gb[i] <- max(s$rss_mb) / 1024
        if (has_max1) bench$peak_worker_gb[i] <- max(s$max1_mb) / 1024
      }
    }
    cat(sprintf("memory joined from %s (%d samples)\n", MEMLOG, nrow(mem)))
  }
} else {
  cat("no memory sampler log found; memory columns left NA\n")
}

## -- derived, and written ---------------------------------------------------
# split() orders its groups as character, so 10 would sort before 2; the
# result is reordered numerically before the serial baseline is taken.
agg <- do.call(rbind, lapply(split(bench, bench$n_workers), function(d)
  data.frame(n_workers = d$n_workers[1], n_rep = nrow(d),
             wall_s = median(d$wall_s), lo = min(d$wall_s), hi = max(d$wall_s),
             startup_s = d$startup_s[1],
             compute_s = median(d$compute_s),
             peak_worker_gb = suppressWarnings(max(d$peak_worker_gb, na.rm = TRUE)),
             stringsAsFactors = FALSE)))
agg <- agg[order(agg$n_workers), ]
agg$peak_worker_gb[!is.finite(agg$peak_worker_gb)] <- NA_real_
agg$speedup    <- agg$wall_s[agg$n_workers == min(agg$n_workers)] / agg$wall_s
agg$efficiency <- agg$speedup / agg$n_workers
# The same ratio with the constant removed from both sides. On a genome-scale
# run the two agree; on a short one they are the difference between "parallel
# made it worse" and "the job was shorter than the start-up".
agg$speedup_compute <- agg$compute_s[agg$n_workers == min(agg$n_workers)] /
  agg$compute_s

csv <- file.path(OUTDIR, "bench_tm_profile.csv")
utils::write.csv(bench, csv, row.names = FALSE)
utils::write.csv(agg, file.path(OUTDIR, "bench_tm_profile_summary.csv"),
                 row.names = FALSE)

cat("\n=== median over repetitions ===\n")
print(format(agg, digits = 4), row.names = FALSE)

# A run that is mostly start-up says nothing about how the function scales,
# and the table above would be read as though it did.
frac <- max(agg$startup_s / agg$wall_s, na.rm = TRUE)
if (is.finite(frac) && frac > 0.25)
  warning(sprintf(paste0("start-up is up to %.0f%% of wall time here, so these ",
                         "timings measure process launch more than the ",
                         "calculation.\n  Each PSOCK worker attaches ",
                         "TmCalculator and the BSgenome for itself, a fixed ",
                         "cost of roughly ten seconds.\n  Sweep a whole ",
                         "genome, or read speedup_compute rather than ",
                         "speedup."), 100 * frac), call. = FALSE)
cat(sprintf("\nwrote %s\n      %s\n", csv,
            file.path(OUTDIR, "bench_tm_profile_summary.csv")))
cat("\nsessionInfo:\n"); print(sessionInfo())
