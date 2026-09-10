#!/usr/bin/env Rscript
# ===========================================================================
# bench_parallel_strategy.R -- how the unit of work is chosen determines the
#                              speedup of a genome-wide TmCalculator run
#
# BACKGROUND
#
# Running one chromosome per worker with the BiocParallel default gives poor
# scaling, and the reason is scheduling, not the calculation. bplapply()
# splits X into as many static chunks as there are workers, in order. With 24
# chromosomes and 5 workers, worker 1 receives chromosomes 1-5 -- the five
# largest -- and worker 5 receives 21, 22, X, Y. Worker 1 then defines the
# makespan while worker 5 has long since finished and sits idle.
#
# This script measures the three strategies side by side, with correct
# per-task timing (the naive instrumentation records time since the worker
# process started, which is cumulative across the tasks in a chunk and is
# easily mistaken for a per-task cost):
#
#   "static"   one chromosome per task, BiocParallel default chunking
#   "dynamic"  one chromosome per task, tasks = length(X) so that the manager
#              hands out chromosomes one at a time; combined with a
#              largest-first ordering this is longest-processing-time-first
#              scheduling, which is within 4/3 of optimal
#   "segment"  fixed-size genomic segments, so tasks are near-equal by
#              construction and peak memory per worker is bounded by the
#              segment rather than by the largest chromosome
#
# Each task reports its own elapsed time and peak resident memory, so the run
# can be decomposed afterwards into
#   sum(task times)              serial-equivalent work
#   sum(task times) / n_workers  the floor achievable by perfect balancing
#   max(task time)               the hard floor, since one task cannot be split
#   wall clock                   what was actually achieved
# The gap between the wall clock and the balancing floor is scheduling loss;
# the gap between the balancing floor and max(task time) is what finer task
# granularity can still recover.
#
# HOW TO RUN
#   Rscript inst/scripts/bench_parallel_strategy.R
# TmCalculator must be INSTALLED (not load_all'ed): SnowParam workers load the
# installed package, and load_all() compiles the C++ core at -O0, which makes
# every timing here about six times too slow.
# ===========================================================================

## -- Configuration ---------------------------------------------------------
# A full sweep is 3 strategies x 4 worker counts x 3 repetitions, each pass
# covering all 24 chromosomes, and takes hours. Every setting is therefore
# overridable from the command line so that a smoke test can be run first:
#
#   Rscript inst/scripts/bench_parallel_strategy.R \
#     --chrs chr21,chr22 --workers 2,3 --reps 1 --outfile smoke.csv
#
# That exercises all three strategies end to end in a few minutes. Only when
# it completes and the window counts agree is it worth starting the real run.
args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}

pkg        <- argval("--pkg", "BSgenome.Hsapiens.UCSC.hg38")
chrs       <- strsplit(argval("--chrs", paste0("chr", c(1:22, "X", "Y"),
                                               collapse = ",")), ",")[[1]]
window     <- as.integer(argval("--window", "200"))
slide      <- as.integer(argval("--slide", as.character(window)))
nn_table   <- argval("--nn-table", "DNA_NN_Breslauer_1986")  # matches Section 2.2
Na         <- as.numeric(argval("--na", "50"))
seg_size   <- as.numeric(argval("--segsize", "50e6"))        # bp per task
# Physical cores, portably. Logical cores share execution units, so a sweep
# over them measures hyperthreading rather than parallel efficiency; the
# optimum is routinely reported one or two workers too high when they are
# counted. macOS exposes the physical count through sysctl and Linux through
# lscpu; parallel::detectCores(logical = FALSE) is the fallback and is not
# reliable on every platform, so it is reported when it is used.
.physical_cores <- function() {
  sysname <- Sys.info()[["sysname"]]
  n <- NA_integer_
  src <- "detectCores(logical = FALSE)"
  if (identical(sysname, "Darwin")) {
    v <- suppressWarnings(as.integer(system2("sysctl", c("-n", "hw.physicalcpu"),
                                             stdout = TRUE, stderr = FALSE)))
    if (length(v) && !is.na(v[1])) { n <- v[1]; src <- "sysctl hw.physicalcpu" }
  } else if (identical(sysname, "Linux")) {
    out <- suppressWarnings(system2("lscpu", stdout = TRUE, stderr = FALSE))
    pick <- function(key) {
      ln <- grep(key, out, value = TRUE)
      if (!length(ln)) return(NA_integer_)
      suppressWarnings(as.integer(sub(".*:\\s*", "", ln[1])))
    }
    cps <- pick("^Core\\(s\\) per socket:"); sk <- pick("^Socket\\(s\\):")
    if (!is.na(cps) && !is.na(sk)) { n <- cps * sk; src <- "lscpu" }
  }
  if (is.na(n) || n < 1L) n <- parallel::detectCores(logical = FALSE)
  if (is.na(n) || n < 1L) { n <- parallel::detectCores(); src <- "detectCores()" }
  list(n = as.integer(n), src = src)
}

phys <- .physical_cores()
# The sweep starts at one worker: it is the baseline every speedup is divided
# by. Note what that point actually measures. BiocParallel runs a single
# worker inside the manager process rather than spawning one, so no packages
# are re-loaded and nothing is serialised; the measured start-up is zero and
# the task runs in a session that is already warm. It is therefore a true
# serial baseline, not "the parallel machinery with one worker", and speedups
# computed against it include the cost of going parallel at all.
worker_set <- as.integer(strsplit(argval("--workers", "1,2,3,4,5,6"), ",")[[1]])
strategies <- strsplit(argval("--strategies", "static,dynamic,segment"), ",")[[1]]
n_rep      <- as.integer(argval("--reps", "3"))
                              # repeat each configuration: the same worker count
                              # has been observed to differ by ~28% between
                              # runs on this machine, so a single measurement
                              # must not be reported as a point value
outfile    <- argval("--outfile", "bench_parallel_strategy.csv")

message("platform   : ", Sys.info()[["sysname"]], ", ", phys$n,
        " physical cores (", phys$src, ")")
message("chromosomes: ", paste(chrs, collapse = ", "))
# The manager process is not idle: it receives and binds the per-task results.
# Asking for as many workers as there are physical cores therefore leaves it
# competing with them, and the sweep flattens or reverses at the top for that
# reason rather than for want of parallelism.
if (max(worker_set) >= phys$n)
  message("NOTE: ", max(worker_set), " workers requested on ", phys$n,
          " physical cores; the manager competes with the workers at the top ",
          "of this sweep.")
message("workers    : ", paste(worker_set, collapse = ", "),
        "   strategies: ", paste(strategies, collapse = ", "),
        "   reps: ", n_rep)

suppressPackageStartupMessages({
  library(TmCalculator)
  library(BiocParallel)
  library(GenomicRanges)
  library(GenomeInfoDb)
})
stopifnot(requireNamespace(pkg, quietly = TRUE),
          requireNamespace("ps", quietly = TRUE))

genome <- BSgenome::getBSgenome(pkg)
sl <- GenomeInfoDb::seqlengths(genome)[chrs]

## -- Task lists ------------------------------------------------------------
# Whole chromosomes, ordered largest first. The ordering is what makes
# dynamic dispatch effective: handing out the long tasks first leaves only
# short tasks to fill the ragged end of the run.
tasks_chrom <- lapply(names(sort(sl, decreasing = TRUE)), function(ch)
  list(chr = ch, start = 1L, end = as.integer(sl[[ch]])))

# Fixed-size segments. Boundaries are placed on multiples of `slide` so that
# the window grid is identical to the one a whole-chromosome run would
# produce; with slide = window (non-overlapping tiling) this makes the
# segmented result identical to the unsegmented one.
make_segments <- function(sl, seg_size, slide) {
  seg_size <- floor(seg_size / slide) * slide
  out <- list()
  for (ch in names(sl)) {
    st <- seq(1, sl[[ch]], by = seg_size)
    for (s in st)
      out[[length(out) + 1L]] <- list(chr = ch, start = as.integer(s),
                                      end = as.integer(min(s + seg_size - 1, sl[[ch]])))
  }
  # largest first, though segments are near-equal by construction
  out[order(vapply(out, function(z) z$end - z$start, numeric(1)), decreasing = TRUE)]
}
tasks_seg <- make_segments(sl, seg_size, slide)

message(sprintf("chromosome tasks: %d   segment tasks: %d",
                length(tasks_chrom), length(tasks_seg)))

## -- The unit of work ------------------------------------------------------
# Everything the task needs is derived from its own arguments; nothing is
# captured from the manager session, so the closure serialised to a SnowParam
# worker stays a few kilobytes rather than carrying the genome with it.
run_task <- function(task, pkg, window, slide, nn_table, Na, mode) {
  t0 <- proc.time()[["elapsed"]]           # per TASK, not per worker process

  # A PSOCK worker attaches the packages on its FIRST task only, and the
  # BSgenome package is several seconds to load. Timed inside the task, that
  # cost lands on whichever task happens to be dispatched to a cold worker,
  # where it can exceed the calculation by an order of magnitude: in a smoke
  # test an 0.8 Mb segment took 10.4 s on a cold worker and 0.6 s on a warm
  # one. Left in, it inflates sum(task times) and therefore the balancing
  # floor and the efficiency, which are the quantities this benchmark exists
  # to report. It is measured separately instead.
  suppressPackageStartupMessages({
    library(TmCalculator)
    library(pkg, character.only = TRUE)
  })
  t_load <- proc.time()[["elapsed"]]

  # trim_N: on a whole chromosome, trim the leading/trailing assembly gaps;
  # on a segment, do NOT trim, or segment boundaries would shift the window
  # grid relative to a whole-chromosome run and the two would not be
  # comparable. Interior N windows are dropped afterwards in both cases.
  bins <- make_genomiccoord(
    bsgenome = pkg, chromosomes = task$chr,
    window = window, slide = slide,
    start = task$start, end = task$end, strand = "+",
    trim_N = if (identical(mode, "segment")) "none" else "ends",
    verbose = FALSE)

  gr <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
                               method = "preload_chr")

  out <- tm_calculate(gr, method = "tm_nn",
                      nn_table = nn_table, Na = Na)$gr   # serial inside a task

  # Drop the sequence columns before returning. They are ~500 MB per large
  # chromosome and serialising them back to the manager would cost more than
  # the calculation; Tm and GC are what the profile consists of.
  out$sequence   <- NULL
  out$complement <- NULL

  attr(out, "bench") <- list(
    chr     = task$chr,
    start   = task$start,
    end     = task$end,
    bp      = task$end - task$start + 1,
    n_win   = length(out),
    secs      = proc.time()[["elapsed"]] - t_load,   # calculation only
    load_s    = t_load - t0,                        # 0 on a warm worker
    secs_wall = proc.time()[["elapsed"]] - t0,
    rss_gb  = ps::ps_memory_info()[["rss"]] / 1e9,
    pid     = Sys.getpid())
  out
}

## -- One configuration -----------------------------------------------------
run_config <- function(strategy, n_workers) {
  tasks <- if (strategy == "segment") tasks_seg else tasks_chrom
  mode  <- if (strategy == "segment") "segment" else "chromosome"

  # `tasks =` is the whole difference between static and dynamic. Left at its
  # default (0) BiocParallel pre-splits X into one chunk per worker; setting
  # it to length(X) makes the manager dispatch one element at a time, so a
  # worker that finishes early immediately takes the next chromosome.
  BPPARAM <- if (strategy == "static")
    SnowParam(workers = n_workers)
  else
    SnowParam(workers = n_workers, tasks = length(tasks))

  wall <- system.time({
    res <- bplapply(tasks, run_task,
                    pkg = pkg, window = window, slide = slide,
                    nn_table = nn_table, Na = Na, mode = mode,
                    BPPARAM = BPPARAM)
  })[["elapsed"]]

  b <- do.call(rbind, lapply(res, function(x)
    as.data.frame(attr(x, "bench"), stringsAsFactors = FALSE)))

  n_win_total <- sum(b$n_win)
  rm(res); invisible(gc())

  # Bringing the workers up is a one-off cost that the workers pay CONCURRENTLY,
  # so its contribution to the wall clock is the longest of them, not their
  # sum. It has to be separated out: `secs` now excludes package loading while
  # the measured wall clock still contains it, and mixing the two understates
  # efficiency by exactly this amount. On a two-chromosome smoke test it is
  # most of the run; across 24 chromosomes it is a rounding error, which is
  # itself worth being able to see rather than assume.
  startup <- if ("load_s" %in% names(b)) max(b$load_s, na.rm = TRUE) else 0
  wall_c  <- max(wall - startup, .Machine$double.eps)

  list(
    summary = data.frame(
      strategy      = strategy,
      n_workers     = n_workers,
      n_tasks       = nrow(b),
      wall_s        = wall,                      # as measured, incl. start-up
      startup_s     = startup,                   # concurrent worker start-up
      wall_compute_s = wall_c,                   # what the scheduling achieved
      work_s        = sum(b$secs),               # serial-equivalent total
      balance_floor = sum(b$secs) / n_workers,   # floor if perfectly balanced
      hard_floor    = max(b$secs),               # longest single task
      speedup       = sum(b$secs) / wall_c,
      efficiency    = sum(b$secs) / (n_workers * wall_c),
      idle_s        = n_workers * wall_c - sum(b$secs),
      max_rss_gb    = max(b$rss_gb),
      n_windows     = n_win_total,
      stringsAsFactors = FALSE),
    tasks = cbind(strategy = strategy, n_workers = n_workers, b))
}

## -- Sweep -----------------------------------------------------------------
sum_rows <- list(); task_rows <- list()
for (rep in seq_len(n_rep)) {
  for (st in strategies) {
    for (nw in worker_set) {
      message(sprintf("[rep %d] %-8s workers = %d ...", rep, st, nw))
      r <- run_config(st, nw)
      r$summary$rep <- rep
      sum_rows[[length(sum_rows) + 1L]]  <- r$summary
      task_rows[[length(task_rows) + 1L]] <- cbind(rep = rep, r$tasks)
      print(format(r$summary[, c("wall_s", "work_s", "balance_floor",
                                 "hard_floor", "efficiency", "idle_s",
                                 "max_rss_gb")], digits = 4),
            row.names = FALSE)
    }
  }
}

S <- do.call(rbind, sum_rows)
TASKS <- do.call(rbind, task_rows)
utils::write.csv(S, outfile, row.names = FALSE)
utils::write.csv(TASKS, sub("\\.csv$", "_tasks.csv", outfile), row.names = FALSE)

## -- Report ----------------------------------------------------------------
agg <- aggregate(cbind(wall_s, startup_s, wall_compute_s, work_s,
                       efficiency, max_rss_gb) ~
                   strategy + n_workers, data = S, FUN = mean)

# Speedup against the serial run of the SAME strategy. The per-configuration
# `speedup` column in the CSV divides by that configuration's own work_s, which
# makes it useless for comparing strategies: segmenting raises the total work
# it has to do, so it can report the higher ratio while finishing later. Here
# the denominator is a wall clock and the numerator is the serial wall clock,
# so the numbers are comparable across strategies and answer the question the
# section asks -- what does adding workers buy, and which unit of work is best?
base <- agg[agg$n_workers == 1L, c("strategy", "wall_compute_s")]
names(base)[2] <- "serial_s"
agg <- merge(agg, base, by = "strategy", all.x = TRUE)
agg$speedup <- agg$serial_s / agg$wall_compute_s
agg <- agg[order(agg$strategy, agg$n_workers), ]

# The serial run does the same work whichever way the tasks are cut, so the
# one-worker times should agree. Where they do not, the difference is
# machine noise, and it goes straight into the speedups because each
# strategy is divided by its own baseline. A strategy whose baseline came
# out slow then reports a flatteringly high speedup for reasons that have
# nothing to do with scheduling.
#
# Segmenting does add a little real work -- an extra task per boundary, and
# sequence retrieved once per segment -- but that is small, and separating
# it from the noise needs repetitions rather than a single pass.
if (length(unique(agg$strategy)) > 1L) {
  ser <- agg[agg$n_workers == 1L, c("strategy", "wall_compute_s", "work_s")]
  cat("\nSerial baselines (should agree; differences propagate into speedup)\n")
  print(format(ser, digits = 4), row.names = FALSE)
  spread <- (max(ser$wall_compute_s) - min(ser$wall_compute_s)) /
    stats::median(ser$wall_compute_s)
  if (spread > 0.05)
    cat(sprintf(
      "WARNING: serial baselines differ by %.0f%%. With one repetition this is\n%s\n",
      100 * spread,
      "noise, not a property of the strategies; do not compare speedups yet."))
}

cat("\nMean over", n_rep, "repetitions\n")
cat("--------------------------------------------------------------\n")
print(format(agg[, c("strategy", "n_workers", "wall_s", "startup_s",
                     "wall_compute_s", "work_s", "speedup", "efficiency",
                     "max_rss_gb")], digits = 4),
      row.names = FALSE)

# Consistency: the three strategies must produce the same profile, otherwise
# a speed comparison is meaningless.
cat("\nWindows produced (must be identical across strategies):\n")
print(unique(S[, c("strategy", "n_windows")]))

cat("\nWritten: ", outfile, " and ", sub("\\.csv$", "_tasks.csv", outfile),
    "\n", sep = "")
