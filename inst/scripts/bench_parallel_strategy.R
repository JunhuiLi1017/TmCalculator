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
pkg        <- "BSgenome.Hsapiens.UCSC.hg38"
chrs       <- paste0("chr", c(1:22, "X", "Y"))
window     <- 200L
slide      <- 200L
nn_table   <- "DNA_NN_SantaLucia_2004"
Na         <- 50
seg_size   <- 50e6            # segment mode: bp per task
worker_set <- c(3L, 4L, 5L, 6L)
strategies <- c("static", "dynamic", "segment")
n_rep      <- 2L              # repeat each configuration; run-to-run spread on
                              # a laptop is large enough that a single
                              # measurement should not be reported
outfile    <- "bench_parallel_strategy.csv"

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

  suppressPackageStartupMessages({
    library(TmCalculator)
    library(pkg, character.only = TRUE)
  })

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
    secs    = proc.time()[["elapsed"]] - t0,
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
  rm(res); invisible(base::gc())          # base:: -- the package exports gc()

  list(
    summary = data.frame(
      strategy      = strategy,
      n_workers     = n_workers,
      n_tasks       = nrow(b),
      wall_s        = wall,
      work_s        = sum(b$secs),               # serial-equivalent total
      balance_floor = sum(b$secs) / n_workers,   # floor if perfectly balanced
      hard_floor    = max(b$secs),               # longest single task
      efficiency    = sum(b$secs) / (n_workers * wall),
      idle_s        = n_workers * wall - sum(b$secs),
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
agg <- aggregate(cbind(wall_s, work_s, efficiency, max_rss_gb) ~
                   strategy + n_workers, data = S, FUN = mean)
agg$speedup <- agg$work_s / agg$wall_s
agg <- agg[order(agg$strategy, agg$n_workers), ]

cat("\nMean over", n_rep, "repetitions\n")
cat("--------------------------------------------------------------\n")
print(format(agg[, c("strategy", "n_workers", "wall_s", "speedup",
                     "efficiency", "max_rss_gb")], digits = 4),
      row.names = FALSE)

# Consistency: the three strategies must produce the same profile, otherwise
# a speed comparison is meaningless.
cat("\nWindows produced (must be identical across strategies):\n")
print(unique(S[, c("strategy", "n_windows")]))

cat("\nWritten: ", outfile, " and ", sub("\\.csv$", "_tasks.csv", outfile),
    "\n", sep = "")
