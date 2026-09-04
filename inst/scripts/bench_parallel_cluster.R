#!/usr/bin/env Rscript
# ===========================================================================
# bench_parallel_cluster.R -- parallel-strategy scaling benchmark, run on a
#                             dedicated compute node
#
# WHY THIS EXISTS SEPARATELY FROM bench_parallel_strategy.R
#
# The laptop version of this benchmark is not reproducible: repeated runs of
# the same configuration differed by nearly 30%, because whatever else the
# machine was doing at the time competed for cores and memory. A scaling
# curve measured that way cannot distinguish a real effect of worker count
# from background load. This version is written for an exclusively allocated
# compute node, and it removes the three things that otherwise make cluster
# timings *worse* than laptop timings:
#
#   1. Shared filesystem. BSgenome reads the .2bit file by random access. If
#      the package lives on NFS/GPFS, several workers reading it at once turn
#      sequence retrieval into a measurement of filesystem latency. The
#      genome package is therefore staged to node-local scratch before any
#      timing starts, and the workers are pointed at that copy.
#   2. Core count. detectCores() reports the whole node, not the allocation.
#      Oversubscribing a cgroup-limited job produces meaningless numbers, so
#      the core count is read from the scheduler.
#   3. Hyperthreading. Worker counts are swept over PHYSICAL cores; logical
#      cores share execution units and inflate the apparent optimum.
#
# Every configuration is repeated, and the report gives the median and the
# full range, so that the spread is visible rather than hidden behind a
# single number.
#
# HOW TO RUN
#   bsub < inst/scripts/bench_parallel_cluster.lsf
# or interactively on an exclusively allocated node:
#   Rscript inst/scripts/bench_parallel_cluster.R --outdir results
#
# TmCalculator must be INSTALLED, not load_all'ed: workers load the installed
# package, and load_all() compiles the C++ core at -O0, which makes every
# timing here roughly six times too slow.
# ===========================================================================

## -- Command line ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
outdir    <- argval("--outdir", "bench_cluster_out")
n_rep     <- as.integer(argval("--reps", "3"))
scratch   <- argval("--scratch", "")
seg_size  <- as.numeric(argval("--segsize", "50e6"))
stage_lib <- !identical(argval("--stage", "yes"), "no")

## -- Configuration ---------------------------------------------------------
pkg        <- "BSgenome.Hsapiens.UCSC.hg38"
chrs       <- paste0("chr", c(1:22, "X", "Y"))
window     <- 200L
slide      <- 200L
nn_table   <- "DNA_NN_SantaLucia_2004"
Na         <- 50
strategies <- c("static", "dynamic", "segment")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## -- 1. How many cores do we actually have? --------------------------------
# LSF sets LSB_DJOB_NUMPROC to the number of slots granted. SLURM's
# equivalent is read too, so the script also works if the cluster migrates.
# Falling back to detectCores() is a last resort and is reported as such,
# because on a shared node it is almost certainly wrong.
alloc_cores <- local({
  for (v in c("LSB_DJOB_NUMPROC", "SLURM_CPUS_ON_NODE", "NSLOTS")) {
    x <- suppressWarnings(as.integer(Sys.getenv(v, "")))
    if (!is.na(x) && x > 0L) return(list(n = x, src = v))
  }
  list(n = parallel::detectCores(), src = "detectCores() [UNRELIABLE]")
})

# Physical cores: logical cores share execution units, so a sweep over
# logical cores measures hyperthreading, not parallel efficiency.
phys_cores <- local({
  x <- suppressWarnings(as.integer(
    system2("lscpu", stdout = TRUE, stderr = FALSE) |>
      grep(pattern = "^Core\\(s\\) per socket:", value = TRUE) |>
      sub(pattern = ".*:\\s*", replacement = "")))
  y <- suppressWarnings(as.integer(
    system2("lscpu", stdout = TRUE, stderr = FALSE) |>
      grep(pattern = "^Socket\\(s\\):", value = TRUE) |>
      sub(pattern = ".*:\\s*", replacement = "")))
  if (length(x) && length(y) && !is.na(x[1]) && !is.na(y[1])) x[1] * y[1] else NA_integer_
})

usable <- if (!is.na(phys_cores)) min(alloc_cores$n, phys_cores) else alloc_cores$n
# Leave one core for the manager process, which is not idle: it receives and
# binds the per-task results.
worker_set <- unique(pmax(2L, seq.int(2L, max(2L, usable - 1L),
                                      length.out = min(6L, max(2L, usable - 1L)))))
worker_set <- sort(unique(as.integer(round(worker_set))))

cpu_model <- local({
  x <- grep("^Model name:", system2("lscpu", stdout = TRUE, stderr = FALSE),
            value = TRUE)
  if (length(x)) trimws(sub(".*:", "", x[1])) else NA_character_
})

env_info <- data.frame(
  host        = Sys.info()[["nodename"]],
  cpu_model   = cpu_model,
  alloc_cores = alloc_cores$n,
  cores_src   = alloc_cores$src,
  phys_cores  = phys_cores,
  worker_set  = paste(worker_set, collapse = ","),
  r_version   = paste0(R.version$major, ".", R.version$minor),
  stringsAsFactors = FALSE)
message("=== environment ===")
print(t(env_info))

## -- 2. Stage the genome package to node-local scratch ---------------------
# This is the single most important step. Left on a shared filesystem, the
# sequence-retrieval phase measures storage contention rather than the
# package.
if (!nzchar(scratch)) {
  for (v in c("TMPDIR", "LSB_JOBTMPDIR", "SCRATCH")) {
    x <- Sys.getenv(v, "")
    if (nzchar(x) && dir.exists(x)) { scratch <- x; break }
  }
  if (!nzchar(scratch)) scratch <- tempdir()
}
message("node-local scratch: ", scratch)

lib_local <- file.path(scratch, "Rlib_bench")
if (stage_lib) {
  dir.create(lib_local, showWarnings = FALSE, recursive = TRUE)
  src <- find.package(pkg)                      # errors early if not installed
  if (!dir.exists(file.path(lib_local, pkg))) {
    message("staging ", pkg, " (",
            round(sum(file.size(list.files(src, recursive = TRUE,
                                           full.names = TRUE)), na.rm = TRUE)/1e9, 2),
            " GB) to node-local scratch ...")
    t_stage <- system.time(
      ok <- file.copy(src, lib_local, recursive = TRUE))[["elapsed"]]
    stopifnot(ok)
    message(sprintf("  staged in %.1f s", t_stage))
  }
  .libPaths(c(lib_local, .libPaths()))
  # SnowParam starts fresh R processes; they inherit this environment
  # variable, and run_task() also sets .libPaths() defensively.
  Sys.setenv(R_LIBS = paste(c(lib_local, .libPaths()), collapse = .Platform$path.sep))
} else {
  lib_local <- .libPaths()[1]
}

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

## -- 3. Task lists ---------------------------------------------------------
# Largest first. With dynamic dispatch this is longest-processing-time-first
# scheduling, whose makespan is provably within 4/3 of optimal; with the
# default static chunking it is the worst possible order, which is exactly
# the effect being measured.
tasks_chrom <- lapply(names(sort(sl, decreasing = TRUE)), function(ch)
  list(chr = ch, start = 1L, end = as.integer(sl[[ch]])))

make_segments <- function(sl, seg_size, slide) {
  seg_size <- floor(seg_size / slide) * slide   # keep the window grid aligned
  out <- list()
  for (ch in names(sl)) {
    for (s in seq(1, sl[[ch]], by = seg_size))
      out[[length(out) + 1L]] <- list(
        chr = ch, start = as.integer(s),
        end = as.integer(min(s + seg_size - 1, sl[[ch]])))
  }
  out[order(vapply(out, function(z) z$end - z$start, numeric(1)),
            decreasing = TRUE)]
}
tasks_seg <- make_segments(sl, seg_size, slide)
message(sprintf("chromosome tasks: %d   segment tasks: %d (%.0f Mb each)",
                length(tasks_chrom), length(tasks_seg), seg_size / 1e6))

## -- 4. The unit of work ---------------------------------------------------
run_task <- function(task, pkg, window, slide, nn_table, Na, mode, lib) {
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  t0 <- proc.time()[["elapsed"]]                # per TASK, not per process

  suppressPackageStartupMessages({
    library(TmCalculator)
    library(pkg, character.only = TRUE)
  })

  # On a whole chromosome, trim the leading/trailing assembly gaps. On a
  # segment, do not: trimming would shift the window grid relative to a
  # whole-chromosome run and the two would no longer be comparable. Interior
  # all-N windows are dropped downstream in both cases.
  bins <- make_genomiccoord(
    bsgenome = pkg, chromosomes = task$chr,
    window = window, slide = slide,
    start = task$start, end = task$end, strand = "+",
    trim_N = if (identical(mode, "segment")) "none" else "ends",
    verbose = FALSE)

  gr  <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
                                method = "preload_chr")
  out <- tm_calculate(gr, method = "tm_nn",
                      nn_table = nn_table, Na = Na)$gr   # serial within a task

  # ~500 MB per large chromosome; returning it would cost more than the
  # calculation. Tm and GC are the profile.
  out$sequence <- NULL; out$complement <- NULL

  attr(out, "bench") <- list(
    chr = task$chr, start = task$start, end = task$end,
    bp = task$end - task$start + 1, n_win = length(out),
    secs = proc.time()[["elapsed"]] - t0,
    rss_gb = ps::ps_memory_info()[["rss"]] / 1e9,
    pid = Sys.getpid())
  out
}

## -- 5. One configuration --------------------------------------------------
run_config <- function(strategy, n_workers) {
  tasks <- if (strategy == "segment") tasks_seg else tasks_chrom
  mode  <- if (strategy == "segment") "segment" else "chromosome"

  # `tasks =` is the entire difference between static and dynamic. At its
  # default (0), bplapply pre-splits X into one contiguous chunk per worker;
  # set to length(X) the manager hands out one element at a time, so a worker
  # that finishes early immediately takes the next piece of work.
  BPPARAM <- if (strategy == "static")
    SnowParam(workers = n_workers)
  else
    SnowParam(workers = n_workers, tasks = length(tasks))

  wall <- system.time({
    res <- bplapply(tasks, run_task, pkg = pkg, window = window, slide = slide,
                    nn_table = nn_table, Na = Na, mode = mode, lib = lib_local,
                    BPPARAM = BPPARAM)
  })[["elapsed"]]

  b <- do.call(rbind, lapply(res, function(x)
    as.data.frame(attr(x, "bench"), stringsAsFactors = FALSE)))
  n_win_total <- sum(b$n_win)
  rm(res); invisible(base::gc())        # base:: -- the package exports gc()

  list(summary = data.frame(
         strategy = strategy, n_workers = n_workers, n_tasks = nrow(b),
         wall_s = wall,
         work_s = sum(b$secs),                    # serial-equivalent total
         balance_floor = sum(b$secs) / n_workers, # floor if perfectly balanced
         hard_floor = max(b$secs),                # longest indivisible task
         efficiency = sum(b$secs) / (n_workers * wall),
         idle_s = n_workers * wall - sum(b$secs),
         max_rss_gb = max(b$rss_gb), n_windows = n_win_total,
         stringsAsFactors = FALSE),
       tasks = cbind(strategy = strategy, n_workers = n_workers, b))
}

## -- 6. Sweep --------------------------------------------------------------
sum_rows <- list(); task_rows <- list()
for (rep in seq_len(n_rep)) {
  for (st in strategies) for (nw in worker_set) {
    message(sprintf("[rep %d/%d] %-8s workers=%2d ...", rep, n_rep, st, nw))
    r <- run_config(st, nw)
    r$summary$rep <- rep
    sum_rows[[length(sum_rows) + 1L]]   <- r$summary
    task_rows[[length(task_rows) + 1L]] <- cbind(rep = rep, r$tasks)
    print(format(r$summary[, c("wall_s", "work_s", "balance_floor",
                               "hard_floor", "efficiency", "max_rss_gb")],
                 digits = 4), row.names = FALSE)
  }
}

S     <- do.call(rbind, sum_rows)
TASKS <- do.call(rbind, task_rows)
S <- cbind(S, env_info[rep(1, nrow(S)), c("host", "cpu_model")], row.names = NULL)
utils::write.csv(S, file.path(outdir, "bench_parallel_cluster.csv"), row.names = FALSE)
utils::write.csv(TASKS, file.path(outdir, "bench_parallel_cluster_tasks.csv"), row.names = FALSE)
utils::write.csv(env_info, file.path(outdir, "bench_environment.csv"), row.names = FALSE)

## -- 7. Report: median and range, never a single number --------------------
agg <- do.call(rbind, lapply(split(S, list(S$strategy, S$n_workers), drop = TRUE),
  function(g) data.frame(
    strategy = g$strategy[1], n_workers = g$n_workers[1], n_rep = nrow(g),
    wall_med = median(g$wall_s), wall_min = min(g$wall_s), wall_max = max(g$wall_s),
    spread_pct = 100 * (max(g$wall_s) - min(g$wall_s)) / median(g$wall_s),
    speedup = median(g$work_s) / median(g$wall_s),
    efficiency = median(g$efficiency), max_rss_gb = max(g$max_rss_gb),
    stringsAsFactors = FALSE)))
agg <- agg[order(agg$strategy, agg$n_workers), ]

cat("\nMedian over", n_rep, "repetitions (spread_pct = range / median)\n")
cat("---------------------------------------------------------------\n")
print(format(agg, digits = 4), row.names = FALSE)

if (any(agg$spread_pct > 10))
  cat("\nWARNING: some configurations vary by >10% across repetitions.\n",
      "The node was probably not idle; do not report these numbers.\n", sep = "")

cat("\nWindows produced (must be identical across all strategies):\n")
print(unique(S[, c("strategy", "n_windows")]))

cat("\nWritten to ", normalizePath(outdir), "\n", sep = "")
