#!/usr/bin/env Rscript
# ===========================================================================
# make_table6.R -- Table 6: parallel strategy sweep on the human genome
#
#   Rscript inst/scripts/make_table6.R
#   Rscript inst/scripts/make_table6.R --csv results/bench_parallel_cluster.csv
#   Rscript inst/scripts/make_table6.R --common-baseline
#
# Reads the summary CSV written by bench_parallel_strategy.R or
# bench_parallel_cluster.R and emits three files in --outdir:
#
#   table6_parallel_strategy.csv   tidy numbers, one row per configuration
#   table6_parallel_strategy.tsv   same, rounded and formatted, for pasting
#   table6_parallel_strategy.md    same again as a Markdown table
#
# SPEEDUP IS RECOMPUTED HERE AND DOES NOT COME FROM THE CSV. The benchmark
# writes `speedup` as work_s / wall_compute_s, the sum of the task times over
# the wall clock. That ratio measures how well a configuration kept its own
# workers busy, which is a useful diagnostic but is not what the section
# asks. Segmenting changes how much total work there is -- each boundary adds
# a task and re-retrieves sequence -- so a strategy that does more work gets
# credited for it, and the three strategies cannot be compared. Here speedup
# is the serial wall clock divided by the configuration's wall clock, which
# is the quantity a user actually experiences, and efficiency follows from it
# as speedup / workers rather than from the task times.
#
# Start-up is excluded from the denominator throughout (wall_compute_s, not
# wall_s). Workers pay it concurrently, so it contributes its longest
# instance rather than its sum, and across 24 chromosomes it is a rounding
# error; leaving it in would understate every efficiency by the same amount
# and make the strategies look more similar than they are.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
has_flag <- function(f) f %in% args

csv    <- argval("--csv", system.file("extdata", "parallel_strategy_bench.csv",
                                      package = "TmCalculator"))
outdir <- argval("--outdir", "tables")
# Whether each strategy is measured against its own one-worker run or against
# a single number shared by all three. Per-strategy is the default because it
# is what bench_parallel_strategy.R prints, so the table agrees with the
# console output the numbers were first read from. --common-baseline uses the
# median of the one-worker runs instead, which is the fairer comparison when
# the three baselines disagree, and the printed spread says whether they do.
per_strategy <- !has_flag("--common-baseline")

if (!nzchar(csv) || !file.exists(csv))
  stop("benchmark summary not found. Run inst/scripts/bench_parallel_strategy.R ",
       "(or bench_parallel_cluster.R) and pass its CSV with --csv.")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

S <- utils::read.csv(csv, stringsAsFactors = FALSE)

need <- c("strategy", "n_workers", "n_tasks", "wall_s", "startup_s",
          "wall_compute_s", "work_s", "hard_floor", "max_rss_gb", "n_windows")
miss <- setdiff(need, names(S))
if (length(miss))
  stop("this summary predates the split of start-up from compute; re-run the ",
       "benchmark. Missing: ", paste(miss, collapse = ", "))

## -- Correctness gate ------------------------------------------------------
# Every configuration must have covered the same windows. Whole chromosomes
# are trimmed at the assembly gaps and segments are not, so a boundary or
# trimming mistake shows up here as a different window count -- and if the
# configurations did not compute the same thing, none of the times below can
# be compared with each other.
nw <- unique(S$n_windows)
if (length(nw) > 1L) {
  print(unique(S[, c("strategy", "n_workers", "n_windows")]))
  stop("configurations cover different numbers of windows (",
       paste(format(range(nw), big.mark = ","), collapse = " to "),
       "); the timings are not comparable.")
}

## -- Median over repetitions ----------------------------------------------
num <- c("wall_s", "startup_s", "wall_compute_s", "work_s", "hard_floor",
         "max_rss_gb", "idle_s")
num <- intersect(num, names(S))
grp <- list(strategy = S$strategy, n_workers = S$n_workers)
agg <- aggregate(S[, num], by = grp, FUN = stats::median)
# Task count is a property of the configuration, not something to take a
# median of; merged by key rather than assigned positionally, because two
# aggregate() calls are only in the same row order by coincidence.
agg <- merge(agg, aggregate(list(n_tasks = S$n_tasks), by = grp,
                            FUN = function(v) v[1]),
             by = c("strategy", "n_workers"))
n_rep <- max(table(paste(S$strategy, S$n_workers)))

## -- Speedup against the serial run ---------------------------------------
ser <- agg[agg$n_workers == 1L, c("strategy", "wall_compute_s")]
if (!nrow(ser))
  stop("the sweep contains no one-worker run, so there is no serial baseline ",
       "to divide by. Re-run with --workers 1,...")
spread <- (max(ser$wall_compute_s) - min(ser$wall_compute_s)) /
  stats::median(ser$wall_compute_s)

if (per_strategy) {
  names(ser)[2] <- "serial_s"
  agg <- merge(agg, ser, by = "strategy", all.x = TRUE)
} else {
  agg$serial_s <- stats::median(ser$wall_compute_s)
}
agg$speedup    <- agg$serial_s / agg$wall_compute_s
agg$efficiency <- agg$speedup / agg$n_workers
# How much more total work the configuration did than the serial run. This is
# the column that explains the turnover: past a few workers the wall clock
# stops improving not because the schedule is ragged but because the run has
# more work in it than it started with.
ser_work <- stats::setNames(agg$work_s[agg$n_workers == 1L],
                            agg$strategy[agg$n_workers == 1L])
agg$work_ratio <- agg$work_s / ser_work[agg$strategy]

ord <- c("static", "dynamic", "segment")
agg$strategy <- factor(agg$strategy, levels = intersect(ord, agg$strategy))
agg <- agg[order(agg$strategy, agg$n_workers), ]
agg$strategy <- as.character(agg$strategy)

## -- Emit ------------------------------------------------------------------
tidy <- data.frame(
  strategy       = agg$strategy,
  n_workers      = agg$n_workers,
  n_tasks        = agg$n_tasks,
  wall_compute_s = round(agg$wall_compute_s, 1),
  startup_s      = round(agg$startup_s, 1),
  work_s         = round(agg$work_s, 1),
  work_ratio     = round(agg$work_ratio, 2),
  longest_task_s = round(agg$hard_floor, 1),
  speedup        = round(agg$speedup, 2),
  efficiency     = round(agg$efficiency, 2),
  max_rss_gb     = round(agg$max_rss_gb, 2),
  stringsAsFactors = FALSE)

stem <- file.path(outdir, "table6_parallel_strategy")
utils::write.csv(tidy, paste0(stem, ".csv"), row.names = FALSE)

hdr <- c("Strategy", "Workers", "Tasks", "Wall time (s)", "Start-up (s)",
         "Total task time (s)", "Work ratio", "Longest task (s)", "Speedup",
         "Efficiency", "Peak RSS per worker (GB)")
body <- as.matrix(format(tidy, trim = TRUE, nsmall = 0))
# The strategy name is written once per block rather than on every row: the
# table is read down the worker column within a strategy, and repeating the
# name eleven times adds width without adding information.
body[, 1][duplicated(tidy$strategy)] <- ""

utils::write.table(rbind(hdr, body), paste0(stem, ".tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, col.names = FALSE)

md <- c(paste0("| ", paste(hdr, collapse = " | "), " |"),
        paste0("|", paste(rep("---", length(hdr)), collapse = "|"), "|"),
        apply(body, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, paste0(stem, ".md"))

## -- Report ----------------------------------------------------------------
cat("\nTable 6  (median of", n_rep, "repetition(s);",
    format(nw, big.mark = ","), "windows)\n\n")
print(`colnames<-`(body, hdr), quote = FALSE, right = TRUE)

cat("\nBaseline: ", if (per_strategy) "each strategy's own one-worker run" else
  sprintf("common, %.1f s (median of the one-worker runs)", agg$serial_s[1]), "\n")
if (spread > 0.05)
  cat(sprintf(
    "WARNING: the one-worker runs differ by %.0f%% across strategies. Some of\n%s\n%s\n",
    100 * spread,
    "  each speedup below is that difference rather than the scheduling. With",
    "  one repetition it is usually noise; confirm with --reps 3."))

best <- tidy[which.max(tidy$speedup), ]
cat(sprintf("\nBest configuration: %s at %d workers, %.0f s, speedup %.2f, %.2f GB per worker\n",
            best$strategy, best$n_workers, best$wall_compute_s, best$speedup,
            best$max_rss_gb))
cat("\nWritten:\n", paste0("  ", stem, c(".csv", ".tsv", ".md"), collapse = "\n"),
    "\n", sep = "")
