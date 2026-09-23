#!/usr/bin/env Rscript
# ===========================================================================
# A short worker-count sweep on macOS, for checking that parallel execution
# behaves on this machine before committing to the full sweep.
#
#   Rscript inst/scripts/quick_parallel_mac.R            # chr1 + chr2, 1..6
#   WORKERS=1,2,4 REGIONS=chr21 Rscript ...              # smaller and faster
#
# This is not the benchmark. bench_tm_calculate_local.sh is: it repeats each
# configuration three times, separates worker start-up from compute time, and
# samples memory from outside the R process. This script runs each count once
# and reports the wall clock, which is enough to see whether the machine
# scales and whether it runs out of memory, and not enough to quote.
#
# WHY WALL CLOCK INCLUDES START-UP. Every call starts its workers and stops
# them again, and a PSOCK worker attaches TmCalculator and the BSgenome for
# itself, which costs on the order of ten seconds each. A user waits through
# that, so it stays in the number; it is also a constant, so on a short region
# it swamps the calculation and the table can read as though parallelism made
# things worse. Use a region big enough that the compute dominates, which is
# why the default is two whole chromosomes rather than chr21.
# ===========================================================================
suppressPackageStartupMessages({
  library(TmCalculator)
  library(BiocParallel)
  library(GenomicRanges)
})

PKG     <- Sys.getenv("GENOME",  "BSgenome.Hsapiens.UCSC.hg38")
REGIONS <- strsplit(Sys.getenv("REGIONS", "chr1,chr2"), ",")[[1]]
WORKERS <- as.integer(strsplit(Sys.getenv("WORKERS", "1,2,3,4,5,6"), ",")[[1]])
SEGSIZE <- as.numeric(Sys.getenv("SEGSIZE", "50e6"))

if (!requireNamespace(PKG, quietly = TRUE))
  stop("Install ", PKG, " first: BiocManager::install(\"", PKG, "\")")

cores <- tryCatch(as.integer(system("sysctl -n hw.physicalcpu", intern = TRUE)),
                  error = function(e) NA_integer_)
memgb <- tryCatch(round(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 2^30),
                  error = function(e) NA_real_)
cat(sprintf("machine   : %s physical cores, %s GB\n",
            ifelse(is.na(cores), "?", cores), ifelse(is.na(memgb), "?", memgb)))
cat(sprintf("genome    : %s\nregions   : %s\nsegment   : %.0f Mb\n\n",
            PKG, paste(REGIONS, collapse = ","), SEGSIZE / 1e6))
if (!is.na(cores) && max(WORKERS) > cores)
  message("NOTE: asking for more workers than physical cores. On this workload ",
          "the extra ones share a core and usually make the run slower.")

ARGS <- list(PKG, regions = REGIONS, window = 200, slide = 200,
             unit = "segment", segment_size = SEGSIZE,
             method = "tm_nn", nn_table = "DNA_NN_Breslauer_1986", Na = 50,
             verbose = FALSE)

ref <- NULL; rows <- list()
for (n in WORKERS) {
  gc(full = TRUE)
  bp <- if (n <= 1L) NULL else SnowParam(workers = n)
  el <- system.time(
    gr <- do.call(tm_calculate, c(ARGS, list(BPPARAM = bp)))$gr
  )[["elapsed"]]

  # Every configuration must return the same profile. Checked against the
  # first run rather than between neighbours, so a drift that accumulates
  # cannot hide.
  if (is.null(ref)) ref <- gr
  same <- identical(start(ref), start(gr)) && isTRUE(all.equal(ref$Tm, gr$Tm))
  if (!same) warning("worker count ", n, " returned a DIFFERENT profile")

  rows[[length(rows) + 1L]] <- data.frame(
    workers = n, wall_s = el, windows = length(gr),
    speedup = NA_real_, identical_to_serial = same)
  cat(sprintf("  %d worker%s : %7.1f s   %s windows   %s\n", n,
              if (n == 1L) " " else "s", el,
              format(length(gr), big.mark = ","),
              if (same) "same profile" else "DIFFERENT PROFILE"))
  rm(gr)
}

res <- do.call(rbind, rows)
res$speedup <- res$wall_s[res$workers == min(res$workers)] / res$wall_s
cat("\n")
print(res, row.names = FALSE, digits = 3)

best <- res$workers[which.min(res$wall_s)]
cat(sprintf("\nfastest   : %d workers, %.1f s (%.2fx)\n",
            best, min(res$wall_s), max(res$speedup)))
if (best < max(WORKERS))
  cat("Adding workers past that point did not help. On a laptop that is\n",
      "usually memory pressure or thermal throttling rather than anything in\n",
      "the software; watch Activity Monitor while it runs.\n", sep = "")
cat("\nFor numbers worth quoting, run inst/scripts/bench_tm_calculate_local.sh:\n",
    "three repetitions per count, start-up separated, memory sampled.\n", sep = "")
