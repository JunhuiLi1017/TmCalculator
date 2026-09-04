#!/usr/bin/env Rscript
# ===========================================================================
# bench_methods.R -- relative cost of the three Tm methods on identical input
#
# Motivation: the nearest-neighbor path runs in compiled code, whereas
# tm_gc() and tm_wallace() loop over sequences in R and call gc(), which
# splits each sequence into a character vector with s2c() and then makes
# five %in% passes over it. It is therefore possible that the "simple"
# composition-based methods are slower per window than the "expensive"
# nearest-neighbor one. This script measures that rather than assuming it,
# in either direction.
#
# Run with the INSTALLED package. devtools::load_all() compiles the C++ core
# at -O0, which would flatter the R paths by roughly a factor of six.
#
#   Rscript inst/scripts/bench_methods.R
# ===========================================================================

suppressPackageStartupMessages(library(TmCalculator))

set.seed(1)
w     <- 200L                        # window width, as in the case study
sizes <- c(2000L, 10000L, 50000L)    # ramp: stop early if a method blows up

## -- Input: one pool string, cut into non-overlapping windows --------------
# Building the pool once and slicing it is far cheaper than calling sample()
# per sequence, and keeps input construction out of the timed section.
make_seqs <- function(n) {
  pool   <- paste(sample(c("A", "C", "G", "T"), n * w, TRUE), collapse = "")
  starts <- seq(1L, by = w, length.out = n)
  substring(pool, starts, starts + w - 1L)
}

rows <- list()
for (n in sizes) {
  seqs <- make_seqs(n)

  # GRanges construction is shared by all three methods and is reported
  # separately so that it does not mask the difference between them.
  t_prep <- system.time(gr <- to_genomic_ranges(seqs))[["elapsed"]]

  for (m in c("tm_nn", "tm_gc", "tm_wallace")) {
    t <- tryCatch(
      system.time(suppressWarnings(tm_calculate(gr, method = m)))[["elapsed"]],
      error = function(e) { message("  ", m, " failed: ", conditionMessage(e)); NA_real_ })
    rows[[length(rows) + 1L]] <- data.frame(
      n = n, method = m, prep_s = t_prep, compute_s = t,
      us_per_window = 1e6 * t / n, stringsAsFactors = FALSE)
    message(sprintf("[%6d] %-11s %7.2f s  (%6.1f us/window)", n, m, t, 1e6 * t / n))
  }
}

res <- do.call(rbind, rows)
cat("\n")
print(format(res, digits = 4), row.names = FALSE)

## -- The comparison the script exists to make ------------------------------
wide <- reshape(res[, c("n", "method", "compute_s")], idvar = "n",
                timevar = "method", direction = "wide")
names(wide) <- sub("^compute_s\\.", "", names(wide))
wide$gc_vs_nn      <- wide$tm_gc / wide$tm_nn
wide$wallace_vs_nn <- wide$tm_wallace / wide$tm_nn
cat("\nCompute time relative to tm_nn (values > 1 mean the composition-based\n",
    "method is slower than the nearest-neighbor one):\n", sep = "")
print(format(wide, digits = 4), row.names = FALSE)

utils::write.csv(res, "bench_methods.csv", row.names = FALSE)
cat("\nWritten to bench_methods.csv\n")
