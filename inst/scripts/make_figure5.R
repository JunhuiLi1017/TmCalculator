#!/usr/bin/env Rscript
# ===========================================================================
# make_figure5.R -- publication-resolution export of Figure 5
#
#   (A) speedup against worker count, three task-partitioning strategies
#   (B) total task time, showing how much work each configuration performed
#   (C) peak resident set size of the heaviest worker
#
#   Rscript inst/scripts/make_figure5.R
#   Rscript inst/scripts/make_figure5.R --csv results/bench_parallel_cluster.csv
#   Rscript inst/scripts/make_figure5.R --common-baseline
#
# The three panels are one argument, not three measurements. Panel A shows
# that every strategy turns over well before the cores run out, which invites
# the usual explanation -- a ragged schedule, one long task holding up the
# end. Panel B rules that out: the total task time grows with the worker
# count, so the later configurations are not dividing a fixed amount of work
# badly, they are doing more of it. Panel C gives the reason and the remedy
# in the same picture, since the strategy that holds its memory down is the
# one whose work grows least.
#
# Speedup is recomputed from the wall clocks rather than read from the CSV's
# `speedup` column, which is work_s / wall_compute_s. That ratio rewards a
# strategy for the extra work segmenting creates and so cannot compare the
# three. See make_table6.R for the same reasoning at more length.
#
# Panel A carries a dashed line of slope one. Without it a curve reaching 3.2
# looks like a good result on its own terms; against the line it is visibly a
# third of what the hardware was asked for, which is the point of the panel.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
has_flag <- function(f) f %in% args

csv    <- argval("--csv", system.file("extdata", "parallel_strategy_bench.csv",
                                      package = "TmCalculator"))
outdir <- argval("--outdir", "figures")
dpi    <- as.numeric(argval("--dpi", "600"))
fig_w  <- as.numeric(argval("--width",  "10.5"))   # three panels, MDPI full width
fig_h  <- as.numeric(argval("--height", "3.9"))
per_strategy <- !has_flag("--common-baseline")

if (!nzchar(csv) || !file.exists(csv))
  stop("benchmark summary not found. Run inst/scripts/bench_parallel_strategy.R ",
       "(or bench_parallel_cluster.R) and pass its CSV with --csv.")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

S <- utils::read.csv(csv, stringsAsFactors = FALSE)
need <- c("strategy", "n_workers", "wall_compute_s", "work_s", "max_rss_gb",
          "n_windows")
miss <- setdiff(need, names(S))
if (length(miss))
  stop("this summary predates the split of start-up from compute; re-run the ",
       "benchmark. Missing: ", paste(miss, collapse = ", "))

# Configurations that did not cover the same windows are not comparable, and
# the difference is a boundary or trimming mistake rather than a measurement.
if (length(unique(S$n_windows)) > 1L)
  stop("configurations cover different numbers of windows; the timings are ",
       "not comparable. Inspect n_windows in ", csv)

## -- Median and range over repetitions -------------------------------------
summ <- function(v) c(med = stats::median(v), lo = min(v), hi = max(v))
agg <- do.call(rbind, lapply(
  split(S, list(S$strategy, S$n_workers), drop = TRUE), function(d) {
    w <- summ(d$wall_compute_s); k <- summ(d$work_s); m <- summ(d$max_rss_gb)
    data.frame(strategy = d$strategy[1], n_workers = d$n_workers[1],
               wall_med = w[["med"]], wall_lo = w[["lo"]], wall_hi = w[["hi"]],
               work_med = k[["med"]], work_lo = k[["lo"]], work_hi = k[["hi"]],
               rss_med  = m[["med"]], rss_lo  = m[["lo"]], rss_hi  = m[["hi"]],
               stringsAsFactors = FALSE)
  }))

ser <- agg[agg$n_workers == 1L, c("strategy", "wall_med")]
if (!nrow(ser))
  stop("the sweep contains no one-worker run, so there is no serial baseline ",
       "to divide by. Re-run with --workers 1,...")
if (per_strategy) {
  base <- stats::setNames(ser$wall_med, ser$strategy)[agg$strategy]
} else {
  base <- stats::median(ser$wall_med)
}
# The range bars follow the same division: a longer wall clock is a smaller
# speedup, so the low and high ends swap.
agg$sp_med <- base / agg$wall_med
agg$sp_lo  <- base / agg$wall_hi
agg$sp_hi  <- base / agg$wall_lo

ord  <- intersect(c("static", "dynamic", "segment"), unique(agg$strategy))
agg  <- agg[order(match(agg$strategy, ord), agg$n_workers), ]
pal  <- stats::setNames(c("#1B5E9C", "#C0392B", "#5C6B73")[seq_along(ord)], ord)
pchs <- stats::setNames(c(16, 17, 15)[seq_along(ord)], ord)

draw_range <- function(x, lo, hi, col) {
  v <- is.finite(lo) & is.finite(hi) & hi / pmax(lo, 1e-12) > 1.02
  if (any(v)) graphics::arrows(x[v], lo[v], x[v], hi[v], code = 3, angle = 90,
                               length = 0.03, col = col)
}
series <- function(ymed, ylo, yhi) {
  for (s in ord) {
    d <- agg[agg$strategy == s, ]
    graphics::lines(d$n_workers, d[[ymed]], col = pal[s], lwd = 2)
    draw_range(d$n_workers, d[[ylo]], d[[yhi]], pal[s])
    graphics::points(d$n_workers, d[[ymed]], col = pal[s], pch = pchs[s],
                     cex = 1)
  }
}

wk <- sort(unique(agg$n_workers))

draw <- function() {
  op <- graphics::par(mfrow = c(1, 3), mar = c(4.4, 4.5, 2.2, 0.8), las = 1,
                      cex = 0.8, mgp = c(2.8, 0.7, 0))
  on.exit(graphics::par(op), add = TRUE)

  ## ---- A: speedup ------------------------------------------------------
  plot(NA, xlim = range(wk), ylim = c(0, max(max(agg$sp_hi), max(wk)) * 1.02),
       bty = "n", xaxt = "n", xlab = "Workers", ylab = "Speedup")
  graphics::axis(1, at = wk)
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey55")
  graphics::text(max(wk), max(wk), "linear", adj = c(1.1, -0.4), cex = 0.75,
                 col = "grey40")
  graphics::mtext("A", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
  series("sp_med", "sp_lo", "sp_hi")
  graphics::legend("topleft", bty = "n", legend = ord, col = pal[ord],
                   pch = pchs[ord], lwd = 2, seg.len = 1.3, cex = 0.85)

  ## ---- B: total task time ----------------------------------------------
  # Drawn on a zero-based axis so that the growth is read as a proportion of
  # the serial run rather than as a shape floating above a cropped baseline.
  plot(NA, xlim = range(wk), ylim = c(0, max(agg$work_hi) * 1.05), bty = "n",
       xaxt = "n", xlab = "Workers", ylab = "Total task time (s)")
  graphics::axis(1, at = wk)
  graphics::abline(h = stats::median(agg$work_med[agg$n_workers == 1L]),
                   lty = 2, col = "grey55")
  graphics::text(max(wk), stats::median(agg$work_med[agg$n_workers == 1L]),
                 "serial", adj = c(1.1, -0.5), cex = 0.75, col = "grey40")
  graphics::mtext("B", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
  series("work_med", "work_lo", "work_hi")

  ## ---- C: peak memory per worker ---------------------------------------
  plot(NA, xlim = range(wk), ylim = c(0, max(agg$rss_hi) * 1.05), bty = "n",
       xaxt = "n", xlab = "Workers",
       ylab = "Peak resident set size per worker (GB)")
  graphics::axis(1, at = wk)
  graphics::mtext("C", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
  series("rss_med", "rss_lo", "rss_hi")
}

## -- Write -----------------------------------------------------------------
stem <- file.path(outdir, "figure5_parallel_strategy")

grDevices::pdf(paste0(stem, ".pdf"), width = fig_w, height = fig_h,
               useDingbats = FALSE)
draw(); grDevices::dev.off()

# LZW is lossless. The panels are line work with narrow range bars, which is
# the first detail JPEG compression inside a TIFF would soften.
grDevices::tiff(paste0(stem, ".tif"), width = fig_w, height = fig_h,
                units = "in", res = dpi, compression = "lzw",
                type = if (capabilities("cairo")) "cairo" else "quartz")
draw(); grDevices::dev.off()

best <- agg[which.max(agg$sp_med), ]
cat(sprintf("\nPeak: %s at %d workers, %.0f s, speedup %.2f, %.2f GB per worker\n",
            best$strategy, best$n_workers, best$wall_med, best$sp_med,
            best$rss_med))
cat(sprintf("Work inflation from 1 to %d workers: %.0f s to %.0f s (%.2fx)\n",
            max(wk), stats::median(agg$work_med[agg$n_workers == 1L]),
            stats::median(agg$work_med[agg$n_workers == max(wk)]),
            stats::median(agg$work_med[agg$n_workers == max(wk)]) /
              stats::median(agg$work_med[agg$n_workers == 1L])))

info <- file.info(paste0(stem, c(".pdf", ".tif")))
cat("\nWritten:\n")
for (i in seq_len(nrow(info)))
  cat(sprintf("  %-44s %6.1f MB\n", rownames(info)[i], info$size[i] / 1e6))
cat(sprintf("\n%.0f x %.0f pixels at %g dpi (%.1f x %.1f in)\n",
            fig_w * dpi, fig_h * dpi, dpi, fig_w, fig_h))
