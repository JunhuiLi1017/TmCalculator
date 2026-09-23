#!/usr/bin/env Rscript
# ===========================================================================
# make_figure5.R -- Figure 5: the hg38 worker-count sweep in two environments
#
#   (A) wall-clock time against worker count, worker start-up included
#   (B) peak memory of the whole job, summed over manager and workers
#
#   Rscript inst/scripts/make_figure5.R                       # shipped CSVs
#   Rscript inst/scripts/make_figure5.R --outdir figures \
#       --laptop results_laptop/bench_tm_calculate_summary.csv \
#       --cluster results_16g/bench_tm_calculate_summary.csv
#
# WHAT REPLACED WHAT. The earlier Figure 5 compared three task-partitioning
# strategies, because at that point the question was which one tm_calculate()
# should use. That was settled: segments, longest first. The strategies are
# gone from the package, so a figure that still showed three curves would be
# documenting a choice the software no longer offers. This one asks the
# question a user actually has, which is how many workers to give it, and
# answers it on two machines.
#
# ONE VISUAL CHANNEL PER VARIABLE. Line type and marker fill are the
# environment and nothing else: solid and filled is the laptop, dashed and
# open is the compute node. Everything is grey, so the figure survives
# greyscale printing and colour-blind readers without a second encoding.
#
# WHY THE WHOLE JOB RATHER THAN THE HEAVIEST WORKER IN PANEL B. Both are
# recorded. The per-worker figure is the noisier of the two, because the
# sampler wakes every ten seconds and catches whichever moment it catches:
# across three repetitions at six workers it ranged 1.60 to 2.81 GB on the
# laptop for identical work. The summed figure is steadier (8.53 to 9.09 GB
# over the same three runs) and is the one that decides whether a machine can
# host the run at all; it is also what a cluster scheduler enforces. Pass
# --panel-b worker to plot the per-worker series instead.
#
# THE TWO SERIES IN PANEL B ARE NOT COMPARABLE WITH EACH OTHER. The laptop is
# sampled with ps and the node with /proc/*/smaps_rollup, and macOS compresses
# memory, so a difference between the curves may be the tools rather than the
# software. Each curve against its own machine's limit is the reading that
# holds; that is why the reference line is drawn and the two are not compared
# in the caption.
#
# WHY WALL CLOCK RATHER THAN COMPUTE TIME. Worker start-up is about nine
# seconds per call and a user waits through it, so it belongs in the number
# the figure reports. The sweep also records compute_s with start-up removed;
# pass --metric compute to plot that instead, which is the right choice only
# when comparing against a machine whose start-up differs.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
outdir  <- argval("--outdir", ".")
metric  <- match.arg(argval("--metric", "wall"), c("wall", "compute"))
fmt     <- match.arg(argval("--format", "tif"), c("tif", "pdf", "png"))
dpi     <- as.numeric(argval("--dpi", "600"))
fig_w   <- as.numeric(argval("--width", "7.2"))
fig_h   <- as.numeric(argval("--height", "3.4"))
panel_b <- match.arg(argval("--panel-b", "job"), c("job", "worker"))
mem_gb  <- as.numeric(argval("--memory-limit", "16"))   # the reference line

shipped <- function(f) {
  p <- system.file("extdata", f, package = "TmCalculator")
  if (nzchar(p)) return(p)
  q <- file.path("inst", "extdata", f)          # running from a checkout
  if (file.exists(q)) return(q)
  stop("cannot find ", f, ": install the package or run from its root")
}
lap_csv <- argval("--laptop",  NA_character_)
clu_csv <- argval("--cluster", NA_character_)
if (is.na(lap_csv)) lap_csv <- shipped("bench_hg38_laptop.csv")
if (is.na(clu_csv)) clu_csv <- shipped("bench_hg38_cluster.csv")

# The columns the figure needs. Checked up front rather than discovered
# halfway through, so a summary from an older version of the benchmark fails
# with a message that names the missing column.
NEED <- c("n_workers", "wall_s", "lo", "hi", "compute_s",
          if (panel_b == "job") "peak_job_gb" else "peak_worker_gb")
read_env <- function(path, env) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  miss <- setdiff(NEED, names(d))
  if (length(miss))
    stop(path, " has no column(s): ", paste(miss, collapse = ", "))
  d$env <- env
  d[order(d$n_workers), ]
}
d <- rbind(read_env(lap_csv,  "Laptop"),
           read_env(clu_csv, "Compute node"))
envs <- c("Laptop", "Compute node")
message("laptop  : ", lap_csv, "  (", sum(d$env == "Laptop"), " worker counts)")
message("cluster : ", clu_csv, "  (", sum(d$env != "Laptop"), " worker counts)")

# compute_s has no measured range: start-up is calibrated once per worker
# count, so lo and hi would be the wall-clock range shifted by a constant.
# Drawing them as though they were measured would overstate what is known.
d$y  <- if (metric == "wall") d$wall_s else d$compute_s
d$lo <- if (metric == "wall") d$lo else NA_real_
d$hi <- if (metric == "wall") d$hi else NA_real_

d$mem <- if (panel_b == "job") d$peak_job_gb else d$peak_worker_gb
ylab_time <- if (metric == "wall") "Wall clock (s)" else "Compute time (s)"
# Parenthesised: at top level R closes the `if` at the end of the line and the
# bare `else` on the next one is a syntax error.
ylab_mem  <- if (panel_b == "job") "Peak memory, whole job (GB)" else
             "Peak resident size per worker (GB)"
GREY <- "grey20"; BAR <- "grey45"

draw <- function() {
  op <- par(mfrow = c(1, 2), mar = c(4.0, 4.3, 2.0, 0.8),
            mgp = c(2.5, 0.7, 0), cex.axis = 0.9, cex.lab = 0.95)
  on.exit(par(op), add = TRUE)

  for (what in c("time", "rss")) {
    y_all <- if (what == "time") c(d$y, d$lo, d$hi) else c(d$mem, mem_gb)
    y_all <- y_all[is.finite(y_all)]
    plot(range(d$n_workers), c(0, max(y_all) * 1.08), type = "n",
         xlab = "Workers", xaxt = "n", yaxs = "i",
         ylab = if (what == "time") ylab_time else ylab_mem)
    axis(1, at = sort(unique(d$n_workers)))
    # The limit is the point of panel B: a curve is read against it, not
    # against the other curve.
    if (what == "rss" && is.finite(mem_gb)) {
      abline(h = mem_gb, lty = 3, col = "grey55")
      text(min(d$n_workers), mem_gb, sprintf("%g GB available", mem_gb),
           adj = c(0, -0.4), cex = 0.75, col = "grey35")
    }
    mtext(if (what == "time") "A" else "B", side = 3, line = 0.4,
          adj = 0, font = 2, cex = 1.1)

    for (e in envs) {
      s <- d[d$env == e, ]
      solid <- e == "Laptop"
      y <- if (what == "time") s$y else s$mem
      if (what == "time" && any(is.finite(s$lo)))
        arrows(s$n_workers, s$lo, s$n_workers, s$hi, angle = 90, code = 3,
               length = 0.03, col = BAR, lwd = 1)
      lines(s$n_workers, y, lty = if (solid) 1 else 2, col = GREY, lwd = 1.6)
      points(s$n_workers, y, pch = if (solid) 19 else 1, col = GREY, cex = 1.1)
    }
    if (what == "time")
      legend("topright", envs, lty = c(1, 2), pch = c(19, 1), lwd = 1.6,
             bty = "n", cex = 0.85, col = GREY)
  }
}

## -- Export ----------------------------------------------------------------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stem <- file.path(outdir, paste0("figure5_hg38_", metric))
if (identical(fmt, "tif")) {
  # On macOS tiff() defaults to type = "quartz", which silently ignores
  # `compression`. At 600 dpi the uncompressed raster is tens of megabytes,
  # which journals reject, so the cairo device is requested when the build
  # has it and LZW is only asked for when it will actually be honoured.
  a <- list(filename = paste0(stem, ".tif"), width = fig_w, height = fig_h,
            units = "in", res = dpi)
  if (isTRUE(unname(capabilities("cairo")))) {
    a$type <- "cairo"; a$compression <- "lzw"
  } else {
    warning("no cairo device; writing an uncompressed TIFF. Convert with ",
            "`tiffcp -c lzw in.tif out.tif` before submission.", call. = FALSE)
  }
  do.call(grDevices::tiff, a)
} else if (identical(fmt, "pdf")) {
  grDevices::pdf(paste0(stem, ".pdf"), width = fig_w, height = fig_h)
} else {
  grDevices::png(paste0(stem, ".png"), width = fig_w, height = fig_h,
                 units = "in", res = dpi)
}
draw(); grDevices::dev.off()
f <- paste0(stem, ".", fmt)
message("wrote ", f, "  (", round(file.info(f)$size / 1e6, 2), " MB)")

## -- The numbers the caption and Section 3.5 quote --------------------------
# Printed rather than left to be read off the figure, so that the text and
# the picture cannot drift apart.
for (e in envs) {
  s <- d[d$env == e, ]
  i <- which.min(s$y)
  cat(sprintf("\n%s\n", e))
  cat(sprintf("  serial            : %.1f s\n", s$y[s$n_workers == 1]))
  cat(sprintf("  fastest           : %.1f s at %d workers (%.2fx)\n",
              s$y[i], s$n_workers[i], s$y[s$n_workers == 1] / s$y[i]))
  cat(sprintf("  widest range      : %.1f s at %d workers\n",
              max(s$hi - s$lo, na.rm = TRUE),
              s$n_workers[which.max(s$hi - s$lo)]))
  cat(sprintf("  peak memory       : %.2f GB at %d workers (%.0f%% of %g GB)\n",
              max(s$mem), s$n_workers[which.max(s$mem)],
              100 * max(s$mem) / mem_gb, mem_gb))
  # Amdahl, fitted on 1/n. The constant is the part no worker count removes.
  fit <- stats::lm(s$y ~ I(1 / s$n_workers))
  cf <- stats::coef(fit)
  cat(sprintf(paste0("  Amdahl            : T(n) = %.1f + %.1f/n  ",
                     "(R2 = %.4f, serial %.1f%%, ceiling %.1fx)\n"),
              cf[1], cf[2], summary(fit)$r.squared,
              100 * cf[1] / (cf[1] + cf[2]), (cf[1] + cf[2]) / cf[1]))
}
