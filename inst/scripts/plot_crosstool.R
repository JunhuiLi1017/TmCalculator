#!/usr/bin/env Rscript
# ===========================================================================
# plot_crosstool.R -- figure and table for the cross-tool comparison
#
#   Rscript inst/scripts/plot_crosstool.R --indir bench200
#
# Reads bench_crosstool.csv written by bench_crosstool.R and produces
#
#   (A) compute time against input size, log-log, with a slope-1 reference.
#       A tool with no fixed cost lies on that slope throughout; a tool whose
#       curve is flat at small n is paying a per-call overhead that has not
#       yet been amortised. The size at which the curves cross is the size
#       below which a benchmark reaches the opposite conclusion, which is the
#       point of showing the whole range rather than one input size.
#
#   (B) peak resident set size against input size. This is the half of the
#       comparison that does not favour a Bioconductor package: an R session
#       with Biostrings and GenomicRanges loaded starts near 0.8 GB whatever
#       the input, whereas an interpreter holding only the sequences grows
#       from almost nothing. Showing both makes the trade-off explicit
#       instead of leaving a reviewer to find it.
#
# Dispersion is the range over repetitions at each size, not a standard
# deviation: with three repetitions an SD is not a meaningful estimate, and
# repeated runs of the same configuration on a laptop have been observed to
# differ by nearly 30%. The centre is the median for the same reason.
#
# Base graphics only; no dependency is added to the package.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
indir  <- argval("--indir", "bench200")
outdir <- argval("--outdir", indir)
device <- argval("--device", "pdf")          # pdf, png or svg

csv <- file.path(indir, "bench_crosstool.csv")
if (!file.exists(csv)) stop("not found: ", csv)
S <- utils::read.csv(csv, stringsAsFactors = FALSE)
S <- S[S$ok & !is.na(S$compute_s), , drop = FALSE]
if (!nrow(S)) stop("no usable rows in ", csv)

# A benchmark file accumulates across invocations, so it can hold rows taken
# before and after a change to the package. Averaging those together produces
# a table that matches no build and a start-up bar that belongs to neither.
# Refuse rather than quietly plot the mixture.
if ("pkg_version" %in% names(S)) {
  vs <- sort(unique(stats::na.omit(S$pkg_version)))
  if (length(vs) > 1L)
    stop("bench_crosstool.csv mixes TmCalculator versions (",
         paste(vs, collapse = ", "), "). Re-run the benchmark with --fresh, ",
         "or keep one version with --version.")
  keep_v <- argval("--version", NA)
  if (!is.na(keep_v)) S <- S[S$pkg_version %in% keep_v, , drop = FALSE]
} else {
  message("NOTE: this file predates version stamping; it may mix builds.")
}

# Sizes that only some tools were measured at are legitimate (a slow tool is
# given a smaller sweep), but a size at which NOTHING was measured in the
# current sweep is a leftover from an earlier gradient and is easy to miss in
# a seven-row table.
tab <- table(S$n, S$tool)
lone <- rownames(tab)[rowSums(tab > 0) == 1L]
if (length(lone))
  message("sizes measured for a single tool only: ", paste(lone, collapse = ", "),
          " (expected if a slow tool has its own sweep; otherwise leftovers ",
          "from a previous gradient)")

## -- Summary: median and range over repetitions ----------------------------
summarise <- function(v) c(med = stats::median(v), lo = min(v), hi = max(v))

agg <- do.call(rbind, lapply(split(S, list(S$tool, S$n), drop = TRUE), function(d) {
  ct <- summarise(d$compute_s)
  rs <- summarise(d$rss_gb)
  st <- summarise(d$startup_s)
  # File I/O exists only because each tool is run as a separate process so
  # that peak memory can be attributed to it. It is measured and reported,
  # never added to start-up, and never drawn: a user passes data in memory.
  io <- if ("io_s" %in% names(d)) summarise(d$io_s) else c(med = NA_real_)
  data.frame(
    tool        = d$tool[1],
    n           = d$n[1],
    reps        = nrow(d),
    compute_med = ct[["med"]], compute_lo = ct[["lo"]], compute_hi = ct[["hi"]],
    # Throughput is derived from the median time rather than averaged over
    # per-run rates: the ratio of medians is what the time axis shows.
    seqs_per_s  = d$n[1] / ct[["med"]],
    rss_med     = rs[["med"]], rss_lo = rs[["lo"]], rss_hi = rs[["hi"]],
    startup_med = st[["med"]], startup_lo = st[["lo"]], startup_hi = st[["hi"]],
    io_med      = io[["med"]],
    stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$tool, agg$n), ]

## -- Fixed cost and marginal throughput ------------------------------------
# Two points at the top of the range separate the constant part of a call
# from the part that scales, which is what makes the small-n numbers
# interpretable rather than merely wrong.
fit <- do.call(rbind, lapply(split(agg, agg$tool), function(d) {
  d <- d[order(d$n), ]
  if (nrow(d) < 2L)
    return(data.frame(tool = d$tool[1], fixed_s = NA_real_,
                      marginal_seqs_per_s = NA_real_, n_points = nrow(d),
                      stringsAsFactors = FALSE))
  k <- nrow(d)
  slope <- (d$compute_med[k] - d$compute_med[k - 1L]) / (d$n[k] - d$n[k - 1L])
  data.frame(tool = d$tool[1],
             fixed_s = d$compute_med[k] - slope * d$n[k],
             marginal_seqs_per_s = 1 / slope,
             n_points = k, stringsAsFactors = FALSE)
}))

cat("\n== Median compute time and range over repetitions ==\n")
print(format(agg[, c("tool", "n", "reps", "compute_med", "compute_lo",
                     "compute_hi", "seqs_per_s", "rss_med", "startup_med",
                     "io_med")],
             digits = 4), row.names = FALSE)

# Start-up must not depend on the input. If it does, something that scales
# with n is being counted as start-up, which is what happened when file I/O
# was folded into it.
for (t in unique(agg$tool)) {
  d <- agg[agg$tool == t, ]
  if (nrow(d) > 1L && max(d$startup_med) / min(d$startup_med) > 1.5)
    message("WARNING: ", t, " start-up varies ",
            sprintf("%.1fx", max(d$startup_med) / min(d$startup_med)),
            " across input sizes; it should be constant.")
}
cat("\n== Fixed cost and marginal throughput (two largest sizes) ==\n")
print(format(fit, digits = 4), row.names = FALSE)
utils::write.csv(agg, file.path(outdir, "crosstool_summary.csv"), row.names = FALSE)

## -- Figure ----------------------------------------------------------------
tools  <- sort(unique(agg$tool))
pal    <- c("#1B5E9C", "#C0392B", "#5C6B73", "#7D3C98")[seq_along(tools)]
names(pal) <- tools
pch    <- c(16, 17, 15, 18)[seq_along(tools)]
names(pch) <- tools

open_dev <- function(file, w, h) {
  switch(device,
         pdf = grDevices::pdf(paste0(file, ".pdf"), width = w, height = h,
                              useDingbats = FALSE),
         png = grDevices::png(paste0(file, ".png"), width = w, height = h,
                              units = "in", res = 300),
         svg = grDevices::svg(paste0(file, ".svg"), width = w, height = h),
         stop("unknown --device: ", device))
}

# Error bars are drawn only where the range is wide enough to be visible;
# a bar shorter than the plotting symbol is noise dressed as information.
draw_range <- function(x, lo, hi, col) {
  vis <- is.finite(lo) & is.finite(hi) & (hi / pmax(lo, .Machine$double.eps) > 1.02)
  if (any(vis))
    graphics::arrows(x[vis], lo[vis], x[vis], hi[vis], code = 3,
                     angle = 90, length = 0.03, col = col, lwd = 1)
}

open_dev(file.path(outdir, "crosstool"), 9, 4.2)
op <- graphics::par(mfrow = c(1, 2), mar = c(4.3, 4.4, 2.2, 0.8),
                    las = 1, cex = 0.85)

## Tools shown on the time axes. A tool whose cost is orders of magnitude
## larger forces every other curve onto the axis, so it is reported in the
## table and in the text rather than being drawn; the memory panel keeps it,
## where the spread is small enough to show.
excl   <- trimws(strsplit(argval("--exclude-time", "rmelting"), ",")[[1]])
excl   <- excl[nzchar(excl)]
tt     <- setdiff(tools, excl)
aggt   <- agg[agg$tool %in% tt, ]

## Panel A -- compute time
# Linear on both axes. Both tools are essentially proportional to the input
# above a few thousand sequences, so on linear axes each is a straight line
# and the ratio of their slopes is exactly the ratio of their throughputs --
# the quantity the comparison is about, read directly off the picture. On
# logarithmic axes that ratio would look the same as every other ratio in the
# plot, including the reversed one at the smallest input.
plot(NA, xlim = c(0, max(aggt$n)), ylim = c(0, max(aggt$compute_hi) * 1.05),
     bty = "n", xlab = "Sequences", ylab = "Compute time (s)", xaxt = "n")
graphics::axis(1, at = pretty(c(0, max(aggt$n))),
               labels = format(pretty(c(0, max(aggt$n))), big.mark = ",",
                               scientific = FALSE, trim = TRUE))
graphics::mtext("A", side = 3, adj = 0, font = 2, line = 0.8, cex = 1.1)
for (t in tt) {
  d <- aggt[aggt$tool == t, ]
  d <- d[order(d$n), ]
  graphics::lines(d$n, d$compute_med, col = pal[t], lwd = 2)
  draw_range(d$n, d$compute_lo, d$compute_hi, pal[t])
  graphics::points(d$n, d$compute_med, col = pal[t], pch = pch[t], cex = 1.05)
}
graphics::legend("topleft", bty = "n", legend = tt, col = pal[tt],
                 pch = pch[tt], lwd = 2, seg.len = 1.4, cex = 0.9)
if (length(excl))
  graphics::mtext(paste(excl, collapse = ", "), side = 3, adj = 1,
                  line = 0.4, cex = 0.72, col = "grey45")

## Panel B -- peak memory
xr <- range(agg$n)
# Expand to whole decades. Taking range() directly puts the extreme points on
# the frame, which crops the symbols and their range bars and leaves the two
# curves looking closer to the edges than they are.
yr2 <- range(c(agg$rss_lo, agg$rss_hi))
yr2 <- c(10^floor(log10(yr2[1])), 10^ceiling(log10(yr2[2])))
plot(NA, xlim = xr, ylim = yr2, log = "xy", bty = "n", yaxt = "n",
     xlab = "Sequences", ylab = "Peak resident set size (GB)")
local({
  ticks <- 10^seq(log10(yr2[1]), log10(yr2[2]))
  graphics::axis(2, at = ticks,
                 labels = format(ticks, scientific = FALSE, drop0trailing = TRUE))
  graphics::axis(2, at = as.numeric(outer(2:9, ticks)), labels = FALSE,
                 tcl = -0.2)
})
graphics::mtext("B", side = 3, adj = 0, font = 2, line = 0.8, cex = 1.1)
for (t in tools) {
  d <- agg[agg$tool == t, ]
  d <- d[order(d$n), ]
  if (nrow(d) > 1L) graphics::lines(d$n, d$rss_med, col = pal[t], lwd = 2)
  draw_range(d$n, d$rss_lo, d$rss_hi, pal[t])
  graphics::points(d$n, d$rss_med, col = pal[t], pch = pch[t], cex = 1.1)
}
graphics::par(op)
grDevices::dev.off()

## -- Second figure: what one invocation actually costs -----------------------
# Panel A answers "how fast is the calculation", which is the right question
# for a genome-scale run, where a process is started once and then works for
# minutes. It is the wrong question for someone invoking a tool from the
# command line on one input, because it excludes start-up: loading R with the
# Bioconductor stack costs seconds whatever the input, and that is larger than
# an entire Biopython run at the smaller sizes.
#
# Quoting only the ratio of compute times would therefore overstate the
# practical advantage, so here start-up and compute are stacked and the size
# at which the TOTAL crosses over is marked. The two numbers differ by more
# than a factor of two, and a reader who takes one for the other will be
# wrong about which tool to use.
# Only tools measured across the whole range appear here. A tool run at two of
# the seven sizes cannot be drawn as a line without inviting the reader to
# interpolate through sizes at which it was never measured, and its cost is
# already reported, with its reason, in the table and in panel A.
keep <- agg$tool[agg$n == max(agg$n)]
drop <- setdiff(unique(agg$tool), keep)
if (length(drop))
  message("wall-time figure omits (not measured at the largest size): ",
          paste(drop, collapse = ", "))

wall <- agg[agg$tool %in% keep, c("tool", "n", "startup_med",
                                  "compute_med", "compute_lo", "compute_hi")]
wall <- wall[order(wall$n, wall$tool), ]
wall$total_med <- wall$startup_med + wall$compute_med
lab  <- paste(wall$tool, wall$n)

m <- rbind(compute = wall$compute_med, `start-up` = wall$startup_med)
colnames(m) <- lab

# One gap between input sizes, none within.
grp   <- cumsum(c(1, diff(wall$n) != 0))
space <- ifelse(c(TRUE, diff(grp) != 0), 1.1, 0.18)

# One tall bar would otherwise flatten every other bar into the axis. Rather
# than switch to a logarithmic scale -- which would make the differences that
# matter, all of them small ratios, visually indistinguishable -- the axis is
# broken: the lower segment keeps its full resolution and the outlier is shown
# above the break, compressed and unmistakably discontinuous.
brk_lo <- suppressWarnings(as.numeric(argval("--break-lo", NA)))
brk_hi <- suppressWarnings(as.numeric(argval("--break-hi", NA)))
broken <- is.finite(brk_lo) && is.finite(brk_hi) && brk_hi > brk_lo &&
          max(colSums(m)) > brk_hi

fill <- c("#1B5E9C", "#BFD3E6")
draw_bars <- function(ylim, xlabels) {
  bp <- graphics::barplot(m, space = space, col = fill, border = NA,
                          names.arg = rep("", ncol(m)), ylab = "",
                          ylim = ylim, xpd = FALSE)
  if (xlabels) {
    u <- graphics::par("usr")
    graphics::text(bp, u[3] - 0.02 * (u[4] - u[3]), labels = wall$tool,
                   srt = 40, adj = 1, xpd = TRUE, cex = 0.72)
    graphics::mtext(format(unique(wall$n), big.mark = ","), side = 1,
                    line = 3.6, at = tapply(bp, grp, mean), cex = 0.85)
    graphics::mtext("Sequences", side = 1, line = 4.7, cex = 0.85)
  }
  invisible(bp)
}
# Diagonal ticks on the axis at the break, so the discontinuity is visible on
# the axis itself and not only in the spacing of the labels.
break_marks <- function(side_up) {
  u <- graphics::par("usr"); h <- 0.018 * (u[4] - u[3])
  y <- if (side_up) u[3] else u[4]
  graphics::segments(u[1] - 0.012 * (u[2] - u[1]), c(y - h, y + h),
                     u[1] + 0.012 * (u[2] - u[1]), c(y + h, y + 3 * h),
                     xpd = TRUE, lwd = 1.2)
}

open_dev(file.path(outdir, "crosstool_walltime"), 7.5, 4.8)
if (broken) {
  graphics::layout(matrix(1:2, ncol = 1), heights = c(1, 2.6))
  op <- graphics::par(mar = c(0.2, 4.4, 2.4, 0.8), las = 1, cex = 0.85)
  draw_bars(c(brk_hi, max(colSums(m)) * 1.06), xlabels = FALSE)
  break_marks(TRUE)
  graphics::par(mar = c(5.6, 4.4, 0.4, 0.8))
  bp <- draw_bars(c(0, brk_lo), xlabels = TRUE)
  break_marks(FALSE)
  graphics::mtext("Elapsed time (s)", side = 2, line = 3, las = 0,
                  cex = 0.85, adj = 0.2)
} else {
  op <- graphics::par(mar = c(5.6, 4.4, 2.4, 0.8), las = 1, cex = 0.85)
  bp <- draw_bars(c(0, max(colSums(m)) * 1.18), xlabels = TRUE)
  graphics::mtext("Elapsed time (s)", side = 2, line = 3, las = 0, cex = 0.85)
}
graphics::legend("topleft", bty = "n", fill = fill, border = NA,
                 legend = c("compute", "start-up"), cex = 0.9)

# Where the TOTAL, not the rate, crosses over.
if (nrow(fit) >= 2L && all(is.finite(fit$fixed_s))) {
  o <- order(fit$marginal_seqs_per_s, decreasing = TRUE)
  a <- fit[o[1], ]; b <- fit[o[2], ]
  su <- vapply(list(a, b), function(z)
    stats::median(agg$startup_med[agg$tool == z$tool]), numeric(1))
  n_eq <- ((a$fixed_s + su[1]) - (b$fixed_s + su[2])) /
          (1 / b$marginal_seqs_per_s - 1 / a$marginal_seqs_per_s)
  if (is.finite(n_eq) && n_eq > 0)
    graphics::mtext(sprintf(
      "including start-up, %s overtakes %s at about %s sequences",
      a$tool, b$tool, format(signif(n_eq, 2), big.mark = ",")),
      side = 3, line = 0.4, adj = 1, cex = 0.78, col = "grey35")
}
graphics::par(op)
grDevices::dev.off()

cat("\nWritten ", file.path(outdir, paste0("crosstool.", device)), ", ",
    file.path(outdir, paste0("crosstool_walltime.", device)), " and ",
    file.path(outdir, "crosstool_summary.csv"), "\n", sep = "")
