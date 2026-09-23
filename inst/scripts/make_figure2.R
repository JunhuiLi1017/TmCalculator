#!/usr/bin/env Rscript
# ===========================================================================
# make_figure2.R -- publication-resolution export of Figure 2
#
#   (A) compute time against input size
#   (B) elapsed time for one invocation, split into start-up and compute
#
#   Rscript inst/scripts/make_figure2.R
#   Rscript inst/scripts/make_figure2.R --csv bench200/crosstool_bench.csv
#   Rscript inst/scripts/make_figure2.R --memory        # add a third panel
#
# Panel A is linear on both axes rather than log-log. Both tools are
# proportional to their input above a few thousand sequences, so on linear
# axes each is a straight line and the ratio of the slopes is the ratio of
# the throughputs, read directly off the picture. Logarithmic axes give every
# ratio the same visual weight, which would make the fourfold separation at
# the largest input look no larger than the reversed ranking at the smallest.
#
# Panel B answers a different question from panel A and the two must not be
# conflated. Panel A is the right comparison for a genome-scale run, where a
# process starts once and then works for minutes. Panel B is the right one
# for a tool invoked from the command line on a single input, because it
# includes the cost of attaching R and its Bioconductor dependencies, which
# is larger than an entire Biopython run at the smaller sizes.
#
# File input and output is measured by the benchmark and deliberately drawn
# in neither panel: it exists only because each tool is run as a separate
# process so that peak memory can be attributed to it, and it scales with the
# input, so including it would inflate both tools by an artefact of the
# harness rather than a property of the software.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
has_flag <- function(f) f %in% args

csv     <- argval("--csv", system.file("extdata", "crosstool_bench.csv",
                                       package = "TmCalculator"))
outdir  <- argval("--outdir", "figures")
dpi     <- as.numeric(argval("--dpi", "600"))
add_mem <- has_flag("--memory")
# Panel B is categorical, so every input size costs two bars and a label.
# Drawing all eight leaves the group labels overlapping; four spanning the
# same range carry the same message legibly. Panel A keeps every size.
bar_sizes <- as.numeric(strsplit(argval("--bar-sizes",
                                        "1000,10000,100000,300000"), ",")[[1]])
fig_w   <- as.numeric(argval("--width",  if (add_mem) "10.5" else "7.5"))
fig_h   <- as.numeric(argval("--height", "3.9"))

if (!nzchar(csv) || !file.exists(csv))
  stop("benchmark results not found. Run inst/scripts/bench_crosstool.R and ",
       "pass the CSV with --csv, or install the package with the results in ",
       "inst/extdata/crosstool_bench.csv.")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

S <- utils::read.csv(csv, stringsAsFactors = FALSE)
S <- S[S$ok & !is.na(S$compute_s), , drop = FALSE]
if (!nrow(S)) stop("no usable rows in ", csv)

# Rows taken before and after a change to the package describe no single
# build; averaging them silently is how a start-up bar ends up belonging to
# neither. Refuse instead.
if ("pkg_version" %in% names(S)) {
  vs <- sort(unique(stats::na.omit(S$pkg_version)))
  if (length(vs) > 1L)
    stop("the benchmark file mixes TmCalculator versions (",
         paste(vs, collapse = ", "), "); re-run with --fresh.")
} else {
  stop("this benchmark file predates the split of start-up from file I/O; ",
       "re-run inst/scripts/bench_crosstool.R.")
}

## -- Median and range over repetitions -------------------------------------
summ <- function(v) c(med = stats::median(v), lo = min(v), hi = max(v))
agg <- do.call(rbind, lapply(split(S, list(S$tool, S$n), drop = TRUE), function(d) {
  ct <- summ(d$compute_s); rs <- summ(d$rss_gb); st <- summ(d$startup_s)
  data.frame(tool = d$tool[1], n = d$n[1],
             compute_med = ct[["med"]], compute_lo = ct[["lo"]],
             compute_hi  = ct[["hi"]],
             rss_med = rs[["med"]], rss_lo = rs[["lo"]], rss_hi = rs[["hi"]],
             startup_med = st[["med"]], stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$tool, agg$n), ]

# A tool measured over only part of the range is left off the time axes: its
# cost is orders of magnitude larger, which flattens the others onto the
# axis, and two or three points cannot be drawn as a line through sizes at
# which nothing was measured.
tt   <- sort(unique(agg$tool[agg$n == max(agg$n)]))
excl <- setdiff(sort(unique(agg$tool)), tt)
aggt <- agg[agg$tool %in% tt, ]

pal <- c("#1B5E9C", "#C0392B", "#5C6B73")[seq_along(sort(unique(agg$tool)))]
names(pal) <- sort(unique(agg$tool))
pch <- c(16, 17, 15)[seq_along(pal)]; names(pch) <- names(pal)

draw_range <- function(x, lo, hi, col) {
  v <- is.finite(lo) & is.finite(hi) & hi / pmax(lo, 1e-12) > 1.02
  if (any(v)) graphics::arrows(x[v], lo[v], x[v], hi[v], code = 3, angle = 90,
                               length = 0.03, col = col)
}
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

draw <- function() {
  op <- graphics::par(mfrow = c(1, if (add_mem) 3 else 2),
                      mar = c(4.8, 4.5, 2.2, 0.8), las = 1, cex = 0.8,
                      mgp = c(2.7, 0.7, 0))
  on.exit(graphics::par(op), add = TRUE)

  ## ---- A: compute time -------------------------------------------------
  plot(NA, xlim = c(0, max(aggt$n)), ylim = c(0, max(aggt$compute_hi) * 1.05),
       bty = "n", xlab = "Sequences", ylab = "Compute time (s)", xaxt = "n")
  graphics::axis(1, at = pretty(c(0, max(aggt$n))),
                 labels = fmt_n(pretty(c(0, max(aggt$n)))))
  graphics::mtext("A", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
  for (t in tt) {
    d <- aggt[aggt$tool == t, ]; d <- d[order(d$n), ]
    graphics::lines(d$n, d$compute_med, col = pal[t], lwd = 2)
    draw_range(d$n, d$compute_lo, d$compute_hi, pal[t])
    graphics::points(d$n, d$compute_med, col = pal[t], pch = pch[t], cex = 1)
  }
  graphics::legend("topleft", bty = "n", legend = tt, col = pal[tt],
                   pch = pch[tt], lwd = 2, seg.len = 1.3, cex = 0.85)

  ## ---- B: elapsed time for one invocation ------------------------------
  # One hue, light for start-up and dark for compute, so the legend needs two
  # entries rather than one per tool. Restricting to four input sizes leaves
  # eight bars, few enough for the tool names to be written under them without
  # colliding; all eight sizes would not fit.
  wall <- aggt[aggt$n %in% bar_sizes, ]
  wall <- wall[order(wall$n, wall$tool), ]
  m <- rbind(compute = wall$compute_med, `start-up` = wall$startup_med)
  fill <- c("#1B5E9C", "#BFD3E6")

  grp   <- cumsum(c(1, diff(wall$n) != 0))
  space <- ifelse(c(TRUE, diff(grp) != 0), 1.0, 0.2)

  bp <- graphics::barplot(m, space = space, col = fill, border = NA,
                          names.arg = rep("", ncol(m)), ylab = "",
                          ylim = c(0, max(colSums(m)) * 1.20))
  graphics::mtext("B", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
  graphics::mtext("Elapsed time (s)", side = 2, line = 3, las = 0, cex = 0.8)
  u <- graphics::par("usr")
  graphics::text(bp, u[3] - 0.02 * (u[4] - u[3]), labels = wall$tool,
                 srt = 40, adj = 1, xpd = TRUE, cex = 0.62)
  graphics::mtext(fmt_n(unique(wall$n)), side = 1, line = 2.7,
                  at = tapply(bp, grp, mean), cex = 0.75)
  graphics::mtext("Sequences", side = 1, line = 3.7, cex = 0.8)
  graphics::legend("topleft", bty = "n", fill = fill, border = NA,
                   legend = c("compute", "start-up"), cex = 0.85)

  ## ---- C: peak memory (optional) ---------------------------------------
  if (add_mem) {
    yr <- range(c(agg$rss_lo, agg$rss_hi))
    yr <- c(10^floor(log10(yr[1])), 10^ceiling(log10(yr[2])))
    plot(NA, xlim = range(agg$n), ylim = yr, log = "xy", bty = "n", yaxt = "n",
         xlab = "Sequences", ylab = "Peak resident set size (GB)")
    ticks <- 10^seq(log10(yr[1]), log10(yr[2]))
    graphics::axis(2, at = ticks,
                   labels = format(ticks, scientific = FALSE,
                                   drop0trailing = TRUE))
    graphics::axis(2, at = as.numeric(outer(2:9, ticks)), labels = FALSE,
                   tcl = -0.2)
    graphics::mtext("C", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
    for (t in names(pal)) {
      d <- agg[agg$tool == t, ]; d <- d[order(d$n), ]
      if (nrow(d) > 1L) graphics::lines(d$n, d$rss_med, col = pal[t], lwd = 2)
      draw_range(d$n, d$rss_lo, d$rss_hi, pal[t])
      graphics::points(d$n, d$rss_med, col = pal[t], pch = pch[t], cex = 1)
    }
    graphics::legend("topleft", bty = "n", legend = names(pal), col = pal,
                     pch = pch, lwd = 2, seg.len = 1.3, cex = 0.85)
  }
}

## -- Write -----------------------------------------------------------------
stem <- file.path(outdir, "figure2_tool_comparison")

grDevices::pdf(paste0(stem, ".pdf"), width = fig_w, height = fig_h,
               useDingbats = FALSE)
draw(); grDevices::dev.off()

# LZW is lossless. JPEG compression inside a TIFF would soften the thin lines
# and the range bars, which are the only fine detail the panels contain.
grDevices::tiff(paste0(stem, ".tif"), width = fig_w, height = fig_h,
                units = "in", res = dpi, compression = "lzw",
                type = if (capabilities("cairo")) "cairo" else "quartz")
draw(); grDevices::dev.off()

if (length(excl))
  message("omitted from the time panels (not measured across the full range): ",
          paste(excl, collapse = ", "))
info <- file.info(paste0(stem, c(".pdf", ".tif")))
cat("\nWritten:\n")
for (i in seq_len(nrow(info)))
  cat(sprintf("  %-42s %6.1f MB\n", rownames(info)[i], info$size[i] / 1e6))
cat(sprintf("\n%.0f x %.0f pixels at %g dpi (%.1f x %.1f in)\n",
            fig_w * dpi, fig_h * dpi, dpi, fig_w, fig_h))
