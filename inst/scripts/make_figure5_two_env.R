#!/usr/bin/env Rscript
# ===========================================================================
# make_figure5_two_env.R -- Figure 5, two environments on one pair of axes
#
#   (A) wall-clock time against worker count, start-up included
#   (B) peak resident set size of the heaviest worker
#
#   Rscript inst/scripts/make_figure5_two_env.R \
#       --csv-hpc results/bench_parallel_cluster.csv \
#       --csv-mac bench_parallel_strategy.csv
#
# Replaces the three-panel single-environment Figure 5 and the two 18-row
# sweep tables. The node is overlaid on the laptop so the one comparison
# that matters is the first thing the eye does: on the laptop two of the
# three curves in A turn upward beyond their optimum, on the node none do.
# B gives the reason and the remedy. The total-task-time panel of the old
# figure is not reproduced here; those numbers are quoted in the text and
# tabulated in full in the vignette.
#
# Colour is the strategy and nothing else; line type and marker fill are the
# environment and nothing else. Solid and filled is the laptop, dashed and
# open is the node.
#
# The time metric defaults to wall, start-up included, which is what the
# revised Table 6 reports (its note 1 says so). See --metric below.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}

csv_hpc <- argval("--csv-hpc", "")
csv_mac <- argval("--csv-mac", "")
outdir  <- argval("--outdir", "figures")
dpi     <- as.numeric(argval("--dpi", "600"))
fmt     <- argval("--format", "tif")
fig_w   <- as.numeric(argval("--width",  "9.0"))    # MDPI full width, 2 panels
fig_h   <- as.numeric(argval("--height", "4.1"))
lab_hpc <- argval("--label-hpc", "Compute node")
lab_mac <- argval("--label-mac", "Laptop")

# wall | compute. Default wall, start-up included: it is the time a user
# waits, and it is what the revised Table 6 reports. `compute` subtracts the
# ~8-10 s of SnowParam start-up and isolates the scaling of the computation
# itself; the previous Table 6 used that convention.
#
# The two must never be mixed between environments. An earlier version of
# this script read wall_compute_s from the cluster CSV while the embedded
# laptop series held wall_s, which shifted every laptop point upward by its
# own start-up and silently corrupted the comparison.
metric <- match.arg(argval("--metric", "wall"), c("wall", "compute"))

if (!nzchar(csv_hpc) || !file.exists(csv_hpc))
  stop("pass the cluster sweep with --csv-hpc ",
       "(results/bench_parallel_cluster.csv)")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## -- Load ------------------------------------------------------------------
# Both time columns are required from every source, so that the choice of
# metric is made here rather than being fixed by whichever column a given
# CSV happens to carry.
need <- c("strategy", "n_workers", "wall_s", "wall_compute_s", "max_rss_gb",
          "n_windows")
# Braces, not a bare multi-line if/else: at top level R closes the `if` at
# the newline and then fails on a stray `else`.
tcol <- if (identical(metric, "wall")) "wall_s" else "wall_compute_s"
tlab <- if (identical(metric, "wall")) {
  "Wall-clock time (s)"
} else {
  "Wall time, start-up excluded (s)"
}

read_env <- function(path, label) {
  S <- utils::read.csv(path, stringsAsFactors = FALSE)
  miss <- setdiff(need, names(S))
  if (length(miss))
    stop(label, ": missing column(s) ", paste(miss, collapse = ", "),
         ". This summary predates the split of start-up from compute; re-run.")
  # Configurations covering different windows are not comparable, and the
  # difference is a boundary mistake rather than a measurement.
  if (length(unique(S$n_windows)) > 1L)
    stop(label, ": configurations cover different numbers of windows (",
         paste(unique(S$n_windows), collapse = ", "), "); timings are not ",
         "comparable.")
  S$env <- label
  S
}

H <- read_env(csv_hpc, lab_hpc)

if (!nzchar(csv_mac) || !file.exists(csv_mac))
  stop("pass the laptop sweep with --csv-mac (bench_parallel_strategy.csv). ",
       "An embedded fallback used to live here; it held compute times under ",
       "the name wall_s and was removed once the real CSV was available.")
M <- read_env(csv_mac, lab_mac)

# The two sweeps must describe the same computation, or the panels compare
# different problems and the divergence means nothing.
if (!identical(unique(H$n_windows), unique(M$n_windows)))
  stop("the two environments cover different numbers of windows (",
       unique(H$n_windows), " vs ", unique(M$n_windows),
       "); they are not comparable.")

S <- rbind(H[, c(need, "env")], M[, c(need, "env")])

## -- Median and range over repetitions -------------------------------------
summ <- function(v) c(med = stats::median(v), lo = min(v), hi = max(v))
agg <- do.call(rbind, lapply(
  split(S, list(S$env, S$strategy, S$n_workers), drop = TRUE), function(d) {
    w <- summ(d[[tcol]]); m <- summ(d$max_rss_gb)
    data.frame(env = d$env[1], strategy = d$strategy[1],
               n_workers = d$n_workers[1], n_rep = nrow(d),
               wall_med = w[["med"]], wall_lo = w[["lo"]], wall_hi = w[["hi"]],
               rss_med  = m[["med"]], rss_lo  = m[["lo"]], rss_hi  = m[["hi"]],
               stringsAsFactors = FALSE)
  }))

# Speedup against each strategy's own one-worker run in the same environment,
# as Table 6 note 3 defines it. The range bars invert: a longer time is a
# smaller speedup.
ser <- agg[agg$n_workers == 1L, c("env", "strategy", "wall_med")]
if (!nrow(ser))
  stop("no one-worker run in the sweep, so there is no serial baseline.")
key <- paste(agg$env, agg$strategy)
base <- stats::setNames(ser$wall_med, paste(ser$env, ser$strategy))[key]
agg$sp_med <- base / agg$wall_med
agg$sp_lo  <- base / agg$wall_hi
agg$sp_hi  <- base / agg$wall_lo

ord  <- intersect(c("static", "dynamic", "segment"), unique(agg$strategy))
envs <- c(lab_hpc, lab_mac)
pal  <- stats::setNames(c("#1B5E9C", "#C0392B", "#5C6B73")[seq_along(ord)], ord)
# One visual channel per variable. Colour is the strategy and nothing else;
# line type and marker fill are the environment and nothing else. Strategy
# was briefly encoded twice, in colour and in marker shape, which costs the
# reader a second key to memorise and buys nothing once the colours are
# distinct. Markers are kept only to show where the six measured worker
# counts fall, so a single shape serves, filled for the laptop and open for
# the node.
#
# The laptop carries the unmarked style, solid and filled, because it is the
# environment the reader met in the earlier sections; the node then reads as
# the comparison drawn against it.
# The dash is given as an on/off hex pattern rather than lty = 2, which is
# "44": four units on, four off, in multiples of lwd. At lwd 2 that is a
# 16-unit cycle, long enough that a short legend key or a steep segment of a
# curve can fall entirely inside one dash and read as solid. "22" halves it.
dash <- argval("--dash", "22")
ltys <- stats::setNames(c(dash, "solid"), envs) # envs = c(node, laptop)
pch_for <- function(e) if (identical(e, lab_mac)) 16 else 1
agg  <- agg[order(match(agg$env, envs), match(agg$strategy, ord),
                  agg$n_workers), ]
wk   <- sort(unique(agg$n_workers))

draw_range <- function(x, lo, hi, col) {
  v <- is.finite(lo) & is.finite(hi) & hi / pmax(lo, 1e-12) > 1.02
  if (any(v)) graphics::arrows(x[v], lo[v], x[v], hi[v], code = 3, angle = 90,
                               length = 0.025, col = col)
}
series <- function(ymed, ylo, yhi) {
  for (e in envs) for (s in ord) {
    d <- agg[agg$env == e & agg$strategy == s, ]
    if (!nrow(d)) next
    graphics::lines(d$n_workers, d[[ymed]], col = pal[s], lwd = 2,
                    lty = ltys[e])
    draw_range(d$n_workers, d[[ylo]], d[[yhi]], pal[s])
    graphics::points(d$n_workers, d[[ymed]], col = pal[s],
                     pch = pch_for(e), cex = 0.95, lwd = 1.6)
  }
}

# One legend, two columns: environments on the left as line keys, strategies
# on the right as colour keys. The two things the reader has to decode are
# then side by side rather than stacked in separate boxes.
#
# legend() fills column-major, so both columns must hold the same number of
# rows or the shorter one spills into the longer. The short column is padded
# with a blank label whose lty, pch and col are all NA, which draws nothing.
add_legends <- function(where = "topright") {
  # Listed solid first, which is the laptop; `envs` keeps the node first
  # because that is the order the companion table wants.
  el <- c(lab_mac, lab_hpc)
  nrow <- max(length(el), length(ord))
  pad  <- function(x, fill = NA) c(x, rep(fill, nrow - length(x)))
  # seg.len must hold at least two whole dash-gap cycles, or the dashed key
  # renders as one short stroke and reads as solid. With the "22" pattern the
  # cycle is half what lty = 2 gives, so the key no longer has to be stretched
  # as far to show it.
  graphics::legend(
    where, bty = "n", cex = 0.78, ncol = 2, seg.len = 2.6,
    x.intersp = 0.9, y.intersp = 1.15,
    legend = c(pad(el, ""),                      pad(ord, "")),
    lty    = c(pad(unname(ltys[el])),            pad(rep("solid", length(ord)))),
    # Strategy keys are colour only; the environment keys carry the markers,
    # since that is the variable the marker fill encodes.
    pch    = c(pad(vapply(el, pch_for, numeric(1))),
               pad(rep(NA, length(ord)))),
    col    = c(pad(rep("grey25", length(el))),   pad(unname(pal[ord]))),
    lwd    = 2)
}

draw <- function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.6, 2.2, 0.8), las = 1,
                      cex = 0.8, mgp = c(3.0, 0.7, 0))
  on.exit(graphics::par(op), add = TRUE)

  ## ---- A: wall clock ---------------------------------------------------
  plot(NA, xlim = range(wk), ylim = c(0, max(agg$wall_hi) * 1.04), bty = "n",
       xaxt = "n", xlab = "Workers", ylab = tlab)
  graphics::axis(1, at = wk)
  graphics::mtext("A", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
  series("wall_med", "wall_lo", "wall_hi")
  add_legends()

  ## ---- B: peak memory per worker ---------------------------------------
  plot(NA, xlim = range(wk), ylim = c(0, max(agg$rss_hi) * 1.05), bty = "n",
       xaxt = "n", xlab = "Workers",
       ylab = "Peak resident set size per worker (GB)")
  graphics::axis(1, at = wk)
  graphics::mtext("B", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
  series("rss_med", "rss_lo", "rss_hi")
}

## -- Export ----------------------------------------------------------------
stem <- file.path(outdir, paste0("figure5_two_env_", metric))
if (identical(fmt, "tif")) {
  # On macOS tiff() defaults to type = "quartz", which silently ignores
  # `compression`. At 600 dpi the uncompressed raster is tens of megabytes,
  # which journals reject, so the cairo device is requested when the build
  # has it and LZW is only asked for when it will actually be honoured.
  args <- list(filename = paste0(stem, ".tif"), width = fig_w,
               height = fig_h, units = "in", res = dpi)
  # capabilities() lives in base, not grDevices.
  if (isTRUE(unname(capabilities("cairo")))) {
    args$type <- "cairo"
    args$compression <- "lzw"
  } else {
    warning("no cairo device; writing an uncompressed TIFF. Convert with ",
            "`tiffcp -c lzw in.tif out.tif` before submission.",
            call. = FALSE)
  }
  do.call(grDevices::tiff, args)
  draw(); grDevices::dev.off()
  message("wrote ", stem, ".tif  (",
          round(file.info(paste0(stem, ".tif"))$size / 1e6, 1), " MB)")
} else if (identical(fmt, "jpg")) {
  grDevices::jpeg(paste0(stem, ".jpg"), width = fig_w, height = fig_h,
                  units = "in", res = dpi, quality = 95)
  draw(); grDevices::dev.off()
  message("wrote ", stem, ".jpg")
} else {
  draw()
}

## -- The six-row companion table -------------------------------------------
# One row per environment and strategy: the best configuration and where it
# occurs. This is what goes in the results section in place of the full
# thirty-six-row sweep, which belongs in the vignette.
best <- do.call(rbind, lapply(split(agg, list(agg$env, agg$strategy),
                                    drop = TRUE), function(d) {
  # Every speedup below is divided by this strategy's own one-worker run, so
  # its absence is an error rather than something to paper over with the
  # fastest configuration available.
  ser <- d$wall_med[d$n_workers == 1L]
  if (!length(ser))
    stop(d$env[1], " / ", d$strategy[1],
         ": no one-worker run, so there is no serial baseline. ",
         "Re-run the sweep with --workers 1,...")
  i <- which.min(d$wall_med)
  data.frame(Environment = d$env[1], Strategy = d$strategy[1],
             `Serial (s)` = round(ser, 1),
             `Best (s)` = sprintf("%.1f [%.1f-%.1f]", d$wall_med[i],
                                  d$wall_lo[i], d$wall_hi[i]),
             Workers = d$n_workers[i],
             Speedup = round(ser / d$wall_med[i], 2),
             `Peak RSS per worker (GB)` = round(d$rss_med[i], 2),
             Reps = d$n_rep[i], Metric = metric, check.names = FALSE,
             stringsAsFactors = FALSE)
}))
best <- best[order(match(best$Environment, envs),
                   match(best$Strategy, ord)), ]
rownames(best) <- NULL
tab_path <- file.path(outdir, paste0("table6_best_config_", metric, ".csv"))
utils::write.csv(best, tab_path, row.names = FALSE)
message("wrote ", tab_path)
print(best, row.names = FALSE)
