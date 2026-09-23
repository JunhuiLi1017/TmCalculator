#!/usr/bin/env Rscript
# ===========================================================================
# make_figure3.R -- publication-resolution export of Figure 3
#
#   (A) circular map of the E. coli K-12 MG1655 chromosome, all tracks
#   (B) the same map restricted to two selected regions
#   (C) a high-magnification locus containing groups of windows with
#       identical GC content but different Tm
#
#   Rscript inst/scripts/make_figure3.R --in-peaks
#   Rscript inst/scripts/make_figure3.R --n-regions 1 --n-groups 4
#   Rscript inst/scripts/make_figure3.R --region "U00096.3:212801-232800"
#
# LAYOUT. One row. Each zoom panel is drawn into its own layout cell with a
# single zoom region: plot_genome_track() calls graphics::layout() itself
# when given more than one region (plot_genome_track.R:806-809), and R's
# layout is global, so letting it do that from inside this figure's layout
# would discard the arrangement of A and B. --n-regions adds cells to the
# row; widen the figure to match if it is raised above one.
#
# Several GC values are shown within the locus rather than one. A single
# shaded pair invites the reading that a particular sequence is peculiar;
# three groups at three different GC contents make it plain that the effect
# is a property of the model rather than of those two windows.
#
# LOCUS SELECTION is scripted, not visual. Every position of a sliding window
# is scored by the summed Tm spread of its best `--n-groups` equal-GC groups;
# the highest-scoring locus is taken, candidates overlapping it are removed,
# and the next highest is taken. With --in-peaks the candidates are further
# restricted to loci centred on a MutL-AR peak, which is what connects these
# panels to the yellow bands in A and B.
#
# Windows are compared at exact GC percentages, never binned. A 200 bp window
# can only take GC values in steps of 0.5%, so treating 49.5% and 50.0% as
# equal would put a real GC difference inside the comparison and the Tm
# difference could then be attributed to GC rather than to base arrangement.
#
# The profiles use the Breslauer et al. (1986) parameter set at 50 mM Na+,
# matching Hasenauer et al. (2025), who computed Tm with this package and
# these settings.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
has_flag <- function(f) f %in% args

outdir    <- argval("--outdir", "figures")
dpi       <- as.numeric(argval("--dpi", "600"))
fmt       <- argval("--format", "tif")            # tif | jpg | none
nn_table  <- argval("--nn-table", "DNA_NN_Breslauer_1986")
Na_mM     <- as.numeric(argval("--na", "50"))
win       <- as.integer(argval("--window", "200"))
locus_kb  <- as.numeric(argval("--locus-kb", "20"))
n_groups  <- as.integer(argval("--n-groups", "3"))
n_regions <- as.integer(argval("--n-regions", "1"))
# Comma-separated to override the search, one entry per zoom panel.
region    <- argval("--region", "")

fig_w  <- as.numeric(argval("--width",  "18"))
fig_h  <- as.numeric(argval("--height", "6.5"))
# The zoom cell is given a little more width than the circular ones: it is
# the only panel whose content is not square.
zoom_w <- as.numeric(argval("--zoom-width", "1.2"))

# One letter per cell. Pass "A,B,C,C" to present two zooms as a single panel
# C, or trailing empties to leave later cells unlabelled.
letters4 <- strsplit(argval("--panel-letters", "A,B,C,D"), ",")[[1]]

# The chromosome bar in the linear panels is labelled with `genome_name`. A
# single space suppresses it: the ideogram is already labelled MutL-AR by its
# own track name, and the coordinates are on the axis. It cannot be the empty
# string, because the linear layout uses this value as the contig name for
# the karyotype and the zoom GRanges. Track data is unaffected either way,
# since the single-contig path overwrites every track's seqnames with the
# karyotype's own (plot_genome_track.R:902-911).
panelc_label <- argval("--panelc-label", " ")
# Text at the centre of the circular panels. In circular mode `genome_name`
# is used for nothing else (plot_genome_track.R:764), so emptying it removes
# the label without touching the tracks. The package name ran out under the
# rings and is already in the caption.
center_label <- argval("--center-label", "")

# Base-position ticks. plot_genome_track() falls back to 500 kb whenever the
# view is under 10 Mb, which places no tick at all inside a 20 kb panel. The
# target is four intervals across the view, matching Figure 4, which spans
# 100 to 300 kb with a tick every 50 kb.
tick_bp   <- as.numeric(argval("--tick", "0"))
n_tick_iv <- as.numeric(argval("--tick-intervals", "4"))

mark_coords    <- has_flag("--mark-coords")
in_peaks       <- has_flag("--in-peaks")
# Require BOTH windows of every shaded pair to lie inside a peak, not merely
# the locus to be centred on one. This is the version that survives the
# question "was the comparison actually made within the target region?".
pairs_in_peaks <- has_flag("--pairs-in-peaks")
if (pairs_in_peaks) in_peaks <- TRUE

leg_cex <- as.numeric(argval("--legend-cex", "0.9"))
# Transparency of the shaded bands. Kept in variables because the legend keys
# are composited at the same alpha: a solid key beside a 40%-alpha band is a
# different colour on the page and the reader has to guess which is which.
grp_alpha  <- as.numeric(argval("--group-alpha", "0.40"))
peak_alpha <- as.numeric(argval("--peak-alpha", "0.18"))
leg_pos <- local({
  v <- argval("--legend-pos", "")
  if (nzchar(v)) as.numeric(strsplit(v, ",")[[1]]) else c(0.75, 1)
})

# plot_genome_track() sizes text for one plot on a default device: axis.cex
# 0.6 circular and 0.5 linear, title.cex 0.7, label.cex 0.6. Those come out
# as a few points once this figure is scaled to a journal column.
axis_cex  <- as.numeric(argval("--axis-cex",  "1.0"))
title_cex <- as.numeric(argval("--title-cex", "1.1"))
lab_cex   <- as.numeric(argval("--label-cex", "1.0"))
panel_cex <- as.numeric(argval("--panel-cex", "1.8"))
cache     <- argval("--cache", file.path(outdir, "figure3_profiles.rds"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(TmCalculator)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

ecoli_pkg  <- "BSgenome.Ecoli.NCBI.ASM584v2"
genome_obj <- "Ecoli"
chr_name   <- "U00096.3"

if (!requireNamespace(ecoli_pkg, quietly = TRUE))
  stop("Install ", ecoli_pkg, " first; see vignette(\"genome_wide_tm_ecoli\").")
suppressPackageStartupMessages(library(ecoli_pkg, character.only = TRUE))

genome     <- base::get(genome_obj, envir = asNamespace(ecoli_pkg))
chr_length <- GenomeInfoDb::seqlengths(genome)[[chr_name]]
data(ecoli_rep_hotspots, package = "TmCalculator")

## -- Profile (cached) ------------------------------------------------------
if (file.exists(cache)) {
  message("using cached profile: ", cache)
  Tm <- readRDS(cache)
} else {
  message("computing the profile over ", win, " bp windows")
  bins <- make_genomiccoord(bsgenome = ecoli_pkg, chromosomes = chr_name,
                            window = win, slide = win, start = 1,
                            end = chr_length, strand = "+", verbose = FALSE)
  gr <- to_genomic_ranges_fast(list(pkg_name = ecoli_pkg, seq = bins))
  nn <- tm_calculate(gr, method = "tm_nn", nn_table = nn_table, Na = Na_mM)$gr
  Tm <- as.data.frame(nn[, c("Tm", "GC")])
  saveRDS(Tm, cache)
  message("cached to ", cache)
}

## -- Windows ---------------------------------------------------------------
# Full-length windows only. A window shortened at the chromosome end, or one
# that contained an ambiguous base, has a GC denominator smaller than `win`
# and would not be comparable with the rest on an exact-GC basis.
W <- Tm[!is.na(Tm$Tm) & !is.na(Tm$GC) & Tm$width == win, ]
W <- W[order(W$start), ]
n_loc <- max(2L, as.integer(round(locus_kb * 1000 / win)))
if (nrow(W) <= n_loc) stop("locus is wider than the profile")

# MutL-AR peaks as an interval set. Only coordinates are used, never sequence
# names: matching on names would fail silently if the peak table calls the
# contig by its accession and the window table by the package name.
peaks_df <- as.data.frame(ecoli_rep_hotspots$all_peaks_IP_mutH)
peaks_ir <- IRanges::IRanges(start = as.numeric(peaks_df$start),
                             end   = as.numeric(peaks_df$end))
W$in_peak <- IRanges::overlapsAny(
  IRanges::IRanges(start = W$start, end = W$end), peaks_ir)

## -- Equal-GC groups -------------------------------------------------------
# The `k` GC values inside a locus whose windows are furthest apart in Tm.
# Each group contributes its two extreme windows; shading every window at
# that GC would fill the panel with stripes and hide the tracks.
top_groups <- function(d, k, peak_only = FALSE) {
  idx <- if (peak_only) which(d$in_peak) else seq_len(nrow(d))
  if (length(idx) < 2L) return(list())
  g <- split(idx, d$GC[idx])
  g <- g[lengths(g) >= 2L]
  if (!length(g)) return(list())
  sp <- vapply(g, function(ix) diff(range(d$Tm[ix])), numeric(1))
  o  <- order(sp, decreasing = TRUE)[seq_len(min(k, length(g)))]
  lapply(o, function(j) {
    ix <- g[[j]]
    lo <- ix[which.min(d$Tm[ix])]; hi <- ix[which.max(d$Tm[ix])]
    list(gc = as.numeric(names(g)[j]), spread = unname(sp[j]),
         n = length(ix), tm_lo = d$Tm[lo], tm_hi = d$Tm[hi],
         rows = c(lo, hi))
  })
}
score_of <- function(d)
  sum(vapply(top_groups(d, n_groups, pairs_in_peaks),
             function(g) g$spread, numeric(1)))

## -- Choose the loci -------------------------------------------------------
# Wrapped in a function rather than written as top-level if/else. Rscript
# parses a script expression by expression, and a top-level `if` block that
# ends before its `else` is a parse error that only appears after everything
# above it has already run.
make_reg <- function(rows) {
  d  <- W[rows, ]
  gs <- top_groups(d, n_groups, pairs_in_peaks)
  span <- max(d$end) - min(d$start) + 1
  list(d = d, gs = gs, span = span,
       zoom = sprintf("%s:%d-%d", chr_name, min(d$start), max(d$end)))
}

pick_regions <- function() {
  if (nzchar(region)) {
    return(lapply(strsplit(region, ",")[[1]], function(r) {
      rr  <- as.numeric(strsplit(sub("^.*:", "", trimws(r)), "-")[[1]])
      sel <- which(W$start >= rr[1] & W$end <= rr[2])
      if (length(sel) < 2L) stop("fewer than two full windows in ", r)
      reg <- make_reg(sel)
      if (!length(reg$gs)) stop("no two windows share a GC value in ", r)
      reg
    }))
  }

  starts <- seq_len(nrow(W) - n_loc + 1L)
  # A locus is `n_loc` CONSECUTIVE ROWS of W, and W has had short and
  # ambiguous windows dropped. Where a run of windows was removed, consecutive
  # rows are not adjacent on the chromosome, and such a locus would be drawn
  # with an axis spanning far more than --locus-kb: the panel would silently
  # stop being a zoom. Those candidates are discarded rather than trusted.
  contig <- (W$start[seq(n_loc, nrow(W))] - W$start[starts]) ==
    (n_loc - 1L) * win
  if (!all(contig))
    message(sprintf("%d of %d candidate loci span a gap and were dropped",
                    sum(!contig), length(contig)))
  starts <- starts[contig]

  if (in_peaks) {
    # Loci whose midpoint falls inside a peak. Requiring full coverage would
    # exclude almost everything, since the peaks are far narrower than 20 kb.
    mid    <- starts + (n_loc %/% 2L)
    starts <- starts[W$in_peak[pmin(mid, nrow(W))]]
    if (!length(starts))
      stop("no candidate locus is centred on a MutL-AR peak; widen ",
           "--locus-kb or drop --in-peaks")
    message("scanning ", length(starts), " loci centred on MutL-AR peaks")
  } else {
    message("scanning ", length(starts), " loci")
  }

  sc    <- vapply(starts, function(i) score_of(W[i:(i + n_loc - 1L), ]),
                  numeric(1))
  avail <- rep(TRUE, length(starts))
  out   <- vector("list", n_regions)
  for (k in seq_len(n_regions)) {
    if (!any(avail))
      stop("only ", k - 1L, " non-overlapping loci available; lower ",
           "--n-regions or --locus-kb")
    j  <- which(avail)[which.max(sc[avail])]
    i0 <- starts[j]
    out[[k]] <- make_reg(i0:(i0 + n_loc - 1L))
    # Loci are drawn as separate panels, so they must not overlap: two views
    # of the same sequence would be presented as two independent examples.
    avail <- avail & (abs(starts - i0) >= n_loc)
  }
  out
}

regs <- pick_regions()
for (r in regs)
  if (!nzchar(region) && r$span > locus_kb * 1000 * 1.001)
    stop(sprintf("a locus spans %.1f kb but --locus-kb is %g; the windows ",
                 r$span / 1000, locus_kb),
         "are not contiguous. Inspect the profile for dropped windows.")

if (tick_bp <= 0) {
  cand    <- c(100, 200, 500, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5)
  tick_bp <- cand[which.min(abs(regs[[1]]$span / cand - n_tick_iv))]
}

# Amber, sea green, purple, teal: distinguishable from each other and from
# the blue GC and red Tm traces, and still distinct in greyscale. The same
# palette is reused in each zoom panel, which is unambiguous because every
# panel carries its own legend.
grp_cols <- c("#B7791F", "#2E8B57", "#7D3C98", "#117A8B", "#A93226")

## -- Tracks ----------------------------------------------------------------
label <- data.frame(
  seqnames = ecoli_pkg,
  start    = c(3925804, 1590777),
  end      = c(3925804, 1590777),
  label    = c("ori", "dif"))

# "Tm" rather than "Melting temp": the long form is truncated at the width
# the circular panel gives its legend. Black rather than navy for the peaks:
# they are narrow marks inside a grey bar and lose contrast when scaled down.
tracks <- list(
  list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
       col = "black", bg.col = "grey", name = "MutL-AR",
       legend_font_col = "black", ideogram = TRUE, height = 0.5),
  list(type = "line", data = Tm, value_col = "GC", name = "GC content",
       col = "#4A90E2", legend_font_col = "#4A90E2"),
  list(type = "line", data = Tm, value_col = "Tm", name = "Tm",
       col = "#E06666", legend_font_col = "#E06666", height = 2),
  list(type = "line", data = ecoli_rep_hotspots$bins_rep, value_col = "count",
       name = "Microsatellites", col = "#2ECC71", legend_font_col = "#2ECC71"),
  list(type = "line", data = ecoli_rep_hotspots$bins_cru, value_col = "count",
       name = "Cruciform", col = "#3B3E6B", legend_font_col = "#3B3E6B"),
  list(data = ecoli_rep_hotspots$ssdna, name = "ssDNA",
       col = "#8E44AD", legend_font_col = "#8E44AD"),
  list(type = "line", data = ecoli_rep_hotspots$bins_gatc, value_col = "count",
       name = "GATC sites", col = "#D35400", legend_font_col = "#D35400"),
  list(type = "highlight", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
       col = "#F1C40F", alpha = peak_alpha))

zoom_tracks <- function(reg) c(
  list(
    # The MutL-AR layer is drawn here too, in the same black, so the three
    # panels agree on what a peak looks like. As the ideogram it also names
    # the bar, which is why the panel needs no contig label of its own.
    list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
         col = "black", bg.col = "grey", name = "MutL-AR",
         legend_font_col = "black", ideogram = TRUE, height = 0.5),
    # "GC" rather than "GC content": the linear panel writes track names down
    # the left margin, where they compete with the axis numbers for width.
    list(type = "line", data = Tm, value_col = "GC", name = "GC",
         col = "#4A90E2", legend_font_col = "#4A90E2", height = 1),
    list(type = "line", data = Tm, value_col = "Tm", name = "Tm",
         col = "#E06666", legend_font_col = "#E06666", height = 1.4)),
  list(list(type = "highlight", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
            col = "#F1C40F", alpha = peak_alpha)),
  # One highlight per GC group, coordinates taken from the search, so a band
  # cannot drift away from the numbers quoted in the caption.
  lapply(seq_along(reg$gs), function(i)
    list(type = "highlight",
         data = reg$d[reg$gs[[i]]$rows, c("seqnames", "start", "end")],
         col  = grp_cols[i], alpha = grp_alpha)))

coord_labels <- function(reg) {
  if (!mark_coords) return(NULL)
  p <- do.call(rbind, lapply(reg$gs, function(g) reg$d[g$rows, ]))
  data.frame(seqnames = panelc_label, start = p$start, end = p$start,
             label = format(p$start, big.mark = ","),
             stringsAsFactors = FALSE)
}

draw <- function() {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  graphics::layout(matrix(seq_len(2L + n_regions), nrow = 1),
                   widths = c(1, 1, rep(zoom_w, n_regions)))

  ## ---- A: whole genome -------------------------------------------------
  graphics::par(mar = c(2, 2, 3, 2), cex = 1)
  plot_genome_track(genome_name = center_label, genome_size = chr_length,
                    track_list = tracks, circular = TRUE, label = label,
                    legend.cex = leg_cex, legend.position = leg_pos,
                    axis.cex = axis_cex, title.cex = title_cex,
                    label.cex = lab_cex)
  graphics::mtext(letters4[1], side = 3, line = 0.6, adj = 0,
                  cex = panel_cex, font = 2)

  ## ---- B: two arcs -----------------------------------------------------
  # The legend is drawn once, on A. Repeating it takes width from the arcs,
  # which are the only thing B adds over A.
  graphics::par(mar = c(2, 2, 3, 2), cex = 1)
  plot_genome_track(genome_name = center_label, genome_size = chr_length,
                    track_list = tracks, circular = TRUE,
                    zoom = c(paste0(chr_name, ":100000-500000"),
                             paste0(chr_name, ":3600000-4500000")),
                    legend.show = FALSE, axis.cex = axis_cex,
                    title.cex = title_cex, label.cex = lab_cex)
  graphics::mtext(letters4[2], side = 3, line = 0.6, adj = 0,
                  cex = panel_cex, font = 2)

  ## ---- C, D: the loci --------------------------------------------------
  for (k in seq_along(regs)) {
    reg <- regs[[k]]
    graphics::par(mar = c(3, 2, 3, 2), cex = 1)
    plot_genome_track(genome_name = panelc_label, genome_size = chr_length,
                      track_list = zoom_tracks(reg), zoom = reg$zoom,
                      track.gap = 0.06, legend.show = FALSE,
                      axis.cex = axis_cex, title.cex = title_cex,
                      label.cex = lab_cex, base.tick.dist = tick_bp,
                      base.tick.units = TRUE, label = coord_labels(reg))
    lt <- if (length(letters4) >= k + 2L) letters4[k + 2L] else ""
    if (nzchar(lt))
      graphics::mtext(lt, side = 3, line = 0.6, adj = 0,
                      cex = panel_cex, font = 2)

    # legend() multiplies its cex by par("cex"), and plot.new() inside
    # karyoploteR resets par("cex") to the value layout() chose for a
    # multi-panel figure. Setting it back is what makes --legend-cex mean the
    # same thing here as in panel A.
    graphics::par(cex = 1)
    # The track key is drawn only in the first zoom panel; the group keys are
    # per panel, because the palette is reused and each panel's colours stand
    # for different GC values. ASCII only: pdf() and tiff() drop characters
    # they cannot map through mbcsToSbcs, silently in the raster output.
    trk_lab <- if (k == 1L) c("GC", "Tm", "MutL-AR peak") else character(0)
    trk_fil <- if (k == 1L)
      c("#4A90E2", "#E06666",
        grDevices::adjustcolor("#F1C40F", alpha.f = peak_alpha)) else character(0)
    trk_txt <- if (k == 1L) c("#4A90E2", "#E06666", "#B7950B") else character(0)
    gi <- seq_along(reg$gs)
    graphics::legend(
      "topright", bty = "n", cex = leg_cex, border = NA,
      legend = c(trk_lab,
                 vapply(gi, function(i)
                   sprintf("GC %.1f%%:  Tm %.1f-%.1f", reg$gs[[i]]$gc,
                           reg$gs[[i]]$tm_lo, reg$gs[[i]]$tm_hi),
                   character(1))),
      fill     = c(trk_fil,
                   grDevices::adjustcolor(grp_cols[gi], alpha.f = grp_alpha)),
      text.col = c(trk_txt, grp_cols[gi]))
  }
}

## -- Write -----------------------------------------------------------------
stem <- file.path(outdir, "figure3_ecoli_circos_zoom")

grDevices::pdf(paste0(stem, ".pdf"), width = fig_w, height = fig_h,
               useDingbats = FALSE)
draw(); grDevices::dev.off()

if (identical(fmt, "tif")) {
  # LZW is lossless. "Insufficient resolution" applies to a 300 dpi JPEG of a
  # figure made of thin concentric lines; JPEG softens exactly that.
  grDevices::tiff(paste0(stem, ".tif"), width = fig_w, height = fig_h,
                  units = "in", res = dpi, compression = "lzw",
                  type = if (capabilities("cairo")) "cairo" else "quartz")
  draw(); grDevices::dev.off()
} else if (identical(fmt, "jpg")) {
  grDevices::jpeg(paste0(stem, ".jpg"), width = fig_w, height = fig_h,
                  units = "in", res = dpi, quality = 95)
  draw(); grDevices::dev.off()
}

## -- Report ----------------------------------------------------------------
for (k in seq_along(regs)) {
  reg <- regs[[k]]
  lt  <- if (length(letters4) >= k + 2L) letters4[k + 2L] else paste0("#", k)
  cat(sprintf("\nPanel %s: %s\n", lt, reg$zoom))
  cat(sprintf("  x axis %s to %s bp (%.1f kb, %d windows, ticks every %s bp)\n",
              format(min(reg$d$start), big.mark = ","),
              format(max(reg$d$end),   big.mark = ","),
              reg$span / 1000, nrow(reg$d),
              format(tick_bp, big.mark = ",", scientific = FALSE)))
  cat(sprintf("  windows overlapping a MutL-AR peak: %d of %d%s\n",
              sum(reg$d$in_peak), nrow(reg$d),
              if (pairs_in_peaks) "   (pairs required inside peaks)"
              else if (in_peaks) "   (locus centred on a peak)" else ""))
  for (i in seq_along(reg$gs)) {
    g <- reg$gs[[i]]; p <- reg$d[g$rows, ]
    cat(sprintf("\n  group %d  colour %s  GC = %.1f%%  (%d windows at this GC)\n",
                i, grp_cols[i], g$gc, g$n))
    for (j in 1:2)
      cat(sprintf("    %s:%s-%s   Tm = %.2f C%s\n", chr_name,
                  format(p$start[j], big.mark = ","),
                  format(p$end[j],   big.mark = ","), p$Tm[j],
                  if (p$in_peak[j]) "   [in MutL-AR peak]" else ""))
    cat(sprintf("    difference at identical GC: %.2f C\n", g$spread))
  }
}
cat(sprintf("\nLargest such difference genome-wide: %.2f C\n",
            max(vapply(split(W$Tm, W$GC),
                       function(v) if (length(v) > 1L) diff(range(v)) else 0,
                       numeric(1)))))

info <- file.info(list.files(outdir, pattern = "^figure3_.*\\.(pdf|tif|jpg)$",
                             full.names = TRUE))
cat("\nWritten:\n")
for (i in seq_len(nrow(info)))
  cat(sprintf("  %-48s %6.1f MB\n", rownames(info)[i], info$size[i] / 1e6))
cat(sprintf("\n%.0f x %.0f pixels at %g dpi\n", fig_w * dpi, fig_h * dpi, dpi))
