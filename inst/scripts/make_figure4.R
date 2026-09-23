#!/usr/bin/env Rscript
# ===========================================================================
# make_figure4.R -- publication-resolution export of Figure 4
#
#   Tm profiles of the 0.1-0.3 Mb interval of the E. coli K-12 MG1655
#   chromosome computed at 50, 100, 200 and 500 bp.
#
#   Rscript inst/scripts/make_figure4.R                 # PDF + 600 dpi TIFF
#   Rscript inst/scripts/make_figure4.R --dpi 300       # smaller raster
#   Rscript inst/scripts/make_figure4.R --gap 0.05 --axis-cex 0.45
#
# The four profiles are cached to an .rds after the first run. Laying out a
# five-track panel takes several attempts and recomputing 172,000 windows for
# each of them wastes minutes; delete the cache file to force a recompute.
#
# Both a vector PDF and a raster TIFF are written. The PDF is the better
# source for a journal that accepts vector art, since the tracks are line
# work and stay sharp at any magnification; the TIFF exists because some
# production pipelines require raster, and is written at 600 dpi with LZW
# compression rather than the 300 dpi minimum, because fine single-window
# fluctuations at 50 bp are the point of the panel and are the first thing
# lost to downsampling.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}

outdir   <- argval("--outdir", "figures")
dpi      <- as.numeric(argval("--dpi", "600"))
fig_w    <- as.numeric(argval("--width", "7.5"))    # inches; MDPI full width
fig_h    <- as.numeric(argval("--height", "6.0"))
gap      <- as.numeric(argval("--gap", "0.03"))
axis_cex <- as.numeric(argval("--axis-cex", "0.55"))
# karyoploteR reserves 50 plot units above and below the data. That is right
# for a single track and leaves a wide empty band once five tracks share the
# panel, so the margins are trimmed here rather than cropped afterwards.
#
# The top can be cut hard because the panel title is suppressed. The bottom
# cannot: the chromosome name and the base-position ticks live there, and
# trimming it to match the top clips them.
topmar   <- as.numeric(argval("--topmargin", "10"))
botmar   <- as.numeric(argval("--bottommargin", "40"))
outmar   <- as.numeric(argval("--outmargin", "8"))
# Ticks every 50 kb label the 0.1-0.3 Mb interval on the axis itself, which
# is where the region is now stated: the caption repeats it, but a figure
# should not depend on its caption to say which coordinates it shows.
tick_bp  <- as.numeric(argval("--tick", "50000"))
region   <- argval("--region", "U00096.3:100000-300000")
cache    <- argval("--cache", file.path(outdir, "figure4_profiles.rds"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(TmCalculator)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
})

ecoli_pkg  <- "BSgenome.Ecoli.NCBI.ASM584v2"
genome_obj <- "Ecoli"
chr_name   <- "U00096.3"
window_sizes <- c(50L, 100L, 200L, 500L)

if (!requireNamespace(ecoli_pkg, quietly = TRUE))
  stop("Install ", ecoli_pkg, " first; see vignette(\"genome_wide_tm_ecoli\").")
suppressPackageStartupMessages(library(ecoli_pkg, character.only = TRUE))

genome     <- base::get(genome_obj, envir = asNamespace(ecoli_pkg))
chr_length <- GenomeInfoDb::seqlengths(genome)[[chr_name]]
data(ecoli_rep_hotspots, package = "TmCalculator")

## -- Profiles (cached) -----------------------------------------------------
if (file.exists(cache)) {
  message("using cached profiles: ", cache)
  sens_dfs <- readRDS(cache)
} else {
  message("computing profiles at ", paste(window_sizes, collapse = ", "), " bp")
  sens_dfs <- lapply(window_sizes, function(w) {
    bins <- make_genomiccoord(
      bsgenome = ecoli_pkg, chromosomes = chr_name,
      window = w, slide = w, start = 1, end = chr_length,
      strand = "+", verbose = FALSE)
    gr <- to_genomic_ranges_fast(list(pkg_name = ecoli_pkg, seq = bins))
    tm <- tm_calculate(gr, method = "tm_nn",
                       nn_table = "DNA_NN_Breslauer_1986", Na = 50)$gr
    as.data.frame(tm[, c("Tm", "GC")])
  })
  saveRDS(sens_dfs, cache)
  message("cached to ", cache)
}

## -- Tracks ----------------------------------------------------------------
# Dark to light with increasing window size, so that the visual weight
# follows the amount of detail rather than the alphabet.
scale_cols_tm <- c("#7B241C", "#CB4335", "#EC7063", "#F5B7B1")

tracks <- c(
  list(list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
            col = "#2C3E50", bg.col = "grey", name = "MutL-AR",
            legend_font_col = "#2C3E50", ideogram = TRUE, height = 0.6)),
  lapply(seq_along(window_sizes), function(i)
    list(type = "line", data = sens_dfs[[i]], value_col = "Tm",
         name = paste0("Tm ", window_sizes[i], " bp"),
         col = scale_cols_tm[i], legend_font_col = scale_cols_tm[i],
         height = 1.2)),
  list(list(type = "highlight",
            data = ecoli_rep_hotspots$all_peaks_IP_mutH,
            col = "#F1C40F", alpha = 0.18))
)

draw <- function() {
  plot_genome_track(
    ## The MutL-AR track is drawn as the ideogram, and karyoploteR labels the
    ## ideogram with the chromosome name it was given. Passing the genome
    ## package name there produced a label that was both uninformative for
    ## this panel and too long to fit; the band is the MutL-AR peaks, so that
    ## is what it is called. The chromosome and assembly belong in the
    ## caption. Only the start and end of `zoom` are used for a single-contig
    ## genome, so the name has no other effect.
    genome_name = "MutL-AR",
    genome_size = chr_length,
    track_list  = tracks,
    zoom        = region,
    track.gap   = gap,
    axis.cex    = axis_cex,
    ## No panel title: the genome and the interval belong in the caption.
    main        = "",
    ## The left-hand track labels already name every track, so a legend
    ## repeats them and consumes width that the tracks can use.
    legend.show = FALSE,
    base.tick.dist = tick_bp,
    plot.params = list(topmargin      = topmar,
                       bottommargin   = botmar,
                       data1outmargin = outmar)
  )
}

## -- Write -----------------------------------------------------------------
stem <- file.path(outdir, "figure4_window_size")

grDevices::pdf(paste0(stem, ".pdf"), width = fig_w, height = fig_h,
               useDingbats = FALSE)
draw(); grDevices::dev.off()

# LZW is lossless and is what journals expect for line art; JPEG compression
# inside a TIFF would soften exactly the single-window detail this panel is
# meant to show.
grDevices::tiff(paste0(stem, ".tif"), width = fig_w, height = fig_h,
                units = "in", res = dpi, compression = "lzw",
                type = if (capabilities("cairo")) "cairo" else "quartz")
draw(); grDevices::dev.off()

info <- file.info(paste0(stem, c(".pdf", ".tif")))
cat("\nWritten:\n")
for (i in seq_len(nrow(info)))
  cat(sprintf("  %-40s %6.1f MB\n", rownames(info)[i], info$size[i] / 1e6))
cat(sprintf("\n%.0f x %.0f pixels at %g dpi (%.1f x %.1f in)\n",
            fig_w * dpi, fig_h * dpi, dpi, fig_w, fig_h))
cat("MDPI asks for at least 1000 px on the long side; this is",
    format(round(max(fig_w, fig_h) * dpi), big.mark = ","), "px.\n")
