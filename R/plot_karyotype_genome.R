#' Plot multi-track karyotype genome view using karyoploteR
#'
#' @description
#' Generate a linear genome plot with flexible track definitions using the
#' karyoploteR package. Each track is defined in a list, allowing dynamic
#' visualization of genomic data such as GC content, Tm, peak regions, and
#' other multi-omics signals. The input format mirrors
#' \code{\link{plot_circos_genome}} for consistency across visualization modes.
#'
#' @param genome_assembly Character string specifying the genome assembly (e.g.,
#'   "hg19", "mm10"), or any genome specification accepted by
#'   \code{\link[karyoploteR]{plotKaryotype}} (GRanges, Seqinfo, etc.).
#' @param track_list List. Each element defines one track. Supported fields
#'   (identical to \code{plot_circos_genome} where applicable):
#'   \itemize{
#'     \item \strong{type}: Character. Drawing type - one of \code{"points"},
#'       \code{"line"}, \code{"bars"}, or \code{"rect"}. Default: \code{"points"}.
#'     \item \strong{data}: GRanges or data.frame with columns seqnames/chr,
#'       start, end, and the value column.
#'     \item \strong{min_bar_width}: Numeric (0-1). For \code{bars} type only.
#'       Minimum bar width as a fraction of chromosome length. Needed when
#'       sequences are short (e.g. primers / probes) whose actual genomic
#'       footprint is sub-pixel at whole-chromosome scale. Default: \code{0.002}
#'       (0.2\% of chr length, ~500 kb on chr1). Set to \code{0} to use
#'       actual range widths only.
#'     \item \strong{value_col}: Character. Metadata column containing numeric
#'       values (required for \code{points}, \code{line}, \code{bars}).
#'     \item \strong{ylim}: Numeric vector of length 2. Y-axis limits. Inferred
#'       from data if \code{NULL}.
#'     \item \strong{col}: Color for data and legend symbol. Auto-assigned if
#'       \code{NULL}.
#'     \item \strong{name}: Character. Track label shown in legend.
#'     \item \strong{r0}, \strong{r1}: Numeric (0-1). Vertical position of the
#'       track within the data panel. Auto-computed if \code{NULL}.
#'     \item \strong{bg.col}: Track background fill color. Default: \code{"white"}.
#'     \item \strong{bg.border}: Track background border color. Default: \code{NA}.
#'     \item \strong{bg.lwd}: Track background border line width. Default: \code{0.5}.
#'     \item \strong{lwd}: Line width (for \code{line} type). Default: \code{1}.
#'     \item \strong{cex}: Point size (for \code{points} type). Default: \code{0.5}.
#'     \item \strong{pch}: Point shape (for \code{points} type). Default: \code{16}.
#'     \item \strong{ytop}, \strong{ybottom}: Top/bottom y values for
#'       \code{rect} tracks within \code{ylim} range. Defaults: \code{0.98} /
#'       \code{0.02}.
#'     \item \strong{border}: Border color for \code{rect} and \code{bars}
#'       tracks. Default: \code{NA}.
#'   }
#' @param chromosomes Character vector of chromosomes to plot, or \code{"auto"}
#'   to derive from the data, or any value accepted by karyoploteR (e.g.,
#'   \code{"canonical"}, \code{"autosomal"}). Default: \code{"auto"}.
#' @param plot_type Integer. karyoploteR plot type (1 = horizontal chromosomes
#'   stacked vertically, 2 = two panels, etc.). Default: \code{1}.
#' @param track.margin Numeric. Fractional gap between adjacent tracks.
#'   Default: \code{0.05}.
#' @param zoom A GRanges or character string (e.g., \code{"chr1:1e6-2e6"})
#'   specifying a region to zoom into. \code{NULL} shows full chromosomes.
#' @param tick_dist Numeric. Distance between x-axis tick marks in base pairs.
#'   Default: \code{10000000}.
#' @param axis_cex Numeric. Axis label text size. Default: \code{0.6}.
#' @param chr_cex Numeric. Chromosome name text size. Default: \code{0.8}.
#' @param legend.show Logical. Whether to display the legend. Default: \code{TRUE}.
#' @param legend.position Character. Legend position passed to
#'   \code{\link[graphics]{legend}} (e.g., \code{"topright"}, \code{"topleft"}).
#'   Default: \code{"topright"}.
#' @param legend.cex Numeric. Legend label text size. Default: \code{0.7}.
#' @param legend.bty Character. Legend box type (\code{"n"} for no box).
#'   Default: \code{"n"}.
#' @param label Optional GRanges whose metadata column \code{label.column}
#'   contains text labels to plot at each range.
#' @param label.column Character or integer. The metadata column in \code{label}
#'   that holds label text.
#' @param label.cex Numeric. Label text size. Default: \code{0.6}.
#' @param title_name Character. Main plot title. \code{NULL} omits the title.
#' @param title.cex Numeric. Title text size. Default: \code{0.8}.
#'
#' @return Invisibly returns the karyoploteR \code{KaryoPlot} object, enabling
#'   further customisation with karyoploteR functions.
#'
#' @details
#' Track vertical positions (\code{r0}/\code{r1}) are auto-computed by dividing
#' the data panel equally among tracks, separated by \code{track.margin}. Supply
#' explicit \code{r0}/\code{r1} values in individual track list entries to
#' override the automatic layout.
#'
#' The \code{track_list} format is intentionally identical to
#' \code{plot_circos_genome} so that the same track definitions can be reused
#' across circular and linear genome representations.
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' library(karyoploteR)
#'
#' # -- Sample data -----------------------------------------------------------
#' # -- Typical workflow: output from tm_calculate() already has Tm AND GC ------
#' # tm_calculate() stores both columns in result$gr, so one GRanges object
#' # is sufficient for all tracks.
#' #
#' # result <- tm_calculate(my_seqs, method = "tm_nn")
#' # gr_tm  <- result$gr   # has $Tm (deg C) and $GC (%) columns
#'
#' # -- Example 1: Tm points + GC bars from the same GRanges ------------------
#' set.seed(42)
#' gr_tm <- GRanges(
#'   seqnames = c(rep("chr1", 50), rep("chr2", 30)),
#'   ranges = IRanges(
#'     start = c(sort(sample(1:249e6, 50)), sort(sample(1:243e6, 30))),
#'     width = sample(50:200, 80, replace = TRUE)
#'   ),
#'   Tm = runif(80, 55, 85),
#'   GC = runif(80, 30, 70)   # gc() returns percentages (0-100)
#' )
#' gr_peaks <- GRanges(
#'   seqnames = c("chr1", "chr1", "chr2"),
#'   ranges = IRanges(start = c(10e6, 100e6, 50e6), width = c(1e6, 2e6, 500e3))
#' )
#'
#' track_list <- list(
#'   list(type = "points", data = gr_tm, value_col = "Tm",
#'        col = "#e41a1c", name = "Tm (deg C)"),
#'   list(type = "bars",   data = gr_tm, value_col = "GC",
#'        col = "#377eb8", name = "GC (%)"),
#'   list(type = "rect",   data = gr_peaks,
#'        col = "#4daf4a", name = "Peaks")
#' )
#' plot_karyotype_genome(
#'   genome_assembly = "hg19",
#'   track_list      = track_list,
#'   title_name      = "Multi-omics - hg19"
#' )
#'
#' # -- Example 2: line track, zoom -------------------------------------------
#' plot_karyotype_genome(
#'   genome_assembly = "hg19",
#'   track_list      = list(
#'     list(type = "line", data = gr_tm, value_col = "Tm", col = "#984ea3", name = "GC (%)")
#'   ),
#'   zoom            = "chr1:1000000-50000000",
#'   title_name      = "Tm - chr1 zoom"
#' )
#' }
#'
#' @importFrom karyoploteR plotKaryotype kpDataBackground kpAxis kpPoints
#'   kpLines kpBars kpRect kpAddBaseNumbers kpText
#' @importFrom GenomicRanges GRanges seqnames makeGRangesFromDataFrame mcols
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics legend title
#'
#' @encoding UTF-8
#' @author Junhui Li
#' @export

plot_karyotype_genome <- function(
    genome_assembly,
    track_list,
    
    # layout
    chromosomes  = "auto",
    plot_type    = 1,
    track.margin = 0.05,
    zoom         = NULL,
    
    # axis
    tick_dist = 10000000,
    axis_cex  = 0.6,
    chr_cex   = 0.8,
    
    # legend
    legend.show     = TRUE,
    legend.position = "topright",
    legend.cex      = 0.7,
    legend.bty      = "n",
    
    # labels
    label        = NULL,
    label.column = NULL,
    label.cex    = 0.6,
    
    # title
    title_name = NULL,
    title.cex  = 0.8
) {
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # -- Internal defaults ------------------------------------------------------
  default_bg.col    <- "white"
  default_bg.border <- NA
  default_bg.lwd    <- 0.5
  default_lwd       <- 1
  default_cex       <- 0.5
  default_pch       <- 16
  
  # -- Default track colors ---------------------------------------------------
  palette      <- c("#67000d", "#fc9272")
  default_cols <- grDevices::colorRampPalette(palette)(length(track_list))
  
  # -- Helper: coerce data to GRanges ----------------------------------------
  .to_gr <- function(data, track_idx) {
    if (inherits(data, "GRanges")) return(data)
    if (is.data.frame(data)) {
      return(GenomicRanges::makeGRangesFromDataFrame(data, keep.extra.columns = TRUE))
    }
    stop(sprintf("Track %d: 'data' must be a GRanges object or data.frame.", track_idx))
  }
  
  # -- Determine chromosomes -------------------------------------------------
  if (identical(chromosomes, "auto")) {
    all_chrs <- unique(unlist(lapply(seq_along(track_list), function(i) {
      trk <- track_list[[i]]
      if (!is.null(trk$data))
        as.character(GenomicRanges::seqnames(.to_gr(trk$data, i)))
      else NULL
    })))
    chromosomes <- if (length(all_chrs) > 0) sort(unique(all_chrs)) else "all"
  }
  
  # -- Auto-compute r0 / r1 -------------------------------------------------
  n  <- length(track_list)
  th <- (1 - track.margin * (n - 1)) / n  # track height as fraction of panel
  
  for (i in seq_along(track_list)) {
    if (is.null(track_list[[i]]$r0))
      track_list[[i]]$r0 <- (i - 1) * (th + track.margin)
    if (is.null(track_list[[i]]$r1))
      track_list[[i]]$r1 <- track_list[[i]]$r0 + th
  }
  
  # -- Initialise karyotype --------------------------------------------------
  kp_args <- list(
    genome      = genome_assembly,
    plot.type   = plot_type,
    chromosomes = chromosomes,
    cex         = chr_cex
  )
  if (!is.null(zoom)) kp_args$zoom <- zoom
  
  kp <- do.call(karyoploteR::plotKaryotype, kp_args)
  
  # -- Draw tracks -----------------------------------------------------------
  for (i in seq_along(track_list)) {
    trk <- track_list[[i]]
    
    trk$type      <- trk$type      %||% "points"
    trk$col       <- trk$col       %||% default_cols[i]
    trk$bg.col    <- trk$bg.col    %||% default_bg.col
    trk$bg.border <- trk$bg.border %||% default_bg.border
    trk$bg.lwd    <- trk$bg.lwd    %||% default_bg.lwd
    trk$lwd       <- trk$lwd       %||% default_lwd
    trk$cex       <- trk$cex       %||% default_cex
    trk$pch       <- trk$pch       %||% default_pch
    
    data_gr <- .to_gr(trk$data, i)
    # r0/r1 set by pre-computation loop; guard against NULL just in case
    r0 <- trk$r0 %||% 0
    r1 <- trk$r1 %||% 1
    
    # -- validate value_col --
    if (trk$type %in% c("points", "line", "bars")) {
      if (is.null(trk$value_col))
        stop(sprintf("Track %d (type='%s'): 'value_col' must be provided.", i, trk$type))
      if (!(trk$value_col %in% names(GenomicRanges::mcols(data_gr))))
        stop(sprintf("Track %d: column '%s' not found in data.", i, trk$value_col))
    }
    
    # -- ylim --
    if (is.null(trk$ylim)) {
      if (trk$type == "bars") {
        # Bar charts must include 0 in the y-range so that bar heights are
        # proportional to the values.  Using range() gives ymin = min(vals),
        # which pushes y0 = 0 below the visible panel and makes every bar the
        # same height (clipped to the track bottom).
        trk$ylim <- c(0, max(GenomicRanges::mcols(data_gr)[[trk$value_col]],
                             na.rm = TRUE))
      } else if (trk$type %in% c("points", "line")) {
        trk$ylim <- range(GenomicRanges::mcols(data_gr)[[trk$value_col]], na.rm = TRUE)
      } else {
        trk$ylim <- c(0, 1)
      }
    }
    ymin <- trk$ylim[1]
    ymax <- trk$ylim[2]
    
    # -- background --
    karyoploteR::kpDataBackground(kp, r0 = r0, r1 = r1, color = trk$bg.col)
    
    # -- per-track y-axis --
    karyoploteR::kpAxis(kp, ymin = ymin, ymax = ymax,
                        r0 = r0, r1 = r1,
                        cex = axis_cex, numticks = 3)
    
    # -- extract values --
    vals <- if (!is.null(trk$value_col))
      GenomicRanges::mcols(data_gr)[[trk$value_col]]
    else NULL
    
    # -- draw data --
    if (trk$type == "points") {
      karyoploteR::kpPoints(kp,
                            data = data_gr, y = vals,
                            ymin = ymin, ymax = ymax,
                            r0 = r0, r1 = r1,
                            col = trk$col, pch = trk$pch, cex = trk$cex)
      
    } else if (trk$type == "line") {
      # sort by position for a connected line
      ord     <- order(as.character(GenomicRanges::seqnames(data_gr)),
                       GenomicRanges::start(data_gr))
      data_gr <- data_gr[ord]
      vals    <- vals[ord]
      karyoploteR::kpLines(kp,
                           data = data_gr, y = vals,
                           ymin = ymin, ymax = ymax,
                           r0 = r0, r1 = r1,
                           col = trk$col, lwd = trk$lwd)
      
    } else if (trk$type == "bars") {
      # Short sequences (primers / probes) are typically 20-300 bp wide on a
      # chromosome hundreds of Mb long.  kpBars draws rectangles from start to
      # end, so each bar is < 1 pixel wide and invisible.
      # Fix: expand every bar to at least `min_bar_width` × chr_length,
      # centred on the range midpoint.
      min_frac <- trk$min_bar_width %||% 0.002   # default 0.2% of chr length
      chr_vec  <- as.character(GenomicRanges::seqnames(data_gr))
      chr_lens <- kp$chromosome.lengths[chr_vec]
      chr_lens[is.na(chr_lens)] <- max(kp$chromosome.lengths, na.rm = TRUE)
      
      min_w   <- pmax(GenomicRanges::width(data_gr),
                      as.integer(chr_lens * min_frac))
      mids    <- (GenomicRanges::start(data_gr) + GenomicRanges::end(data_gr)) / 2L
      x0_vis  <- pmax(1L, as.integer(mids - min_w / 2))
      x1_vis  <- pmin(chr_lens, as.integer(mids + min_w / 2))
      
      karyoploteR::kpBars(kp,
                          chr = chr_vec,
                          x0  = x0_vis, x1 = x1_vis,
                          y1  = vals, y0 = 0,
                          ymin = ymin, ymax = ymax,
                          r0 = r0, r1 = r1,
                          col = trk$col, border = trk$border %||% NA)
      
    } else if (trk$type == "rect") {
      n_reg <- length(data_gr)
      karyoploteR::kpRect(kp,
                          data    = data_gr,
                          y0      = rep(trk$ybottom %||% ymin, n_reg),
                          y1      = rep(trk$ytop    %||% ymax, n_reg),
                          ymin    = ymin, ymax = ymax,
                          r0 = r0, r1 = r1,
                          col = trk$col, border = trk$border %||% NA)
      
    } else {
      warning(sprintf("Track %d: unknown type '%s'. Skipping.", i, trk$type))
    }
    
    # -- track name annotation (right margin) --
    # Use kp$chromosomes (chromosomes actually rendered, respects zoom) rather
    # than the pre-zoom `chromosomes` vector, which may include chromosomes
    # outside the zoomed region and would produce zero-length coordinate vectors.
    trk_name    <- trk$name %||% paste0("Track ", i)
    label_chr   <- kp$chromosomes[length(kp$chromosomes)]
    tryCatch(
      karyoploteR::kpText(kp,
                          chr     = label_chr,
                          x       = Inf, y = (ymin + ymax) / 2,
                          labels  = trk_name,
                          ymin    = ymin, ymax = ymax,
                          r0 = r0, r1 = r1,
                          cex = axis_cex * 0.9, pos = 4, offset = 0.2,
                          clipping = FALSE),
      error = function(e) {
        warning(sprintf("Track %d: kpText label skipped (%s).", i, e$message))
      }
    )
    
    track_list[[i]] <- trk
  }
  
  # -- Base numbers ----------------------------------------------------------
  tryCatch(
    karyoploteR::kpAddBaseNumbers(kp, tick.dist = tick_dist, cex = axis_cex),
    error = function(e) warning("kpAddBaseNumbers: ", e$message)
  )
  
  # -- Legend ----------------------------------------------------------------
  if (legend.show) {
    legend_names <- sapply(seq_along(track_list), function(i)
      track_list[[i]]$name %||% paste0("Track ", i))
    legend_cols <- sapply(track_list, function(t) t$col)
    
    graphics::legend(
      legend.position,
      legend = legend_names,
      fill   = legend_cols,
      cex    = legend.cex,
      bty    = legend.bty
    )
  }
  
  # -- Labels ----------------------------------------------------------------
  if (!is.null(label)) {
    if (is.character(label) && !inherits(label, "GRanges")) {
      # User passed a plain character string (e.g. label = "My label").
      # The `label` parameter expects a GRanges object whose metadata column
      # `label.column` holds the text to draw at each range.  Ignore silently
      # (the title_name parameter is the right place for a plain string title).
      warning("`label` must be a GRanges object, not a plain character string. ",
              "Use `title_name` for plot titles. `label` ignored.")
    } else if (!is.null(label.column)) {
      label_gr   <- if (inherits(label, "GRanges")) label else
        GenomicRanges::makeGRangesFromDataFrame(label, keep.extra.columns = TRUE)
      label_text <- GenomicRanges::mcols(label_gr)[[label.column]]
      tryCatch(
        karyoploteR::kpText(kp,
                            data   = label_gr,
                            labels = label_text,
                            cex    = label.cex,
                            r0 = 0, r1 = 1,
                            ymin = 0, ymax = 1,
                            y = 0.5),
        error = function(e) {
          warning("kpText (labels): ", e$message)
        }
      )
    }
  }
  
  # -- Title -----------------------------------------------------------------
  if (!is.null(title_name)) {
    graphics::title(main = title_name, cex.main = title.cex)
  }
  
  invisible(kp)
}
