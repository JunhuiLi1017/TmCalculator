#' Plot genome tracks in linear or circular layout
#'
#' @description
#' Visualise genomic tracks as either a linear karyoploteR plot or a circular
#' plot (base R graphics).
#' Set \code{circular = TRUE} to switch to the circular layout; the same
#' \code{track_list} works for both views without modification.
#'
#' Supported track types in each list element:
#' \itemize{
#'   \item type = "rect" : data as GRanges/data.frame, drawn as rectangles
#'   \item type = "line" : data as GRanges/data.frame with value_col, drawn as lines
#'   \item type = "area" : same as line but filled (linear only)
#'   \item type = "region": piled-up regions (linear only)
#'   \item type = "coverage": coverage area plot (linear only)
#' }
#'
#' A track with \code{ideogram = TRUE} is drawn inside the chromosome bar
#' (linear) or as the outermost ring with a grey background (circular).
#' Only \code{"rect"} type is supported for ideogram tracks.
#'
#' A per-track \code{height} field (numeric, default 1) controls relative
#' sizing. In linear mode heights are proportional fractions of the data
#' panel; in circular mode they scale the radial thickness of each ring.
#'
#' A per-track \code{highlight} field (list or list-of-lists with
#' data/col/alpha/border) draws highlight bands within that track only.
#'
#' Entries with \code{type = "highlight"} draw translucent bands spanning
#' \emph{all} real tracks. These accept: data, col (default "#F1C40F"),
#' alpha (default 0.18), border (default NA), and min.degree (circular only,
#' default 0.4).
#'
#' Each list element may also contain: name, col, bg.col, border, ylim, lwd,
#' alpha (0--1 transparency), legend_font_col, ideogram, height, highlight.
#'
#' @param genome_name Character. Used for the title.
#' @param genome_size Numeric. Total genome length (single-chromosome genomes).
#' @param genome Optional GRanges, or a karyoploteR genome string (e.g.
#'   \code{"hg38"}). Ignored in circular mode.
#' @param track_list List of track specs (see above).
#' @param circular Logical. If \code{TRUE}, render as a circular plot using
#'   base R graphics instead of a karyoploteR linear plot.
#' @param label Optional data.frame/GRanges with labels.
#' @param chromosomes Character vector to restrict/reorder chromosomes (linear only).
#' @param zoom A GRanges object, a character string, or a character vector
#'   specifying one or more regions to zoom into (e.g.
#'   \code{"chr1:1e6-2e6"} or
#'   \code{c("chr1:1e6-1.2e6", "chr1:3e6-3.2e6")}).
#'   In linear mode each region is drawn as a separate stacked panel.
#'   In circular mode the regions are concatenated around the circle with
#'   small gaps between them.
#'   \code{NULL} (default) shows the full genome.
#' @param plot.type karyoploteR plot.type (linear only).
#' @param track.gap Relative gap between tracks (linear only, 0 to ~0.05).
#' @param legend.show Logical.
#' @param legend.position Legend position. Character for linear (e.g. "topright"),
#'   numeric vector of length 2 for circular (e.g. c(0.75, 1)).
#' @param legend.cex,legend.bty,legend.border See \code{\link[graphics]{legend}}.
#' @param legend.font.col Character. Default legend text colour.
#' @param legend.lwd,legend.seg.len,legend.box.lwd,legend.x.intersp,legend.lty
#'   Additional legend parameters.
#' @param title.cex Title size.
#' @param axis.cex Axis label size.
#' @param base.tick.dist Numeric or NULL (linear only).
#' @param start.degree Numeric. Starting angle in degrees for the circular
#'   layout (default 34.47, matching the circlize convention).
#' @param track.height,gap.after,cell.padding,track.margin
#'   Kept for backward compatibility; ignored in the base R circular layout.
#' @param circle.margin Numeric vector (length 2 or 4) controlling plot margins
#'   as fractions of the device size. Length 2 is recycled as
#'   \code{c(bottom, left, bottom, left)}. Default \code{c(0.001, 0.001)}.
#' @param canvas.xlim,canvas.ylim Numeric length-2 vectors controlling the
#'   plotting window (circular only). Adjust to pan or zoom into a portion
#'   of the circle.
#' @param axis.unit,axis.step,axis.show.unit Circular axis parameters.
#' @param label.column,label.niceFacing,label.cex,label.side,label.labels_height,label.connection_height,label.line_lwd
#'   Circular label parameters.
#'
#' @return Invisibly returns the KaryoPlot object (linear) or \code{NULL}
#'   (circular).
#'
#' @examples
#' \dontrun{
#' data(ecoli_rep_hotspots)
#' 
#' library("BSgenome.Ecoli.NCBI.ASM584v2")
#' genome_name <- "BSgenome.Ecoli.NCBI.ASM584v2"
#' chr_name    <- "U00096.3"
#' genome <- get(genome_name, envir = asNamespace(genome_name))
#' chr_length  <- length(genome[[chr_name]])

#' genome_name="BSgenome.Ecoli.NCBI.ASM584v2"
#' bins_gc <- make_genomiccoord(
#'   bsgenome    = genome_name,
#'   chromosomes = chr_name,
#'   window      = 200L,
#'   slide       = 200L,
#'   start       = 1,
#'   end         = chr_length,
#'   strand      = "+"
#' )
#' input_new <- list(pkg_name = genome_name, seq = bins_gc)
#' gr_batch <- to_genomic_ranges_fast(input_new)

#' tm_ASM584v2 <- tm_calculate(
#'   gr_batch,
#'   method   = "tm_nn"
#' )
#' Tm <- as.data.frame(tm_ASM584v2$gr[, c("Tm", "GC")])
#'
#' tracks <- list(
#'   list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
#'        col = "#2C3E50", bg.col = "grey", name = "MutL-AR",
#'        legend_font_col = "#2C3E50", ideogram = TRUE, height = 0.5),
#'   list(type = "line", data = Tm, value_col = "GC",
#'        name = "GC content", col = "#4A90E2",
#'        legend_font_col = "#4A90E2"),
#'   list(type = "line", data = Tm, value_col = "Tm",
#'        name = "Melting temp", col = "#E06666",
#'        legend_font_col = "#E06666", height = 2),
#'   list(type = "line", data = ecoli_rep_hotspots$bins_rep,
#'        value_col = "count", name = "Microsatellites", col = "#2ECC71",
#'        legend_font_col = "#2ECC71"),
#'   list(type = "line", data = ecoli_rep_hotspots$bins_cru,
#'        value_col = "count", name = "Cruciform", col = "#3B3E6B",
#'        legend_font_col = "#3B3E6B"),
#'   list(data = ecoli_rep_hotspots$ssdna, name = "ssDNA",
#'        col = "#8E44AD", legend_font_col = "#8E44AD"),
#'   list(type = "line", data = ecoli_rep_hotspots$bins_gatc,
#'        value_col = "count", name = "GATC sites", col = "#D35400",
#'        legend_font_col = "#D35400"),
#'   list(type = "highlight", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
#'        col = "#F1C40F", alpha = 0.18)
#' )
#'
#' # Circular
#' plot_genome_track("E. coli", genome_size = 4641652,
#'                   track_list = tracks, circular = TRUE)
#'
#' # Linear
#' plot_genome_track("E. coli", genome_size = 4641652,
#'                   track_list = tracks)
#'
#' # Linear zoom (single region)
#' plot_genome_track("E. coli", genome_size = 4641652,
#'                   track_list = tracks,
#'                   zoom = "U00096.3:1000000-2000000")
#'
#' # Linear zoom (multiple regions — stacked panels)
#' plot_genome_track("E. coli", genome_size = 4641652,
#'                   track_list = tracks,
#'                   zoom = c("U00096.3:1000000-1200000",
#'                            "U00096.3:3000000-3200000"))
#'
#' # Circular zoom (multiple regions — concatenated arcs)
#' plot_genome_track("E. coli", genome_size = 4641652,
#'                   track_list = tracks, circular = TRUE,
#'                   zoom = c("U00096.3:1000000-1200000",
#'                            "U00096.3:3000000-3200000"))
#'
#' # Circular with canvas panning
#' plot_genome_track("E. coli", genome_size = 4641652,
#'                   track_list = tracks, circular = TRUE,
#'                   canvas.xlim = c(0.5, 1), canvas.ylim = c(0, 1),
#'                   circle.margin = c(0.05, 0.05))
#' }
#'
#' @importFrom GenomicRanges GRanges seqnames start end mcols
#' @importFrom IRanges IRanges
#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom graphics legend text rect polygon lines segments par plot.new
#'   plot.window
#' @importFrom BiocGenerics start end
#' @importFrom utils tail
#' @export
plot_genome_track <- function(
    genome_name,
    genome_size = NULL,
    genome      = NULL,
    track_list,
    circular    = FALSE,

    label        = NULL,
    chromosomes  = NULL,
    zoom         = NULL,
    plot.type    = 1,
    track.gap    = 0.01,

    # legend
    legend.show       = TRUE,
    legend.position   = NULL,
    legend.cex        = 0.6,
    legend.bty        = "n",
    legend.border     = NA,
    legend.font.col   = "black",
    legend.lwd        = 1,
    legend.seg.len    = -0.5,
    legend.box.lwd    = 0.25,
    legend.x.intersp  = 1,
    legend.lty        = 1,

    # axis / title
    axis.cex       = NULL,
    title.cex      = NULL,
    base.tick.dist = NULL,

    # ---- circular-only params ----
    start.degree   = 34.47,
    track.height   = 0.05,
    gap.after      = 0,
    cell.padding   = c(0, 0, 0, 0),
    track.margin   = c(0.005, 0.005),
    circle.margin  = c(0.001, 0.001),
    canvas.xlim    = c(-1, 1),
    canvas.ylim    = c(-1, 1),
    axis.unit      = "Mb",
    axis.step      = NULL,
    axis.show.unit = FALSE,

    # ---- circular label params ----
    label.column            = 4,
    label.niceFacing        = TRUE,
    label.cex               = 0.6,
    label.side              = "inside",
    label.labels_height     = 0.01,
    label.connection_height = 0.03,
    label.line_lwd          = 0.25
) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- sensible defaults that depend on layout ----
  if (is.null(axis.cex))        axis.cex  <- if (circular) 0.6 else 0.5
  if (is.null(title.cex))       title.cex <- if (circular) 0.7 else 0.9
  if (is.null(legend.position)) legend.position <- if (circular) c(0.75, 1) else "topright"

  # ===== shared pre-processing =====
  is_hl          <- vapply(track_list, function(t) identical(t$type, "highlight"), logical(1))
  highlight_list <- track_list[is_hl]
  track_list     <- track_list[!is_hl]

  is_ideo        <- vapply(track_list, function(t) isTRUE(t$ideogram), logical(1))
  ideogram_list  <- track_list[is_ideo]
  track_list     <- track_list[!is_ideo]

  palette      <- c("#67000d", "#fc9272")
  n_all        <- length(ideogram_list) + length(track_list)
  default_cols <- grDevices::colorRampPalette(palette)(max(n_all, 1))

  # ---- resolve genome_size (needed for both layouts) ----
  if (is.null(genome_size) && !is.null(genome) && inherits(genome, "GRanges")) {
    genome_size <- sum(as.numeric(GenomicRanges::width(genome)))
  }

  # ---- parse zoom into a list of (start, end) regions ----
  zoom_regions <- NULL
  if (!is.null(zoom)) {
    parse_one_zoom <- function(z) {
      if (inherits(z, "GRanges")) {
        data.frame(start = as.numeric(GenomicRanges::start(z)),
                   end   = as.numeric(GenomicRanges::end(z)))
      } else if (is.character(z)) {
        do.call(rbind, lapply(z, function(s) {
          range_part <- sub("^[^:]+:", "", s)
          parts <- strsplit(range_part, "-")[[1]]
          data.frame(start = as.numeric(parts[1]),
                     end   = as.numeric(parts[2]))
        }))
      } else {
        stop("zoom must be a GRanges or character vector (e.g. 'chr1:1e6-2e6')")
      }
    }
    zoom_regions <- parse_one_zoom(zoom)
    # extract chromosome name from the first element (used for linear path)
    if (is.character(zoom) && grepl(":", zoom[1])) {
      zoom_chr <- sub(":.*", "", zoom[1])
    } else {
      zoom_chr <- NULL
    }
  }

  # ===================================================================
  #  CIRCULAR LAYOUT (base R graphics)
  # ===================================================================
  if (circular) {
    if (is.null(genome_size))
      stop("genome_size is required for circular layout.")

    # ---- apply zoom (multi-region) ----
    # When zoomed, data positions are remapped into a virtual coordinate
    # space where the zoom regions are concatenated end-to-end.  Small
    # gaps are drawn between regions to visually separate them.
    zoom_offsets <- NULL  # virtual start of each region
    zoom_has_gap <- NULL  # logical: does gap i (after region i) have missing data?
    zoom_gap_deg <- 2     # visual gap between regions (degrees)
    if (!is.null(zoom_regions)) {
      zr <- zoom_regions
      region_sizes <- zr$end - zr$start
      n_reg <- nrow(zr)

      # determine which inter-region gaps have missing genomic data
      zoom_has_gap <- logical(n_reg)
      for (gi in seq_len(n_reg)) {
        if (gi < n_reg) {
          zoom_has_gap[gi] <- zr$start[gi + 1L] != zr$end[gi]
        } else {
          # wrap-around: missing = (genome_size - last_end) + first_start
          zoom_has_gap[gi] <- (genome_size - zr$end[gi]) + zr$start[1L] > 0
        }
      }
      n_gaps <- sum(zoom_has_gap)

      # reserve angular space for real gaps only
      gap_total_frac <- n_gaps * zoom_gap_deg / 360
      genome_size <- sum(region_sizes) / (1 - gap_total_frac)
      gap_bp <- if (n_gaps > 0) genome_size * zoom_gap_deg / 360 else 0
      zoom_offsets <- numeric(n_reg)
      # half wrap-around gap before first region (only if wrap gap exists)
      cursor_bp <- if (zoom_has_gap[n_reg]) gap_bp / 2 else 0
      for (zi in seq_len(n_reg)) {
        zoom_offsets[zi] <- cursor_bp
        cursor_bp <- cursor_bp + region_sizes[zi]
        if (zoom_has_gap[zi]) cursor_bp <- cursor_bp + gap_bp
      }
    }

    # helper: clip + remap a data.frame into virtual coordinates
    clip_zoom <- function(df) {
      if (is.null(zoom_regions)) return(df)
      parts <- vector("list", nrow(zoom_regions))
      for (zi in seq_len(nrow(zoom_regions))) {
        zs <- zoom_regions$start[zi]; ze <- zoom_regions$end[zi]
        keep <- df$end >= zs & df$start <= ze
        if (!any(keep)) next
        sub_df <- df[keep, , drop = FALSE]
        sub_df$start <- pmax(sub_df$start, zs) - zs + zoom_offsets[zi]
        sub_df$end   <- pmin(sub_df$end,   ze) - zs + zoom_offsets[zi]
        parts[[zi]] <- sub_df
      }
      out <- do.call(rbind, parts)
      if (is.null(out)) out <- df[0L, , drop = FALSE]  # 0-row data.frame
      out
    }

    # ---- coordinate helpers ----
    # Map genomic position to angle (radians). Position 0 starts at
    # start.degree and advances clockwise (decreasing angle).
    pos2rad <- function(pos) {
      (start.degree - pos / genome_size * 360) * pi / 180
    }

    # Draw an arc-shaped polygon between two angles at given radii.
    draw_arc <- function(t1, t2, r_in, r_out,
                         col = NA, border = NA, lwd = 0.25) {
      span <- abs(t1 - t2)
      n <- max(ceiling(span * 180 / pi), 4L)
      ang <- seq(t1, t2, length.out = n)
      xs <- c(r_out * cos(ang), r_in * cos(rev(ang)))
      ys <- c(r_out * sin(ang), r_in * sin(rev(ang)))
      polygon(xs, ys, col = col, border = border, lwd = lwd)
    }

    # Draw a full ring (annulus).
    draw_ring <- function(r_in, r_out,
                          col = NA, border = "black", lwd = 0.25) {
      ang <- seq(0, 2 * pi, length.out = 361L)
      xs <- c(r_out * cos(ang), r_in * cos(rev(ang)))
      ys <- c(r_out * sin(ang), r_in * sin(rev(ang)))
      polygon(xs, ys, col = col, border = border, lwd = lwd)
    }

    # Draw ring arcs only over zoom regions (leaving gaps blank).
    draw_ring_arcs <- function(r_in, r_out,
                               col = NA, border = "black", lwd = 0.25) {
      for (zi in seq_len(nrow(zoom_regions))) {
        reg_start <- zoom_offsets[zi]
        reg_end   <- zoom_offsets[zi] +
                     (zoom_regions$end[zi] - zoom_regions$start[zi])
        draw_arc(pos2rad(reg_start), pos2rad(reg_end),
                 r_in, r_out, col = col, border = border, lwd = lwd)
      }
    }

    # Draw zigzag break marks only at gap boundaries with missing data.
    draw_break_marks <- function(r_in, r_out, col = "grey40", lwd = 1.2,
                                 n_teeth = 5, amplitude = 0.012) {
      n_reg <- nrow(zoom_regions)
      edge_positions <- numeric(0)
      for (zi in seq_len(n_reg)) {
        reg_size <- zoom_regions$end[zi] - zoom_regions$start[zi]
        # edge after this region (gap zi)
        if (zoom_has_gap[zi])
          edge_positions <- c(edge_positions, zoom_offsets[zi] + reg_size)
        # edge before this region (gap from previous region)
        prev_gi <- if (zi == 1L) n_reg else zi - 1L
        if (zoom_has_gap[prev_gi])
          edge_positions <- c(edge_positions, zoom_offsets[zi])
      }
      for (edge_bp in edge_positions) {
        th_center <- pos2rad(edge_bp)
        r_pts <- seq(r_in, r_out, length.out = n_teeth * 2 + 1)
        xs <- ys <- numeric(length(r_pts))
        for (k in seq_along(r_pts)) {
          offset <- if (k %% 2 == 0) amplitude else -amplitude
          ang <- th_center + offset
          xs[k] <- r_pts[k] * cos(ang)
          ys[k] <- r_pts[k] * sin(ang)
        }
        lines(xs, ys, col = col, lwd = lwd)
      }
    }

    # Coerce GRanges / data.frame to a plain data.frame with start/end.
    to_df <- function(d) {
      if (inherits(d, "GRanges")) d <- as.data.frame(d)
      if (!is.data.frame(d))
        stop("Unsupported track data class: ", class(d)[1])
      # coerce start/end to plain numeric (IRanges may leave S4 integer)
      if (!is.null(d$start)) d$start <- as.numeric(d$start)
      if (!is.null(d$end))   d$end   <- as.numeric(d$end)
      d
    }

    # ---- radial layout ----
    r_outer    <- 0.92
    r_inner    <- 0.28
    radial_gap <- 0.008

    all_draw <- c(ideogram_list, track_list)
    n_draw   <- length(all_draw)
    heights  <- vapply(all_draw, function(t) t$height %||% 1, numeric(1))
    total_rgap <- max(n_draw - 1, 0) * radial_gap
    avail_r    <- r_outer - r_inner - total_rgap
    norm_h     <- heights / sum(heights) * avail_r

    # Outside-in: first track (ideogram) is outermost
    ring_bounds <- vector("list", n_draw)
    cursor <- r_outer
    for (i in seq_len(n_draw)) {
      r1i <- cursor; r0i <- cursor - norm_h[i]
      ring_bounds[[i]] <- c(r0 = r0i, r1 = r1i)
      cursor <- r0i - radial_gap
    }

    # ---- init plot ----
    # circle.margin: c(bottom, left, top, right) as fraction of device,
    # converted to lines for par(mar).  Default c(0.001, 0.001) is
    # interpreted as c(bottom, left) with top = bottom, right = left.
    cm <- rep_len(circle.margin, 4L)
    dev_h <- grDevices::dev.size("in")[2]
    lpi   <- par("cin")[2]            # line height in inches
    mar_lines <- cm * dev_h / lpi     # fraction → lines
    opar <- par(mar = mar_lines, xpd = TRUE)
    on.exit(par(opar), add = TRUE)
    plot.new()
    plot.window(xlim = canvas.xlim, ylim = canvas.ylim, asp = 1)

    # ---- axis ticks & labels ----
    if (is.null(axis.step)) {
      axis.step <- if (axis.unit == "Mb") {
        ifelse(genome_size < 1e7, 5e5, 1e6)
      } else if (axis.unit == "kb") { 1e4 } else { 1e6 }
    }

    r_tick  <- r_outer + 0.015
    r_label <- r_outer + 0.045

    if (is.null(zoom_regions)) {
      # ---- full-genome axis ----
      brk <- seq(0, genome_size, by = axis.step)
      if (length(brk) > 1 && tail(brk, 1) == genome_size)
        brk <- brk[-length(brk)]
      ax_lab <- switch(axis.unit,
        "Mb" = round(brk / 1e6, 2),
        "kb" = round(brk / 1e3, 0),
        "bp" = brk
      )
      if (axis.show.unit) ax_lab <- paste0(ax_lab, " ", axis.unit)
      for (k in seq_along(brk)) {
        th <- pos2rad(brk[k])
        segments(r_outer * cos(th), r_outer * sin(th),
                 r_tick  * cos(th), r_tick  * sin(th), lwd = 0.5)
        hadj <- 0.5 - 0.5 * cos(th)
        text(r_label * cos(th), r_label * sin(th),
             labels = ax_lab[k], cex = axis.cex, adj = c(hadj, 0.5))
      }
    } else {
      # ---- per-region axis with absolute positions ----
      gap_bp <- genome_size * zoom_gap_deg / 360
      for (zi in seq_len(nrow(zoom_regions))) {
        zs <- zoom_regions$start[zi]; ze <- zoom_regions$end[zi]
        reg_size <- ze - zs
        # ticks within this region (absolute coords)
        brk_abs <- seq(ceiling(zs / axis.step) * axis.step, ze, by = axis.step)
        if (length(brk_abs) == 0) brk_abs <- zs  # at least show start
        # map to virtual coords
        brk_virt <- brk_abs - zs + zoom_offsets[zi]
        ax_lab <- switch(axis.unit,
          "Mb" = round(brk_abs / 1e6, 2),
          "kb" = round(brk_abs / 1e3, 0),
          "bp" = brk_abs
        )
        if (axis.show.unit) ax_lab <- paste0(ax_lab, " ", axis.unit)
        for (k in seq_along(brk_virt)) {
          th <- pos2rad(brk_virt[k])
          segments(r_outer * cos(th), r_outer * sin(th),
                   r_tick  * cos(th), r_tick  * sin(th), lwd = 0.5)
          hadj <- 0.5 - 0.5 * cos(th)
          text(r_label * cos(th), r_label * sin(th),
               labels = ax_lab[k], cex = axis.cex, adj = c(hadj, 0.5))
        }
        # (break marks are drawn after all tracks)
      }
    }

    # ---- draw ideogram tracks (outermost rings) ----
    n_ideo <- length(ideogram_list)
    for (j in seq_along(ideogram_list)) {
      ideo        <- ideogram_list[[j]]
      ideo_col    <- ideo$col    %||% default_cols[j]
      ideo_bg     <- ideo$bg.col %||% "grey"
      ideo_border <- ideo$border %||% NA
      ideo_nm     <- ideo$name   %||% paste0("Ideogram", j)
      if (!is.null(ideo$alpha))
        ideo_col <- grDevices::adjustcolor(ideo_col, alpha.f = ideo$alpha)
      rb          <- ring_bounds[[j]]

      if (!is.null(zoom_regions) && nrow(zoom_regions) > 1L) {
        draw_ring_arcs(rb["r0"], rb["r1"], col = ideo_bg, border = "black", lwd = 0.25)
      } else {
        draw_ring(rb["r0"], rb["r1"], col = ideo_bg, border = "black", lwd = 0.25)
      }
      df <- clip_zoom(to_df(ideo$data))
      for (ri in seq_len(nrow(df))) {
        draw_arc(pos2rad(df$start[ri]), pos2rad(df$end[ri]),
                 rb["r0"], rb["r1"], col = ideo_col, border = ideo_border)
      }
      ideogram_list[[j]]$col  <- ideo_col
      ideogram_list[[j]]$name <- ideo_nm
    }

    # ---- draw regular tracks ----
    for (i in seq_along(track_list)) {
      trk      <- track_list[[i]]
      trk$type <- trk$type %||% "rect"
      hi       <- n_ideo + i
      rb       <- ring_bounds[[hi]]

      trk$col    <- trk$col    %||% default_cols[hi]
      trk$bg.col <- trk$bg.col %||% "white"
      trk$border <- trk$border %||% NA
      trk$lwd    <- trk$lwd    %||% 0.5
      trk$name   <- trk$name   %||% paste0("Track", i)
      if (!is.null(trk$alpha))
        trk$col <- grDevices::adjustcolor(trk$col, alpha.f = trk$alpha)

      # background ring (arcs with gaps when zoomed)
      if (!is.null(zoom_regions) && nrow(zoom_regions) > 1L) {
        draw_ring_arcs(rb["r0"], rb["r1"],
                       col = trk$bg.col, border = "black", lwd = 0.25)
      } else {
        draw_ring(rb["r0"], rb["r1"],
                  col = trk$bg.col, border = "black", lwd = 0.25)
      }

      df <- clip_zoom(to_df(trk$data))

      if (trk$type == "line") {
        if (is.null(trk$value_col))
          stop("type='line' needs value_col for '", trk$name, "'")
        vals <- df[[trk$value_col]]
        if (is.null(vals))
          stop("Column '", trk$value_col, "' not found for '", trk$name, "'")
        ylim   <- trk$ylim %||% range(vals, na.rm = TRUE)
        y_norm <- (vals - ylim[1]) / diff(ylim)
        r_vals <- rb["r0"] + y_norm * (rb["r1"] - rb["r0"])
        mid_pos <- (df$start + df$end) / 2
        thetas  <- pos2rad(mid_pos)
        if (!is.null(zoom_regions) && nrow(zoom_regions) > 1L) {
          # draw lines per-region to avoid connecting across gaps
          for (zli in seq_len(nrow(zoom_regions))) {
            reg_s <- zoom_offsets[zli]
            reg_e <- zoom_offsets[zli] +
                     (zoom_regions$end[zli] - zoom_regions$start[zli])
            in_reg <- mid_pos >= reg_s & mid_pos <= reg_e
            if (sum(in_reg) < 2L) next
            sub_ord <- order(mid_pos[in_reg])
            r_sub  <- r_vals[in_reg][sub_ord]
            th_sub <- thetas[in_reg][sub_ord]
            lines(r_sub * cos(th_sub), r_sub * sin(th_sub),
                  col = trk$col, lwd = trk$lwd)
          }
        } else {
          ord <- order(mid_pos)
          lines(r_vals[ord] * cos(thetas[ord]),
                r_vals[ord] * sin(thetas[ord]),
                col = trk$col, lwd = trk$lwd)
        }
      } else {
        # rect type
        ytop <- trk$ytop    %||% 0.98
        ybot <- trk$ybottom %||% 0.02
        r_top <- rb["r0"] + ytop * (rb["r1"] - rb["r0"])
        r_bot <- rb["r0"] + ybot * (rb["r1"] - rb["r0"])
        for (ri in seq_len(nrow(df))) {
          draw_arc(pos2rad(df$start[ri]), pos2rad(df$end[ri]),
                   r_bot, r_top, col = trk$col, border = trk$border)
        }
      }

      # ---- per-track highlight ----
      if (!is.null(trk$highlight)) {
        hl_specs <- trk$highlight
        if (!is.null(hl_specs$data)) hl_specs <- list(hl_specs)
        for (hs in hl_specs) {
          hs_fill   <- grDevices::adjustcolor(hs$col %||% "#F1C40F",
                                              alpha.f = hs$alpha %||% 0.18)
          hs_border <- hs$border %||% NA
          min_deg   <- hs$min.degree %||% 0.4
          hs_df     <- clip_zoom(to_df(hs$data))
          if (nrow(hs_df) == 0L) next
          # Merge overlapping/adjacent regions to avoid sub-pixel arcs
          hs_df <- hs_df[order(hs_df$start), , drop = FALSE]
          merged_start <- hs_df$start[1]
          merged_end   <- hs_df$end[1]
          merged <- data.frame(start = numeric(0), end = numeric(0))
          for (ri in seq_len(nrow(hs_df))) {
            if (hs_df$start[ri] <= merged_end) {
              merged_end <- max(merged_end, hs_df$end[ri])
            } else {
              merged <- rbind(merged,
                              data.frame(start = merged_start, end = merged_end))
              merged_start <- hs_df$start[ri]
              merged_end   <- hs_df$end[ri]
            }
          }
          merged <- rbind(merged,
                          data.frame(start = merged_start, end = merged_end))
          for (ri in seq_len(nrow(merged))) {
            t1 <- pos2rad(merged$start[ri])
            t2 <- pos2rad(merged$end[ri])
            span_deg <- abs(t1 - t2) * 180 / pi
            if (span_deg < min_deg) {
              mid_th <- (t1 + t2) / 2
              half   <- min_deg * pi / 360
              t1 <- mid_th + half; t2 <- mid_th - half
            }
            draw_arc(t1, t2, rb["r0"], rb["r1"],
                     col = hs_fill, border = hs_border)
          }
        }
      }

      track_list[[i]] <- trk
    }

    # ---- global highlights ----
    if (length(highlight_list) > 0 && n_draw > 0) {
      hl_r_out <- ring_bounds[[1]]["r1"]
      hl_r_in  <- ring_bounds[[n_draw]]["r0"]
      for (hl in highlight_list) {
        hl_fill   <- grDevices::adjustcolor(hl$col %||% "#F1C40F",
                                            alpha.f = hl$alpha %||% 0.18)
        hl_border <- hl$border %||% NA
        min_deg   <- hl$min.degree %||% 0.4
        hl_df     <- clip_zoom(to_df(hl$data))
        if (nrow(hl_df) == 0L) next
        # Merge overlapping/adjacent regions
        hl_df <- hl_df[order(hl_df$start), , drop = FALSE]
        m_start <- hl_df$start[1]; m_end <- hl_df$end[1]
        merged <- data.frame(start = numeric(0), end = numeric(0))
        for (ri in seq_len(nrow(hl_df))) {
          if (hl_df$start[ri] <= m_end) {
            m_end <- max(m_end, hl_df$end[ri])
          } else {
            merged <- rbind(merged, data.frame(start = m_start, end = m_end))
            m_start <- hl_df$start[ri]; m_end <- hl_df$end[ri]
          }
        }
        merged <- rbind(merged, data.frame(start = m_start, end = m_end))
        for (ri in seq_len(nrow(merged))) {
          t1 <- pos2rad(merged$start[ri])
          t2 <- pos2rad(merged$end[ri])
          span_deg <- abs(t1 - t2) * 180 / pi
          if (span_deg < min_deg) {
            mid_th <- (t1 + t2) / 2
            half   <- min_deg * pi / 360
            t1 <- mid_th + half; t2 <- mid_th - half
          }
          draw_arc(t1, t2, hl_r_in, hl_r_out,
                   col = hl_fill, border = hl_border)
        }
      }
    }

    # ---- zigzag break marks at zoom gaps ----
    if (!is.null(zoom_has_gap) && any(zoom_has_gap) && n_draw > 0) {
      draw_break_marks(ring_bounds[[n_draw]]["r0"], ring_bounds[[1]]["r1"])
    }

    # ---- legend ----
    if (legend.show) {
      all_legend   <- c(ideogram_list, track_list)
      legend_names <- vapply(all_legend, function(t) t$name %||% "", character(1))
      legend_cols  <- vapply(all_legend, function(t) t$col  %||% "black", character(1))
      legend_tcols <- vapply(all_legend, function(t) {
        t$legend_font_col %||% legend.font.col
      }, character(1))
      if (is.numeric(legend.position) && length(legend.position) == 2) {
        graphics::legend(
          x = legend.position[1], y = legend.position[2],
          legend = legend_names, fill = legend_cols,
          text.col = legend_tcols, cex = legend.cex,
          bty = legend.bty, border = legend.border)
      } else {
        graphics::legend(
          legend.position, legend = legend_names, fill = legend_cols,
          text.col = legend_tcols, cex = legend.cex,
          bty = legend.bty, border = legend.border)
      }
    }

    # ---- labels ----
    if (!is.null(label)) {
      lab_df  <- clip_zoom(to_df(label))
      lab_col <- if ("label" %in% colnames(lab_df)) "label"
                 else colnames(lab_df)[label.column]
      r_conn <- r_outer + 0.02
      r_lab  <- r_outer + 0.08
      for (ri in seq_len(nrow(lab_df))) {
        mid <- (lab_df$start[ri] + lab_df$end[ri]) / 2
        th  <- pos2rad(mid)
        segments(r_outer * cos(th), r_outer * sin(th),
                 r_conn  * cos(th), r_conn  * sin(th),
                 lwd = label.line_lwd)
        hadj <- 0.5 - 0.5 * cos(th)
        text(r_lab * cos(th), r_lab * sin(th),
             labels = lab_df[[lab_col]][ri],
             cex = label.cex, adj = c(hadj, 0.5))
      }
    }

    # ---- title ----
    graphics::text(0, 0, genome_name, cex = title.cex)
    return(invisible(NULL))
  }

  # ===================================================================
  #  LINEAR LAYOUT (karyoploteR)
  # ===================================================================

  if (!requireNamespace("karyoploteR", quietly = TRUE))
    stop("Package 'karyoploteR' is required for linear genome plots. ",
         "Install it with: BiocManager::install(\"karyoploteR\")", call. = FALSE)

  # ---- resolve genome argument ----
  genome_is_string <- is.character(genome) && length(genome) == 1
  if (is.null(genome)) {
    if (is.null(genome_size))
      stop("Provide either `genome` (GRanges/character) or `genome_size` (numeric).")
    genome <- GenomicRanges::GRanges(
      seqnames = genome_name,
      ranges   = IRanges::IRanges(start = 1, end = genome_size)
    )
  }

  # ---- build zoom GRanges list (one per region) ----
  zoom_gr_list <- NULL
  if (!is.null(zoom_regions)) {
    genome_seqnames <- as.character(GenomicRanges::seqnames(genome))
    single_g <- length(unique(genome_seqnames)) == 1L
    z_chr <- if (single_g) unique(genome_seqnames) else zoom_chr
    if (is.null(z_chr))
      stop("zoom must include a chromosome name for multi-contig genomes.")
    zoom_gr_list <- lapply(seq_len(nrow(zoom_regions)), function(zi) {
      GenomicRanges::GRanges(
        seqnames = z_chr,
        ranges   = IRanges::IRanges(start = zoom_regions$start[zi],
                                    end   = zoom_regions$end[zi]))
    })
  }

  n_zoom <- if (is.null(zoom_gr_list)) 1L else length(zoom_gr_list)

  # ---- multi-region layout ----
  if (n_zoom > 1L) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::layout(matrix(seq_len(n_zoom), nrow = n_zoom))
  }

  kp_last <- NULL
  for (zi in seq_len(n_zoom)) {

    # ---- init karyoplot ----
    pp <- karyoploteR::getDefaultPlotParams(plot.type)
    if (length(ideogram_list) > 0) {
      ideo_h <- ideogram_list[[1]]$height
      if (!is.null(ideo_h)) pp$ideogramheight <- pp$ideogramheight * ideo_h
    }

    panel_title <- if (n_zoom == 1L) genome_name else {
      zr <- zoom_regions[zi, ]
      paste0(genome_name, "  [",
             format(zr$start, big.mark = ","), " - ",
             format(zr$end,   big.mark = ","), "]")
    }

    kp_args <- list(
      genome      = genome,
      chromosomes = chromosomes,
      plot.type   = plot.type,
      plot.params = pp,
      cex         = axis.cex,
      main        = panel_title,
      cex.main    = title.cex
    )
    if (!is.null(zoom_gr_list)) kp_args$zoom <- zoom_gr_list[[zi]]
    kp <- do.call(karyoploteR::plotKaryotype, kp_args)
    kp_last <- kp

    # base numbers
    if (is.null(base.tick.dist)) {
      gsize <- if (!is.null(zoom_gr_list)) {
        as.numeric(GenomicRanges::width(zoom_gr_list[[zi]]))
      } else if (genome_is_string) {
        sum(as.numeric(GenomicRanges::width(kp$genome)))
      } else {
        sum(as.numeric(GenomicRanges::width(genome)))
      }
      btd <- if (gsize < 1e7) 5e5 else 1e6
    } else {
      btd <- base.tick.dist
    }
    karyoploteR::kpAddBaseNumbers(kp, tick.dist = btd, cex = axis.cex)

    # ---- layout ----
    n_tracks <- length(track_list)
    heights  <- vapply(track_list, function(t) t$height %||% 1, numeric(1))
    total_gap <- (n_tracks - 1) * track.gap
    avail     <- 1 - total_gap
    norm_h    <- heights / sum(heights) * avail
    slot_bounds <- vector("list", n_tracks)
    cursor <- 0
    for (i in seq_len(n_tracks)) {
      r0 <- cursor; r1 <- cursor + norm_h[i]
      slot_bounds[[i]] <- c(r0 = r0, r1 = r1)
      cursor <- r1 + track.gap
    }

    # ---- helper: coerce data.frame to GRanges ----
    kp_seqnames   <- as.character(GenomicRanges::seqnames(kp$genome))
    single_contig <- length(unique(kp_seqnames)) == 1L

    to_gr <- function(d) {
      if (inherits(d, "GRanges")) {
        if (single_contig) GenomeInfoDb::seqlevels(d) <- kp_seqnames[1]
        return(d)
      }
      if (is.data.frame(d)) {
        sn <- d$seqnames %||% d$chr %||% d$chrom
        if (is.null(sn)) {
          if (single_contig) { sn <- rep(kp_seqnames[1], nrow(d))
          } else { stop("track data.frame needs a seqnames/chr column.") }
        } else if (single_contig) { sn <- rep(kp_seqnames[1], length(sn)) }
        gr <- GenomicRanges::GRanges(
          seqnames = sn,
          ranges   = IRanges::IRanges(start = d$start, end = d$end)
        )
        mcols_df <- d[, setdiff(colnames(d),
                       c("seqnames","chr","chrom","start","end","width","strand")),
                      drop = FALSE]
        if (ncol(mcols_df)) GenomicRanges::mcols(gr) <- mcols_df
        return(gr)
      }
      stop("Unsupported track data class: ", class(d)[1])
    }

    # ---- ideogram tracks ----
    n_ideo <- length(ideogram_list)
    for (j in seq_along(ideogram_list)) {
      ideo        <- ideogram_list[[j]]
      ideo_col    <- ideo$col    %||% default_cols[j]
      ideo_border <- ideo$border %||% NA
      ideo_nm     <- ideo$name   %||% paste0("Ideogram", j)
      if (!is.null(ideo$alpha))
        ideo_col <- grDevices::adjustcolor(ideo_col, alpha.f = ideo$alpha)
      karyoploteR::kpRect(
        kp, data = to_gr(ideo$data), y0 = 0, y1 = 1,
        col = ideo_col, border = ideo_border, data.panel = "ideogram"
      )
      ideogram_list[[j]]$col  <- ideo_col
      ideogram_list[[j]]$name <- ideo_nm
    }

    # ---- draw each track ----
    for (i in seq_len(n_tracks)) {
      trk <- track_list[[i]]
      trk$type <- trk$type %||% "rect"
      r0 <- slot_bounds[[i]]["r0"]; r1 <- slot_bounds[[i]]["r1"]
      col    <- trk$col    %||% default_cols[n_ideo + i]
      bg.col <- trk$bg.col %||% "white"
      border <- trk$border %||% NA
      lwd    <- trk$lwd    %||% 0.5
      nm     <- trk$name   %||% paste0("Track", i)
      if (!is.null(trk$alpha))
        col <- grDevices::adjustcolor(col, alpha.f = trk$alpha)
      gr     <- to_gr(trk$data)

      karyoploteR::kpAddLabels(kp, labels = nm, r0 = r0, r1 = r1,
                               cex = axis.cex, label.margin = 0.035)

      switch(trk$type,
        "rect" = {
          ytop <- trk$ytop %||% 1; ybottom <- trk$ybottom %||% 0
          karyoploteR::kpRect(kp, data = gr, y0 = ybottom, y1 = ytop,
                              col = col, border = border, r0 = r0, r1 = r1)
        },
        "region" = {
          karyoploteR::kpPlotRegions(kp, data = gr, col = col, border = border,
                                     r0 = r0, r1 = r1)
        },
        "coverage" = {
          karyoploteR::kpPlotCoverage(kp, data = gr, col = col, r0 = r0, r1 = r1)
        },
        "line" = {
          if (is.null(trk$value_col)) stop("type='line' needs value_col for '", nm, "'")
          vals <- GenomicRanges::mcols(gr)[[trk$value_col]]
          if (is.null(vals)) stop("Column '", trk$value_col, "' not found for '", nm, "'")
          ylim <- trk$ylim %||% range(vals, na.rm = TRUE)
          karyoploteR::kpAxis(kp, ymin = ylim[1], ymax = ylim[2],
                              r0 = r0, r1 = r1, numticks = 2, cex = axis.cex * 0.9)
          y_norm <- (vals - ylim[1]) / diff(ylim)
          karyoploteR::kpLines(kp,
            chr = as.character(GenomicRanges::seqnames(gr)),
            x   = (GenomicRanges::start(gr) + GenomicRanges::end(gr)) / 2,
            y   = y_norm, col = col, lwd = lwd, r0 = r0, r1 = r1)
        },
        "area" = {
          if (is.null(trk$value_col)) stop("type='area' needs value_col for '", nm, "'")
          vals <- GenomicRanges::mcols(gr)[[trk$value_col]]
          ylim <- trk$ylim %||% c(0, max(vals, na.rm = TRUE))
          y_norm <- (vals - ylim[1]) / diff(ylim)
          karyoploteR::kpAxis(kp, ymin = ylim[1], ymax = ylim[2],
                              r0 = r0, r1 = r1, numticks = 2, cex = axis.cex * 0.9)
          karyoploteR::kpPlotRibbon(kp, data = gr, y0 = 0, y1 = y_norm,
                                    col = col, border = border, r0 = r0, r1 = r1)
        },
        stop("Unknown track type: ", trk$type)
      )

      # per-track highlight
      if (!is.null(trk$highlight)) {
        hl_specs <- trk$highlight
        if (!is.null(hl_specs$data)) hl_specs <- list(hl_specs)
        for (hs in hl_specs) {
          hs_fill <- grDevices::adjustcolor(hs$col %||% "#F1C40F", alpha.f = hs$alpha %||% 0.18)
          karyoploteR::kpRect(kp, data = to_gr(hs$data), y0 = 0, y1 = 1,
                              col = hs_fill, border = hs$border %||% NA, r0 = r0, r1 = r1)
        }
      }

      track_list[[i]]$col  <- col
      track_list[[i]]$bg.col <- bg.col
      track_list[[i]]$name <- nm
    }

    # ---- global highlight regions ----
    if (length(highlight_list) > 0) {
      hl_r0 <- slot_bounds[[1]]["r0"]
      hl_r1 <- slot_bounds[[n_tracks]]["r1"]
      for (hl in highlight_list) {
        hl_fill <- grDevices::adjustcolor(hl$col %||% "#F1C40F", alpha.f = hl$alpha %||% 0.18)
        karyoploteR::kpRect(kp, data = to_gr(hl$data), y0 = 0, y1 = 1,
                            col = hl_fill, border = hl$border %||% NA,
                            r0 = hl_r0, r1 = hl_r1)
      }
    }

    # ---- labels ----
    if (!is.null(label)) {
      lab_gr <- to_gr(label)
      if (!"label" %in% colnames(GenomicRanges::mcols(lab_gr)))
        stop("`label` must have a 'label' column.")
      karyoploteR::kpPlotMarkers(kp, data = lab_gr, labels = lab_gr$label,
                                 text.orientation = "horizontal",
                                 r1 = 1.05, cex = axis.cex * 1.2, label.margin = 5)
    }

    # ---- legend (first panel only for multi-zoom) ----
    if (legend.show && zi == 1L) {
      all_legend   <- c(ideogram_list, track_list)
      legend_names <- vapply(all_legend, function(t) t$name %||% "", character(1))
      legend_cols  <- vapply(all_legend, function(t) t$col  %||% "black", character(1))
      legend_tcols <- vapply(all_legend, function(t) {
        t$legend_font_col %||% legend.font.col
      }, character(1))
      graphics::legend(
        legend.position,
        legend   = legend_names,
        col      = legend_cols,
        fill     = legend_cols,
        text.col = legend_tcols,
        border   = legend.border,
        bty      = legend.bty,
        cex      = legend.cex
      )
    }

  }  # end zoom panel loop

  invisible(kp_last)
}
