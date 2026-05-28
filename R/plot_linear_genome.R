#' Plot linear genome tracks using karyoploteR
#'
#' Linear equivalent of plot_circos_genome(). Takes the same track_list
#' format so the same input can produce either view.
#'
#' Supported track types in each list element:
#'   - type = "rect" : data as GRanges/data.frame, drawn as rectangles
#'   - type = "line" : data as GRanges/data.frame with value_col, drawn as lines
#'   - type = "area" : same as line but filled (ribbon from 0 to value)
#'   - type = "region": piled-up regions via kpPlotRegions (good for "ssDNA"-style tracks)
#'   - type = "coverage": coverage area plot from a set of regions
#'
#' Each list element may also contain: name, col, bg.col, border, ylim, lwd
#'
#' @param genome_name Character. Used for the title.
#' @param genome_size Numeric. Total length (single-chromosome genomes like E. coli).
#'        Alternatively pass `genome` directly (GRanges) for multi-chromosome plots.
#' @param genome Optional GRanges describing the genome. If NULL, built from
#'        genome_name + genome_size as a single contig.
#' @param track_list List of track specs (see above).
#' @param label data.frame/GRanges with seqnames/start/end/label (like circos version).
#' @param chromosomes Optional character vector to restrict/reorder chromosomes.
#' @param plot.type karyoploteR plot.type (1 = tracks above only, 2 = above+below).
#' @param track.gap Relative gap between tracks (0 to ~0.05).
#' @param legend.show Logical.
#' @param legend.position Passed to graphics::legend (e.g. "topright", or x/y).
#' @param legend.cex,legend.bty,legend.border See graphics::legend.
#' @param title.cex Title size.
#' @param axis.cex Base numbers size.
#' @param base.tick.dist Numeric or \code{NULL}. Distance between base-number
#'   ticks passed to \code{karyoploteR::kpAddBaseNumbers}. When \code{NULL},
#'   a default is chosen based on total genome size.
#'
#' @return Invisibly returns the KaryoPlot object.
#'
#' @importFrom karyoploteR plotKaryotype kpAddBaseNumbers kpRect kpLines
#'   kpPlotRegions kpPlotCoverage kpPlotRibbon kpPlotMarkers kpAxis kpAddLabels
#' @importFrom GenomicRanges GRanges seqnames start end mcols
#' @importFrom IRanges IRanges
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics legend text
#' @export
plot_linear_genome <- function(
    genome_name,
    genome_size = NULL,
    genome      = NULL,
    track_list,
    
    label        = NULL,
    chromosomes  = NULL,
    plot.type    = 1,
    track.gap    = 0.01,
    
    # legend
    legend.show       = TRUE,
    legend.position   = "topright",
    legend.cex        = 0.6,
    legend.bty        = "n",
    legend.border     = NA,
    
    # axis / title
    axis.cex   = 0.5,
    title.cex  = 0.9,
    base.tick.dist = NULL
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # ---- build genome GRanges if not provided ----
  if (is.null(genome)) {
    if (is.null(genome_size))
      stop("Provide either `genome` (GRanges) or `genome_size` (numeric).")
    genome <- GenomicRanges::GRanges(
      seqnames = genome_name,
      ranges   = IRanges::IRanges(start = 1, end = genome_size)
    )
  }
  
  # ---- default palette matches the circos version ----
  palette      <- c("#67000d", "#fc9272")
  default_cols <- grDevices::colorRampPalette(palette)(length(track_list))
  
  # ---- init karyoplot ----
  kp <- karyoploteR::plotKaryotype(
    genome      = genome,
    chromosomes = chromosomes,
    plot.type   = plot.type,
    cex         = axis.cex,
    main        = genome_name,
    cex.main    = title.cex
  )
  
  # base numbers (genomic axis)
  if (is.null(base.tick.dist)) {
    gsize <- sum(as.numeric(GenomicRanges::width(genome)))
    base.tick.dist <- if (gsize < 1e7) 5e5 else 1e6
  }
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = base.tick.dist, cex = axis.cex)
  
  # ---- layout: split data panel 1 into N equal slots with small gaps ----
  n_tracks <- length(track_list)
  slot_h   <- (1 - (n_tracks - 1) * track.gap) / n_tracks
  # draw from top (furthest from ideogram) down; reverse so track_list[[1]] is outermost,
  # matching the circos convention where the first track is the outer ring.
  slot_bounds <- lapply(seq_len(n_tracks), function(i) {
    r1 <- 1 - (i - 1) * (slot_h + track.gap)
    r0 <- r1 - slot_h
    c(r0 = r0, r1 = r1)
  })
  
  # ---- helper: coerce data.frame to GRanges if needed ----
  to_gr <- function(d) {
    if (inherits(d, "GRanges")) return(d)
    if (is.data.frame(d)) {
      sn <- d$seqnames %||% d$chr %||% d$chrom
      if (is.null(sn)) stop("track data.frame needs a seqnames/chr column.")
      gr <- GenomicRanges::GRanges(
        seqnames = sn,
        ranges   = IRanges::IRanges(start = d$start, end = d$end)
      )
      mcols_df <- d[, setdiff(colnames(d), c("seqnames","chr","chrom","start","end","width","strand")),
                    drop = FALSE]
      if (ncol(mcols_df)) GenomicRanges::mcols(gr) <- mcols_df
      return(gr)
    }
    stop("Unsupported track data class: ", class(d)[1])
  }
  
  # ---- draw each track ----
  for (i in seq_len(n_tracks)) {
    trk <- track_list[[i]]
    trk$type <- trk$type %||% "rect"
    r0 <- slot_bounds[[i]]["r0"]
    r1 <- slot_bounds[[i]]["r1"]
    
    col    <- trk$col    %||% default_cols[i]
    bg.col <- trk$bg.col %||% "white"
    border <- trk$border %||% NA
    lwd    <- trk$lwd    %||% 0.5
    nm     <- trk$name   %||% paste0("Track", i)
    
    gr <- to_gr(trk$data)
    
    # track label on the left margin
    karyoploteR::kpAddLabels(kp, labels = nm, r0 = r0, r1 = r1,
                             cex = axis.cex, label.margin = 0.035)
    
    switch(trk$type,
           "rect" = {
             ytop    <- trk$ytop    %||% 1
             ybottom <- trk$ybottom %||% 0
             karyoploteR::kpRect(kp, data = gr, y0 = ybottom, y1 = ytop,
                                 col = col, border = border,
                                 r0 = r0, r1 = r1)
           },
           "region" = {
             karyoploteR::kpPlotRegions(kp, data = gr, col = col, border = border,
                                        r0 = r0, r1 = r1)
           },
           "coverage" = {
             karyoploteR::kpPlotCoverage(kp, data = gr, col = col,
                                         r0 = r0, r1 = r1)
           },
           "line" = {
             if (is.null(trk$value_col))
               stop("type='line' needs value_col for track '", nm, "'")
             vals <- GenomicRanges::mcols(gr)[[trk$value_col]]
             if (is.null(vals))
               stop("Column '", trk$value_col, "' not found for track '", nm, "'")
             ylim <- trk$ylim %||% range(vals, na.rm = TRUE)
             karyoploteR::kpAxis(kp, ymin = ylim[1], ymax = ylim[2],
                                 r0 = r0, r1 = r1, numticks = 2, cex = axis.cex * 0.9)
             # normalize values to [0,1] within the slot
             y_norm <- (vals - ylim[1]) / diff(ylim)
             karyoploteR::kpLines(kp,
                                  chr = as.character(GenomicRanges::seqnames(gr)),
                                  x   = (GenomicRanges::start(gr) + GenomicRanges::end(gr)) / 2,
                                  y   = y_norm,
                                  col = col, lwd = lwd,
                                  r0 = r0, r1 = r1)
           },
           "area" = {
             if (is.null(trk$value_col))
               stop("type='area' needs value_col for track '", nm, "'")
             vals <- GenomicRanges::mcols(gr)[[trk$value_col]]
             ylim <- trk$ylim %||% c(0, max(vals, na.rm = TRUE))
             y_norm <- (vals - ylim[1]) / diff(ylim)
             karyoploteR::kpAxis(kp, ymin = ylim[1], ymax = ylim[2],
                                 r0 = r0, r1 = r1, numticks = 2, cex = axis.cex * 0.9)
             karyoploteR::kpPlotRibbon(kp, data = gr, y0 = 0, y1 = y_norm,
                                       col = col, border = border,
                                       r0 = r0, r1 = r1)
           },
           stop("Unknown track type: ", trk$type)
    )
    
    # stash back for legend
    track_list[[i]]$col    <- col
    track_list[[i]]$bg.col <- bg.col
    track_list[[i]]$name   <- nm
  }
  
  # ---- labels (markers) ----
  if (!is.null(label)) {
    lab_gr <- to_gr(label)
    if (!"label" %in% colnames(GenomicRanges::mcols(lab_gr)))
      stop("`label` must have a 'label' column.")
    karyoploteR::kpPlotMarkers(kp, data = lab_gr,
                               labels = lab_gr$label,
                               text.orientation = "horizontal",
                               r1 = 1.05, cex = axis.cex * 1.2,
                               label.margin = 5)
  }
  
  # ---- legend ----
  if (legend.show) {
    legend_names <- vapply(track_list, function(t) t$name %||% "", character(1))
    legend_cols  <- vapply(track_list, function(t) t$col  %||% "black", character(1))
    graphics::legend(
      legend.position,
      legend  = legend_names,
      col     = legend_cols,
      fill    = legend_cols,
      border  = legend.border,
      bty     = legend.bty,
      cex     = legend.cex
    )
  }
  
  invisible(kp)
}