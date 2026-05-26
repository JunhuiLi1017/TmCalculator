#' Plot circular genome tracks using circlize
#'
#' @description
#' Generate a circular genome plot with flexible track definitions using the
#' circlize package. Each track is defined in a list, allowing dynamic
#' visualization of genomic data such as GC content, Tm, and peak regions.
#'
#' @param genome_name Character. Genome name.
#' @param genome_size Numeric. Genome size.
#' @param track_list List. Each element defines a track with fields such as:
#'   \itemize{
#'     \item type: "line" or "rect"
#'     \item data: data.frame or GRanges
#'     \item ylim: numeric vector of length 2
#'     \item value_col: column name (for line tracks)
#'     \item col: color
#'     \item bg.border, bg.col, bg.lwd: optional track background settings
#'   }
#' @param start.degree Numeric. Starting angle of the plot.
#' @param track.height Numeric. Height of each track.
#' @param gap.after Numeric. Gap between sectors.
#' @param cell.padding Numeric vector of length 4.
#' @param track.margin Numeric vector of length 2.
#' @param circle.margin Numeric vector of length 2.
#' @param canvas.xlim Numeric vector.
#' @param canvas.ylim Numeric vector.
#' @param axis.unit Character. Axis unit.
#' @param axis.step Numeric. Axis step.
#' @param axis.cex Numeric. Axis label size.
#' @param axis.show.unit Logical. Whether to show axis unit.
#' @param legend.show Logical. Whether to show legend.
#' @param legend.position Numeric vector. Legend position.
#' @param legend.cex Numeric. Legend label size.
#' @param legend.lwd Numeric. Legend line width.
#' @param legend.seg.len Numeric. Legend segment length.
#' @param legend.box.lwd Numeric. Legend box line width.
#' @param legend.x.intersp Numeric. Legend x interspacing.
#' @param legend.bty Character. Legend box type.
#' @param legend.border Character. Legend box border color.
#' @param legend.lty Numeric. Legend line type.
#' @param label Optional GRanges/data.frame for labels.
#' @param label.column Numeric. Label column.
#' @param label.niceFacing Logical. Whether to face labels.
#' @param label.cex Numeric. Label size.
#' @param label.side Character. Label side.
#' @param label.labels_height Numeric. Label labels height.
#' @param label.connection_height Numeric. Label connection height.
#' @param label.line_lwd Numeric. Label line width.
#' @param title.cex Numeric. Title size.
#'
#' @return A circular plot.
#'
#' @importFrom circlize circos.clear circos.par circos.initializeCircularGenome
#' @importFrom circlize circos.track circos.axis circos.genomicLines
#' @importFrom circlize circos.genomicRect circos.genomicLabels
#' @importFrom graphics text
#' @importFrom grDevices colorRampPalette
#'
#' @encoding UTF-8
#' @author Junhui Li
#' @export

plot_circos_genome <- function(
    genome_name,
    genome_size,
    track_list,
    
    # layout
    start.degree = 34.47,
    track.height = 0.05,
    gap.after = 0,
    cell.padding = c(0,0,0,0),
    track.margin = c(0.005,0.005),
    circle.margin = c(0.001,0.001),
    canvas.xlim = c(-1,1),
    canvas.ylim = c(-1,1),
    
    # axis
    axis.unit = "Mb",
    axis.step = NULL,
    axis.cex = 0.6,
    axis.show.unit = FALSE,
    
    # legend
    legend.show = TRUE,
    legend.position = c(.75, 1),
    legend.cex = 0.6,
    legend.lwd = 1,
    legend.seg.len = -0.5,
    legend.box.lwd = 0.25,
    legend.x.intersp = 1,
    legend.bty = "n",
    legend.border = NA,
    legend.lty=1,
    
    # labels (1:1 with circlize)
    label = NULL,
    label.column = 4,
    label.niceFacing = TRUE,
    label.cex = 0.6,
    label.side = "inside",
    label.labels_height = 0.01,
    label.connection_height = 0.03,
    label.line_lwd = 0.25,
    
    # title
    title.cex = 0.7
) {
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # ===== internal defaults =====
  default_bg.col    <- "white"
  default_bg.border <- "black"
  default_bg.lwd    <- 0.25
  default_line.lwd  <- 0.25
  default_rect.border <- NA
  
  # ===== colors =====
  palette <- c("#67000d", "#fc9272")
  default_cols <- grDevices::colorRampPalette(palette)(length(track_list))
  
  # ===== init =====
  circlize::circos.clear()
  
  circlize::circos.par(
    start.degree = start.degree,
    track.height = track.height,
    gap.after = gap.after,
    cell.padding = cell.padding,
    track.margin = track.margin,
    circle.margin = circle.margin,
    canvas.xlim = canvas.xlim,
    canvas.ylim = canvas.ylim
  )
  
  circlize::circos.initializeCircularGenome(
    name = genome_name,
    genome_size = genome_size,
    plotType = "none"
  )
  
  # ===== axis =====
  if (is.null(axis.step)) {
    axis.step <- if (axis.unit == "Mb") {
      ifelse(genome_size < 1e7, 5e5, 1e6)
    } else if (axis.unit == "kb") {
      1e4
    } else {
      1e6
    }
  }
  
  brk <- seq(0, genome_size, by = axis.step)
  
  lab <- switch(
    axis.unit,
    "Mb" = round(brk / 1e6, 2),
    "kb" = round(brk / 1e3, 0),
    "bp" = brk
  )
  
  if (axis.show.unit) {
    lab <- paste0(lab, " ", axis.unit)
  }
  
  # ===== helper =====
  draw_track <- function(trk, idx, is_first = FALSE) {
    
    trk$type <- trk$type %||% "rect"
    
    # --- validation ---
    if (trk$type == "line") {
      if (is.null(trk$value_col)) {
        stop("For type = 'line', 'value_col' must be provided.")
      }
      if (!(trk$value_col %in% colnames(trk$data))) {
        stop(paste0("Column '", trk$value_col, "' not found in data."))
      }
    }
    
    # --- ylim ---
    if (is.null(trk$ylim)) {
      if (trk$type == "line") {
        trk$ylim <- range(trk$data[[trk$value_col]], na.rm = TRUE)
      } else {
        trk$ylim <- c(0,1)
      }
    }
    
    # --- rect band ---
    if (trk$type == "rect") {
      trk$ytop    <- trk$ytop    %||% 0.98
      trk$ybottom <- trk$ybottom %||% 0.02
    }
    
    # --- styling ---
    trk$col       <- trk$col       %||% default_cols[idx]
    trk$bg.col    <- trk$bg.col    %||% default_bg.col
    trk$bg.border <- trk$bg.border %||% default_bg.border
    trk$bg.lwd    <- trk$bg.lwd    %||% default_bg.lwd
    trk$lwd       <- trk$lwd       %||% default_line.lwd
    trk$border    <- trk$border    %||% default_rect.border
    
    # store back for legend
    track_list[[idx]] <<- trk
    
    # --- track ---
    circlize::circos.track(
      ylim = trk$ylim,
      bg.border = trk$bg.border,
      bg.col    = trk$bg.col,
      bg.lwd    = trk$bg.lwd
    )
    
    # --- axis ---
    if (is_first) {
      circlize::circos.axis(
        h = "top",
        major.at = brk,
        labels = lab,
        labels.cex = axis.cex,
        direction = "outside"
      )
    }
    
    # --- draw ---
    if (trk$type == "line") {
      
      circlize::circos.genomicLines(
        region = trk$data,
        value = trk$data[[trk$value_col]],
        col = trk$col,
        lwd = trk$lwd
      )
      
    } else {
      
      circlize::circos.genomicRect(
        trk$data,
        ytop = trk$ytop,
        ybottom = trk$ybottom,
        col = trk$col,
        border = trk$border
      )
    }
  }
  
  # ===== draw tracks =====
  for (i in seq_along(track_list)) {
    draw_track(track_list[[i]], i, is_first = (i == 1))
  }
  
  # ===== legend =====
  if (legend.show) {
    
    legend_names <- sapply(seq_along(track_list), function(i) {
      track_list[[i]]$name %||% paste0("Track", i)
    })
    
    legend_cols <- sapply(track_list, function(trk) trk$col)
    legend_fill <- sapply(track_list, function(trk) trk$bg.col)
    
    graphics::legend(
      x = legend.position[1],
      y = legend.position[2],
      legend = legend_names,
      col = legend_cols,
      fill = legend_fill,
      lty = legend.lty,
      lwd = legend.lwd,
      seg.len = legend.seg.len,
      box.lwd = legend.box.lwd,
      x.intersp = legend.x.intersp,
      cex = legend.cex,
      bty = legend.bty,
      border = legend.border
    )
  }
  
  # ===== labels =====
  if (!is.null(label)) {
    circos.genomicLabels(
      bed = label,
      labels.column = label.column,
      niceFacing = label.niceFacing,
      cex = label.cex,
      side = label.side,
      labels_height = label.labels_height,
      connection_height = label.connection_height,
      line_lwd = label.line_lwd
    )
  }
  
  # ===== title =====
  graphics::text(0, 0, genome_name, cex = title.cex)
}
