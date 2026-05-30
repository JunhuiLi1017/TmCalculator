#' Plot Tm values distributed across GRanges regions on the x-axis
#'
#' @description
#' Creates a ggplot2-based linear plot where the x-axis represents individual
#' GRanges regions rather than full chromosomal coordinates. This addresses the
#' sparsity problem of whole-genome views: when probe or primer sequences are
#' scattered across a large genome, plotting against chromosomal position
#' compresses all data into a tiny fraction of the x-axis. \code{plot_tm_linear}
#' instead places each region sequentially on the x-axis so that Tm values are
#' visible at the scale of the input data.
#'
#' @param gr A \code{GRanges} object containing a \code{Tm} metadata column
#'   with numeric melting temperature values.
#' @param x_axis Character. How to represent regions on the x-axis:
#'   \itemize{
#'     \item \code{"index"} (default): Regions are numbered 1, 2, ..., N in
#'       chromosomal order (chr then start). Clean and compact.
#'     \item \code{"label"}: Each tick is labelled \code{"seqname:start-end"}.
#'       Informative but verbose; use \code{x_label_angle} to rotate labels.
#'     \item \code{"position"}: The x-axis shows genomic midpoints; a
#'       \code{facet_wrap} by chromosome is applied automatically so that
#'       each panel is zoomed to its data range.
#'   }
#' @param color_by Character. Variable mapped to point colour:
#'   \itemize{
#'     \item \code{"chromosome"} (default): Each chromosome gets a distinct
#'       colour from the chosen viridis palette.
#'     \item \code{"Tm"}: A continuous viridis colour gradient over Tm values.
#'     \item Any other metadata column name present in \code{gr} (e.g.,
#'       \code{"group"}): treated as a discrete categorical variable - each
#'       unique level gets a distinct colour from the chosen viridis palette.
#'       Useful for comparing categories such as \code{"mutation regions"} vs.
#'       \code{"non-mutation"}.
#'   }
#' @param color_palette Character. Viridis palette to use for colouring.
#'   One of \code{"viridis"} (default), \code{"magma"}, \code{"plasma"},
#'   \code{"inferno"}, or \code{"cividis"}.
#' @param facet_by_chr Logical. Split the plot into one facet per chromosome.
#'   Forced to \code{TRUE} when \code{x_axis = "position"}.
#'   Default: \code{FALSE}.
#' @param sort_by Character. Order in which regions appear on the x-axis
#'   (\code{x_axis = "index"} or \code{"label"} only):
#'   \itemize{
#'     \item \code{"position"} (default): chromosomal order (seqname then start).
#'     \item \code{"Tm"}: ascending Tm value.
#'   }
#' @param add_line Logical. Connect points with a thin line (useful for
#'   showing trends). Default: \code{FALSE}.
#' @param point_size Numeric. Size passed to \code{geom_point}. Default: \code{2}.
#' @param x_label_angle Numeric. Rotation angle for x-axis text. Most useful
#'   when \code{x_axis = "label"}. Default: \code{45}.
#' @param title_name Character. Plot title. Defaults to a generic string when
#'   \code{NULL}.
#' @param ylab Character. Y-axis label. Default: \code{"Tm (deg C)"}.
#' @param show_legend Logical. Display the colour legend. Default: \code{TRUE}.
#'
#' @return A \code{ggplot} object. Print it to render, or convert to an
#'   interactive view with \code{plotly::ggplotly()}.
#'
#' @details
#' Regions are sorted by chromosome (natural order) then start position before
#' indexing unless \code{sort_by = "Tm"}. The returned ggplot uses
#' \code{theme_bw()} and can be further customised with standard ggplot2 layers.
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#'
#' # -- Sample data -----------------------------------------------------------
#' set.seed(1)
#' gr <- GRanges(
#'   seqnames = c(rep("chr1", 40), rep("chr2", 30), rep("chr3", 20)),
#'   ranges = IRanges(
#'     start = c(sort(sample(1:249e6, 40)),
#'               sort(sample(1:243e6, 30)),
#'               sort(sample(1:198e6, 20))),
#'     width = sample(50:150, 90, replace = TRUE)
#'   ),
#'   Tm = runif(90, 55, 85)
#' )
#'
#' # -- Example 1: Default - index on x, colour by chromosome -----------------
#' plot_tm_linear(gr)
#'
#' # -- Example 2: Region labels on x-axis ------------------------------------
#' plotly::ggplotly(plot_tm_linear(gr, x_axis = "label", x_label_angle = 60))
#'
#' # -- Example 3: Colour by Tm, sorted by Tm ---------------------------------
#' plotly::ggplotly(plot_tm_linear(gr, color_by = "chromosome", sort_by = "Tm", add_line = TRUE))
#'
#' # -- Example 4: Per-chromosome positional view -----------------------------
#' plot_tm_linear(gr, x_axis = "position", color_by = "Tm",
#'                color_palette = "magma")
#'
#' # -- Example 5: Faceted by chromosome, index x-axis ------------------------
#' plot_tm_linear(gr, facet_by_chr = TRUE, color_by = "Tm",
#'                color_palette = "plasma")
#'
#' # -- Example 6: Colour by a categorical metadata column ("group") ----------
#' gr$group <- sample(c("mutation regions", "non-mutation"), 90, replace = TRUE)
#' plot_tm_linear(gr, color_by = "group")
#' plotly::ggplotly(plot_tm_linear(gr, color_by = "group", add_line = FALSE))
#'
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_color_viridis_c
#'   scale_color_manual scale_x_continuous scale_x_discrete labs theme_bw
#'   theme element_text facet_wrap
#' @importFrom GenomicRanges seqnames start end mcols
#' @importFrom viridis viridis
#' @importFrom rlang .data
#'
#' @encoding UTF-8
#' @author Junhui Li
#' @export

plot_tm_linear <- function(
    gr,
    x_axis         = c("index", "label", "position"),
    color_by       = "chromosome",
    color_palette  = c("viridis", "magma", "plasma", "inferno", "cividis"),
    facet_by_chr   = FALSE,
    sort_by        = c("position", "Tm"),
    add_line       = FALSE,
    point_size     = 2,
    x_label_angle  = 45,
    title_name     = NULL,
    ylab           = "Tm (\u00B0C)",
    show_legend    = TRUE
) {
  
  # -- Input validation -------------------------------------------------------
  if (!inherits(gr, "GRanges"))
    stop("Input 'gr' must be a GRanges object.")
  if (!"Tm" %in% names(GenomicRanges::mcols(gr)))
    stop("GRanges object must have a 'Tm' metadata column.")
  
  x_axis        <- match.arg(x_axis)
  color_palette <- match.arg(color_palette)
  sort_by       <- match.arg(sort_by)
  
  # color_by: accept "chromosome", "Tm", or any metadata column name
  builtin_color <- c("chromosome", "Tm")
  meta_cols     <- names(GenomicRanges::mcols(gr))
  if (!color_by %in% c(builtin_color, meta_cols)) {
    stop(sprintf(
      "'color_by' must be \"chromosome\", \"Tm\", or a metadata column in gr.\n",
      "  Available metadata columns: %s",
      paste(meta_cols, collapse = ", ")
    ))
  }
  # Determine whether this is a continuous (Tm) or discrete mapping
  color_continuous <- identical(color_by, "Tm")
  
  # -- Build data frame -------------------------------------------------------
  df <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(gr)),
    start      = GenomicRanges::start(gr),
    end        = GenomicRanges::end(gr),
    midpoint   = (GenomicRanges::start(gr) + GenomicRanges::end(gr)) / 2L,
    Tm         = GenomicRanges::mcols(gr)$Tm,
    stringsAsFactors = FALSE
  )
  
  # Pull in any extra metadata column used for coloring
  if (!color_by %in% c("chromosome", "Tm")) {
    df[[color_by]] <- as.character(GenomicRanges::mcols(gr)[[color_by]])
  }
  
  df$label <- paste0(df$chromosome, ":", df$start, "-", df$end)
  
  # -- Sort -------------------------------------------------------------------
  if (sort_by == "position") {
    df <- df[order(df$chromosome, df$start), ]
  } else {
    df <- df[order(df$Tm), ]
  }
  df$index <- seq_len(nrow(df))
  
  # -- Set x variable ---------------------------------------------------------
  if (x_axis == "index") {
    df$x_var      <- df$index
    x_lab         <- "Region Index"
    discrete_x    <- FALSE
    
  } else if (x_axis == "label") {
    df$x_var      <- factor(df$label, levels = df$label)
    x_lab         <- "Region"
    discrete_x    <- TRUE
    
  } else {  # position
    df$x_var      <- df$midpoint
    x_lab         <- "Genomic Position (bp)"
    discrete_x    <- FALSE
    facet_by_chr  <- TRUE          # force per-chromosome faceting
  }
  
  if (is.null(title_name))
    title_name <- "Tm Values across GRanges Regions"
  
  # -- Build aesthetic mapping ------------------------------------------------
  hover_text <- paste0(
    "Region: ", df$label,
    "\nTm: ", round(df$Tm, 2), "\u00B0C",
    "\nChr: ", df$chromosome
  )
  # Append the group column value to the hover tooltip when coloring by metadata
  if (!color_by %in% c("chromosome", "Tm")) {
    hover_text <- paste0(hover_text, "\n", color_by, ": ", df[[color_by]])
  }
  df$hover <- hover_text
  
  if (color_continuous) {
    # Continuous: color = Tm numeric gradient
    aes_map <- ggplot2::aes(
      x     = .data$x_var,
      y     = .data$Tm,
      color = .data$Tm,
      text  = .data$hover
    )
  } else if (color_by == "chromosome") {
    aes_map <- ggplot2::aes(
      x     = .data$x_var,
      y     = .data$Tm,
      color = .data$chromosome,
      text  = .data$hover
    )
  } else {
    # Discrete metadata column (e.g. "group")
    aes_map <- ggplot2::aes(
      x     = .data$x_var,
      y     = .data$Tm,
      color = .data[[color_by]],
      text  = .data$hover
    )
  }
  
  # -- Assemble plot ----------------------------------------------------------
  p <- ggplot2::ggplot(df, aes_map) +
    ggplot2::geom_point(size = point_size)
  
  if (add_line)
    p <- p + ggplot2::geom_line(linewidth = 0.35, alpha = 0.5)
  
  # -- Color scale ------------------------------------------------------------
  if (color_continuous) {
    # Continuous gradient for Tm
    p <- p + ggplot2::scale_color_viridis_c(
      option = color_palette,
      name   = "Tm (\u00B0C)"
    )
  } else if (color_by == "chromosome") {
    n_chrs     <- length(unique(df$chromosome))
    chr_colors <- stats::setNames(
      viridis::viridis(n_chrs, option = color_palette),
      sort(unique(df$chromosome))
    )
    p <- p + ggplot2::scale_color_manual(
      values = chr_colors,
      name   = "Chromosome"
    )
  } else {
    # Arbitrary discrete metadata column (e.g. "group")
    lvls       <- sort(unique(df[[color_by]]))
    grp_colors <- stats::setNames(
      viridis::viridis(length(lvls), option = color_palette),
      lvls
    )
    p <- p + ggplot2::scale_color_manual(
      values = grp_colors,
      name   = color_by
    )
  }
  
  # -- X-axis scale -----------------------------------------------------------
  if (!discrete_x) {
    p <- p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.02))
  } else {
    p <- p + ggplot2::scale_x_discrete()
  }
  
  # -- Faceting ---------------------------------------------------------------
  if (facet_by_chr)
    p <- p + ggplot2::facet_wrap(~ chromosome, scales = "free_x")
  
  # -- Theme and labels -------------------------------------------------------
  legend_pos <- if (show_legend) "right" else "none"
  
  p <- p +
    ggplot2::labs(title = title_name, x = x_lab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(
        angle = x_label_angle,
        hjust = if (x_label_angle == 0) 0.5 else 1,
        size  = 7),
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = legend_pos,
      strip.text      = ggplot2::element_text(face = "bold")
    )
  
  return(p)
}
