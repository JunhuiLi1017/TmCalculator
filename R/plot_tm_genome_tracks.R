#' Plot Tm values as Genome Browser Tracks using Gviz
#'
#' This function generates Gviz plots displaying Tm values as DataTracks
#' alongside genome axes and ideograms for specified chromosomes. The display
#' type (gradient lines, points, connected line, or bars) is configurable, and
#' coloring can follow a continuous Tm gradient or any categorical metadata
#' column (e.g., a \code{"group"} annotation). Multiple zoom windows on the
#' same chromosome can be compared side by side.
#'
#' @param gr A GRanges object. Must contain a metadata column \code{Tm} with
#'   numeric melting temperature values.
#' @param chromosome_to_plot Character string. Chromosome to visualise; must
#'   exist in \code{gr}.
#' @param genome_assembly Character string specifying the genome assembly
#'   (e.g., \code{"hg19"}, \code{"mm10"}).
#' @param tm_track_title Character. Track / axis title.
#'   Default: \code{"Melting Temperature (\u00B0C)"}.
#' @param color_palette Character. Viridis palette: \code{"viridis"} (default),
#'   \code{"magma"}, \code{"plasma"}, \code{"inferno"}, or \code{"cividis"}.
#' @param show_ideogram Logical. Display chromosome ideogram. Default: \code{TRUE}.
#' @param zoom Character string or character vector specifying genomic region(s)
#'   to zoom into, e.g. \code{"chr1:10e6-20e6"} or
#'   \code{c(Mut = "chr1:100-1000", WT = "chr1:5e6-5.001e6")}.
#'   When the vector is \emph{named} the names become region labels.
#'   Multiple entries produce side-by-side panels in \code{"genome"} mode and
#'   filtered, labelled data in \code{"regions"} mode.
#'   \code{NULL} (default) shows the full chromosome extent.
#' @param x_axis Character. X-axis mode:
#'   \itemize{
#'     \item \code{"genome"} (default): Standard Gviz genome-browser view.
#'     \item \code{"regions"}: ggplot2 view with region index on the x-axis,
#'       suitable for sparse data spread across a large chromosome.
#'   }
#' @param track_type Character. How Tm values are drawn:
#'   \itemize{
#'     \item \code{"gradient"} (default): Vertical lines coloured by a
#'       continuous Tm gradient (lollipop view). Not compatible with discrete
#'       \code{color_by}; automatically switches to \code{"points"} when a
#'       categorical \code{color_by} is requested.
#'     \item \code{"points"}: Scatter points, coloured per-point by Tm gradient
#'       or by category.
#'     \item \code{"line"}: Connected line through all regions.
#'     \item \code{"bars"}: Bar / histogram display.
#'   }
#' @param color_by Character. Variable mapped to colour. Options:
#'   \itemize{
#'     \item \code{"Tm"} (default): continuous viridis gradient over Tm values.
#'     \item \code{"zoom_region"}: one colour per zoom window  --  requires
#'       \code{zoom} to have two or more entries.
#'     \item Any metadata column name present in \code{gr} (e.g.
#'       \code{"group"}): treated as a discrete categorical variable.
#'   }
#' @param facet_by_zoom Logical. In \code{"regions"} mode with multiple zoom
#'   windows, split into one ggplot2 facet per window. Default: \code{FALSE}.
#' @param regions_point_size Numeric. Point size in \code{"regions"} mode.
#'   Default: \code{2}.
#' @param regions_show_bars Logical. Deprecated; superseded by
#'   \code{track_type}. Retained for backward compatibility.
#'
#' @return In \code{"regions"} mode: a \code{ggplot} object.
#'   In \code{"genome"} mode: \code{invisible(NULL)} (Gviz renders directly).
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' set.seed(123)
#' gr <- GRanges(
#'   seqnames = rep("chr1", 100),
#'   ranges   = IRanges(start = sort(sample(1:249250621, 100)),
#'                      width = sample(50:200, 100, replace = TRUE)),
#'   Tm    = runif(100, 60, 80),
#'   group = sample(c("mutation", "wildtype"), 100, replace = TRUE)
#' )
#' # Standard gradient view
#' plot_tm_genome_tracks(gr, "chr1", "hg19",
#'                       zoom = "chr1:10062800-20000000")
#'
#' # Points colored by Tm gradient (fixed)
#' plot_tm_genome_tracks(gr, "chr1", "hg19",
#'                       track_type = "points",
#'                       zoom = "chr1:10062800-20000000")
#'
#' # Compare two regions side by side, colored by region
#' plot_tm_genome_tracks(gr, "chr1", "hg19",
#'   zoom      = c(Region1 = "chr1:7634457-133482943",
#'                 Region2 = "chr1:135756721-248931747"),
#'   track_type = "points",
#'   color_by   = "zoom_region")
#'
#' # Color by metadata group column
#' plot_tm_genome_tracks(gr, "chr1", "hg19",
#'   track_type = "points",
#'   color_by   = "group")
#'
#' # Regions mode - index on x-axis, multiple zoom windows
#' p <- plot_tm_genome_tracks(gr, "chr1", "hg19",
#'   x_axis   = "regions",
#'   zoom     = c(Region1 = "chr1:100-200000000",
#'                Region2 = "chr1:200000000-249250621"),
#'   color_by = "zoom_region")
#' print(p)
#' }
#'
#' @importFrom Gviz IdeogramTrack GenomeAxisTrack DataTrack OverlayTrack
#'   plotTracks
#' @importFrom GenomicRanges GRanges seqnames start end strand mcols width
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths seqlevels seqlevelsInUse genome Seqinfo
#' @importFrom viridis viridis magma plasma inferno cividis
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_col geom_segment
#'   scale_color_gradientn scale_color_manual scale_fill_gradientn
#'   scale_fill_manual facet_wrap labs theme_bw theme element_text
#' @importFrom rlang .data
#' @importFrom grid grid.newpage grid.layout viewport pushViewport popViewport
#'   unit
#' @export

plot_tm_genome_tracks <- function(
    gr,
    chromosome_to_plot,
    genome_assembly    = NULL,
    tm_track_title     = "Melting Temperature (\u00B0C)",
    color_palette      = c("viridis", "magma", "plasma", "inferno", "cividis"),
    show_ideogram      = TRUE,
    zoom               = NULL,
    x_axis             = c("genome", "regions"),
    track_type         = c("gradient", "points", "line", "bars"),
    color_by           = "Tm",
    facet_by_zoom      = FALSE,
    regions_point_size = 2,
    regions_show_bars  = TRUE   # back-compat only
) {
  
  x_axis        <- match.arg(x_axis)
  track_type    <- match.arg(track_type)
  color_palette <- match.arg(color_palette)
  
  # -- Validate gr ------------------------------------------------------------
  if (!inherits(gr, "GRanges"))
    stop("'gr' must be a GRanges object.")
  meta_cols <- names(GenomicRanges::mcols(gr))
  if (!"Tm" %in% meta_cols)
    stop("GRanges object must have a 'Tm' metadata column.")
  if (length(chromosome_to_plot) > 1) {
    warning("Only the first chromosome will be plotted.")
    chromosome_to_plot <- chromosome_to_plot[1]
  }
  if (!chromosome_to_plot %in% unique(as.character(GenomicRanges::seqnames(gr))))
    stop(sprintf("Chromosome '%s' not found in gr.", chromosome_to_plot))
  
  # -- Validate color_by ------------------------------------------------------
  valid_color_by <- c("Tm", "chromosome", "zoom_region", meta_cols)
  if (!color_by %in% valid_color_by) {
    stop(sprintf(
      paste0("'color_by' must be \"Tm\", \"chromosome\", \"zoom_region\", ",
             "or a metadata column in gr.\nAvailable metadata columns: %s"),
      paste(meta_cols, collapse = ", ")
    ))
  }
  color_continuous <- identical(color_by, "Tm")
  
  # gradient + discrete color_by -> auto-switch to points
  if (!color_continuous && identical(track_type, "gradient")) {
    message("'gradient' track_type only supports color_by = 'Tm'; switching to track_type = 'points'.")
    track_type <- "points"
  }
  
  # -- Parse zoom into a list of {chr, from, to, label} ----------------------
  zoom_list <- NULL
  
  if (!is.null(zoom)) {
    zoom_labels <- if (!is.null(names(zoom))) names(zoom) else
      paste0("Region ", seq_along(zoom))
    zoom_list <- lapply(seq_along(zoom), function(i) {
      zs    <- zoom[[i]]
      parts <- strsplit(zs, ":")[[1]]
      if (length(parts) < 2)
        stop(sprintf("zoom entry '%s': format must be 'chr:start-end'.", zs))
      pos   <- as.numeric(strsplit(parts[2], "-")[[1]])
      list(chr   = parts[1],
           from  = pos[1],
           to    = pos[2],
           label = zoom_labels[i])
    })
    bad_chrs <- sapply(zoom_list, function(z) z$chr) != chromosome_to_plot
    if (any(bad_chrs))
      stop(sprintf("zoom region(s) '%s' chromosome does not match chromosome_to_plot = '%s'.",
                   paste(zoom[bad_chrs], collapse = ", "), chromosome_to_plot))
  }
  
  # Validate zoom_region color_by
  if (identical(color_by, "zoom_region")) {
    if (is.null(zoom_list) || length(zoom_list) < 2) {
      warning("color_by = 'zoom_region' requires two or more zoom entries. Falling back to color_by = 'Tm'.")
      color_by <- "Tm"
      color_continuous <- TRUE
    }
  }
  
  # -- Color palette helpers --------------------------------------------------
  pal_256 <- viridis::viridis(256, option = color_palette)
  
  .tm_point_colors <- function(tm_vals, tm_range) {
    nm <- pmin(1, pmax(0, (tm_vals - tm_range[1]) / diff(tm_range)))
    pal_256[pmax(1L, round(nm * 255) + 1L)]
  }
  
  .group_palette <- function(group_vals) {
    lvls <- sort(unique(as.character(group_vals)))
    stats::setNames(viridis::viridis(length(lvls), option = color_palette), lvls)
  }
  
  .get_viridis_gviz_colors <- function() {
    viridis::viridis(2, option = color_palette)
  }
  
  # -- Helper: set seqlengths -------------------------------------------------
  .set_default_seqlengths <- function(gr) {
    if (!is.null(GenomeInfoDb::seqlengths(gr)) &&
        all(!is.na(GenomeInfoDb::seqlengths(gr)))) {
      message("seqlengths already defined in gr.")
      return(gr)
    }
    if (is.null(genome_assembly))
      stop("Provide 'genome_assembly' to set default seqlengths.")
    tryCatch({
      fl  <- GenomeInfoDb::seqlengths(GenomeInfoDb::Seqinfo(genome = genome_assembly))
      ok  <- intersect(unique(as.character(GenomicRanges::seqnames(gr))), names(fl))
      if (length(ok) == 0) {
        warning("No matching chromosome lengths found; seqlengths not set.")
        return(gr)
      }
      GenomeInfoDb::seqlengths(gr) <- fl[ok]
      message("Set seqlengths for: ", paste(ok, collapse = ", "))
    }, error = function(e)
      stop("Failed to fetch seqlengths for '", genome_assembly, "': ", e$message))
    gr
  }
  
  # ==========================================================================
  # REGIONS MODE (ggplot2)
  # ==========================================================================
  if (x_axis == "regions") {
    
    gr_chr <- gr[GenomicRanges::seqnames(gr) == chromosome_to_plot]
    
    plot_df <- data.frame(
      start    = GenomicRanges::start(gr_chr),
      end      = GenomicRanges::end(gr_chr),
      midpoint = (GenomicRanges::start(gr_chr) + GenomicRanges::end(gr_chr)) / 2L,
      Tm       = GenomicRanges::mcols(gr_chr)$Tm,
      chromosome = chromosome_to_plot,
      stringsAsFactors = FALSE
    )
    
    # Pull group metadata column if used
    if (!color_by %in% c("Tm", "chromosome", "zoom_region") &&
        color_by %in% meta_cols) {
      plot_df[[color_by]] <- as.character(GenomicRanges::mcols(gr_chr)[[color_by]])
    }
    
    plot_df <- plot_df[order(plot_df$start), ]
    plot_df$label <- paste0(chromosome_to_plot, ":", plot_df$start, "-", plot_df$end)
    
    # Filter by zoom windows; tag each row with its zoom region label
    if (!is.null(zoom_list)) {
      region_tag <- rep(NA_character_, nrow(plot_df))
      for (zi in seq_along(zoom_list)) {
        in_w <- plot_df$start >= zoom_list[[zi]]$from &
          plot_df$end   <= zoom_list[[zi]]$to
        region_tag[in_w] <- zoom_list[[zi]]$label
      }
      keep       <- !is.na(region_tag)
      plot_df    <- plot_df[keep, ]
      plot_df$zoom_region <- region_tag[keep]
    }
    
    if (nrow(plot_df) == 0)
      stop("No data points fall within the specified zoom region(s).")
    
    plot_df$region_index <- seq_len(nrow(plot_df))
    
    tm_range <- range(plot_df$Tm, na.rm = TRUE)
    
    # Hover tooltip
    hover_text <- paste0("Region: ", plot_df$label,
                         "\nTm: ", round(plot_df$Tm, 2), "\u00B0C")
    if (!color_continuous && color_by %in% names(plot_df))
      hover_text <- paste0(hover_text, "\n", color_by, ": ", plot_df[[color_by]])
    plot_df$hover <- hover_text
    
    # Decide the color aesthetic column
    color_col <- if (color_continuous) "Tm" else color_by
    
    # -- Build ggplot ---------------------------------------------------------
    if (track_type == "bars") {
      # bars use fill aesthetic
      fill_aes <- if (color_continuous)
        ggplot2::aes(fill = .data$Tm)
      else
        ggplot2::aes(fill = .data[[color_col]])
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$region_index, y = .data$Tm,
                     text = .data$hover)) +
        ggplot2::geom_col(fill_aes, color = NA, width = 0.8)
    } else {
      base_aes <- ggplot2::aes(
        x     = .data$region_index,
        y     = .data$Tm,
        color = .data[[color_col]],
        text  = .data$hover
      )
      p <- ggplot2::ggplot(plot_df, base_aes)
      
      if (track_type == "gradient") {
        p <- p +
          ggplot2::geom_segment(
            ggplot2::aes(xend = .data$region_index, yend = tm_range[1]),
            linewidth = 0.4, alpha = 0.6
          ) +
          ggplot2::geom_point(size = regions_point_size)
      } else if (track_type == "points") {
        p <- p + ggplot2::geom_point(size = regions_point_size)
      } else if (track_type == "line") {
        p <- p +
          ggplot2::geom_line(linewidth = 0.6, alpha = 0.8) +
          ggplot2::geom_point(size = regions_point_size * 0.6)
      }
    }
    
    # -- Color / fill scales ---------------------------------------------------
    if (color_continuous) {
      if (track_type == "bars") {
        p <- p + ggplot2::scale_fill_gradientn(
          colors = pal_256, limits = tm_range, name = tm_track_title)
      } else {
        p <- p + ggplot2::scale_color_gradientn(
          colors = pal_256, limits = tm_range, name = tm_track_title)
      }
    } else {
      lvls     <- sort(unique(plot_df[[color_col]]))
      grp_cols <- stats::setNames(viridis::viridis(length(lvls), option = color_palette), lvls)
      if (track_type == "bars") {
        p <- p + ggplot2::scale_fill_manual(values = grp_cols, name = color_by)
      } else {
        p <- p + ggplot2::scale_color_manual(values = grp_cols, name = color_by)
      }
    }
    
    # -- Optional facet by zoom region -----------------------------------------
    if (facet_by_zoom && "zoom_region" %in% names(plot_df) &&
        length(zoom_list) > 1) {
      p <- p + ggplot2::facet_wrap(~ zoom_region, scales = "free_x")
    }
    
    p <- p +
      ggplot2::labs(
        title = paste(tm_track_title, "\u2014", chromosome_to_plot),
        x     = paste0("Region index (", chromosome_to_plot, ", by position)"),
        y     = tm_track_title
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
    
    return(p)
  }
  
  # ==========================================================================
  # GENOME MODE (Gviz)
  # ==========================================================================
  gr      <- .set_default_seqlengths(gr)
  chr     <- chromosome_to_plot
  gr_chr  <- gr[GenomicRanges::seqnames(gr) == chr]
  GenomicRanges::strand(gr_chr) <- "*"
  
  tm_range <- range(GenomicRanges::mcols(gr_chr)$Tm, na.rm = TRUE)
  if (any(is.infinite(tm_range))) {
    warning("Tm range undefined; defaulting to 50-90.")
    tm_range <- c(50, 90)
  }
  
  gviz_pal  <- .get_viridis_gviz_colors()
  col_low   <- gviz_pal[1]
  col_high  <- gviz_pal[2]
  gviz_type <- switch(track_type,
                      gradient = "gradient",
                      points   = "p",
                      line     = "l",
                      bars     = "histogram"
  )
  
  # -- Compute per-feature colors --------------------------------------------
  # Returns a color vector the same length as gr_chr
  .feature_colors <- function(gr_filt) {
    tm_v <- GenomicRanges::mcols(gr_filt)$Tm
    if (color_continuous) {
      return(.tm_point_colors(tm_v, tm_range))
    }
    if (identical(color_by, "zoom_region")) {
      starts <- GenomicRanges::start(gr_filt)
      ends   <- GenomicRanges::end(gr_filt)
      tags   <- rep("other", length(gr_filt))
      for (zi in seq_along(zoom_list))
        tags[starts >= zoom_list[[zi]]$from & ends <= zoom_list[[zi]]$to] <-
        zoom_list[[zi]]$label
      # build palette from all zoom labels (consistent across panels)
      all_labels <- sapply(zoom_list, `[[`, "label")
      grp_pal    <- stats::setNames(
        viridis::viridis(length(all_labels), option = color_palette),
        all_labels
      )
      grp_pal["other"] <- "#AAAAAA"
      return(grp_pal[tags])
    }
    # metadata column
    grp_vals <- as.character(GenomicRanges::mcols(gr_filt)[[color_by]])
    # Use all levels from full gr for consistent palette
    all_grp  <- as.character(GenomicRanges::mcols(gr_chr)[[color_by]])
    .group_palette(all_grp)[grp_vals]
  }
  
  # -- Build DataTrack(s) ----------------------------------------------------
  .make_single_dt <- function(gr_filt, name) {
    pt_cols <- .feature_colors(gr_filt)
    
    args <- list(
      range            = gr_filt,
      data             = GenomicRanges::mcols(gr_filt)$Tm,
      genome           = genome_assembly,
      chromosome       = chr,
      name             = name,
      type             = gviz_type,
      cex.axis         = 0.8,
      cex.title        = 0.9,
      ylim             = tm_range,
      background.panel = "#FFFFFF",
      background.title = "gray",
      col.axis         = "black"
    )
    
    if (track_type == "gradient") {
      # gradient uses its own color scheme; ignore pt_cols
      args$gradient     <- c(col_low, col_high)
      args$color.scheme <- list(
        "gradient",
        limits = tm_range,
        color  = grDevices::colorRampPalette(c(col_low, col_high))(length(gr_filt))
      )
    } else if (track_type == "points") {
      args$pch       <- 19
      args$cex.point <- 0.8
      # Continuous Tm only; discrete color_by uses OverlayTrack (see .build_dt).
      args$col <- unname(pt_cols)
    } else if (track_type == "line") {
      # Lines get the dominant (most frequent) group color or midpoint Tm color
      args$col <- if (color_continuous) col_high else pt_cols[1]
      args$lwd <- 1.2
    } else if (track_type == "bars") {
      # Per-bar fill via col.histogram; Gviz accepts a vector here
      args$col.histogram <- pt_cols
      args$fill          <- pt_cols
    }
    
    do.call(Gviz::DataTrack, args)
  }
  
  .overlay_style <- function(col_v) {
    if (track_type == "bars") {
      list(col.histogram = col_v, fill = col_v)
    } else if (track_type == "points") {
      list(col = col_v, pch = 19, cex.point = 0.8)
    } else {
      list(col = col_v, lwd = 1.2)
    }
  }

  # Discrete color_by: one DataTrack per category -> OverlayTrack.
  # gr_filt is the already-zoom-filtered GRanges for this panel.
  .make_overlay_dt <- function(gr_filt) {
    if (identical(color_by, "zoom_region") && !is.null(zoom_list)) {
      # Each sub-track covers one zoom window within gr_filt
      all_labels <- sapply(zoom_list, `[[`, "label")
      grp_pal    <- stats::setNames(
        viridis::viridis(length(all_labels), option = color_palette),
        all_labels
      )
      dt_list <- lapply(zoom_list, function(zl) {
        gr_z <- gr_filt[GenomicRanges::start(gr_filt) >= zl$from &
                          GenomicRanges::end(gr_filt)   <= zl$to]
        if (length(gr_z) == 0) return(NULL)
        col_v   <- grp_pal[zl$label]
        col_arg <- .overlay_style(col_v)
        do.call(Gviz::DataTrack, c(list(
          range = gr_z, data = GenomicRanges::mcols(gr_z)$Tm,
          genome = genome_assembly, chromosome = chr,
          name = zl$label, type = gviz_type,
          cex.axis = 0.8, cex.title = 0.9, ylim = tm_range,
          background.panel = "#FFFFFF", background.title = "gray"
        ), col_arg))
      })
      dt_list <- Filter(Negate(is.null), dt_list)
    } else {
      # Metadata column groups  --  split gr_filt by group value
      grp_vals   <- as.character(GenomicRanges::mcols(gr_filt)[[color_by]])
      # Build palette from all levels in full gr_chr for color consistency
      all_grp    <- as.character(GenomicRanges::mcols(gr_chr)[[color_by]])
      grp_levels <- sort(unique(all_grp))
      grp_pal    <- stats::setNames(
        viridis::viridis(length(grp_levels), option = color_palette),
        grp_levels
      )
      dt_list <- lapply(grp_levels, function(g) {
        gr_g <- gr_filt[grp_vals == g]
        if (length(gr_g) == 0) return(NULL)
        col_v   <- grp_pal[g]
        col_arg <- .overlay_style(col_v)
        do.call(Gviz::DataTrack, c(list(
          range = gr_g, data = GenomicRanges::mcols(gr_g)$Tm,
          genome = genome_assembly, chromosome = chr,
          name = g, type = gviz_type,
          cex.axis = 0.8, cex.title = 0.9, ylim = tm_range,
          background.panel = "#FFFFFF", background.title = "gray"
        ), col_arg))
      })
      dt_list <- Filter(Negate(is.null), dt_list)
    }
    Gviz::OverlayTrack(trackList = dt_list, name = tm_track_title)
  }
  
  # Helper: filter gr_chr to a zoom window
  .filter_zoom <- function(zl) {
    gr_chr[GenomicRanges::start(gr_chr) >= zl$from &
             GenomicRanges::end(gr_chr)   <= zl$to]
  }
  
  # Build DataTrack(s) from zoom-filtered data only. Never pass full gr_chr
  # and rely on plotTracks(from=, to=) to subset internally  --  Gviz still
  # validates stored metadata against the post-filter row count and crashes.
  needs_overlay <- !color_continuous &&
    track_type %in% c("line", "bars", "points")
  .build_dt <- function(gr_filt) {
    if (length(gr_filt) == 0) {
      stop("No data points fall within the specified zoom region.")
    }
    if (needs_overlay) .make_overlay_dt(gr_filt) else
      .make_single_dt(gr_filt, tm_track_title)
  }
  
  # IdeogramTrack + GenomeAxisTrack (shared across panels)
  iTrack <- NULL
  if (show_ideogram && !is.null(genome_assembly))
    iTrack <- Gviz::IdeogramTrack(genome = genome_assembly, chromosome = chr)
  gTrack <- Gviz::GenomeAxisTrack()
  
  .make_tracks <- function(dt) {
    if (!is.null(iTrack)) list(iTrack, gTrack, dt) else list(gTrack, dt)
  }
  
  # -- Plot: no zoom / single zoom / multi-panel -----------------------------
  if (is.null(zoom_list)) {
    tracks <- .make_tracks(.build_dt(gr_chr))
    Gviz::plotTracks(tracks, chromosome = chr,
                     from = min(GenomicRanges::start(gr_chr)),
                     to   = max(GenomicRanges::end(gr_chr)),
                     title.width = 1.5)
    
  } else if (length(zoom_list) == 1) {
    gr_z   <- .filter_zoom(zoom_list[[1]])
    tracks <- .make_tracks(.build_dt(gr_z))
    Gviz::plotTracks(tracks, chromosome = chr,
                     from = zoom_list[[1]]$from,
                     to   = zoom_list[[1]]$to,
                     title.width = 1.5)
    
  } else {
    # Side-by-side panels: each panel gets its own DataTrack pre-filtered to
    # that zoom window, so groups / data lengths are always consistent.
    n_zoom <- length(zoom_list)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(1, n_zoom,
                                 widths = grid::unit(rep(1, n_zoom), "null"))
    ))
    for (i in seq_along(zoom_list)) {
      gr_z     <- .filter_zoom(zoom_list[[i]])
      tracks_i <- .make_tracks(.build_dt(gr_z))
      grid::pushViewport(
        grid::viewport(layout.pos.row = 1, layout.pos.col = i))
      Gviz::plotTracks(tracks_i, chromosome = chr,
                       from        = zoom_list[[i]]$from,
                       to          = zoom_list[[i]]$to,
                       title.width = 1.5,
                       main        = zoom_list[[i]]$label,
                       cex.main    = 0.85,
                       add         = TRUE)
      grid::popViewport()
    }
    grid::popViewport()
  }
  
  invisible(NULL)
}
