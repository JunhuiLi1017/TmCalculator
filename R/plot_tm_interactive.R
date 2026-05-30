#' Convert Tm plots to interactive plotly versions
#'
#' These functions convert the standard Tm plots to interactive plotly versions
#' that can be used in Shiny applications or R Markdown documents.
#' Each function mirrors the parameters of its non-interactive counterpart and
#' adds an \code{x_axis} argument that switches between genome-coordinate and
#' region-indexed views.
#'
#' @param gr A GRanges object containing the Tm values (must have a \code{Tm}
#'   metadata column).
#' @param genome_assembly A string specifying the genome assembly (e.g., "hg19").
#' @param chromosome_to_plot A string specifying the chromosome to plot.
#' @param plot_type A string specifying the plot type (\code{"karyogram"} or
#'   \code{"faceted"}).
#' @param color_palette A string specifying the viridis color palette.
#' @param title_name A string specifying the plot title.
#' @param zoom A string specifying the zoom region (genome mode only).
#' @param x_axis Character. \code{"genome"} for chromosomal coordinate view;
#'   \code{"regions"} for region-indexed view where the x-axis shows each
#'   GRanges entry sequentially, making Tm distributions across sparse regions
#'   clearly legible. Passed to the underlying base function.
#' @param regions_facet Logical. (heatmap only) Facet the regions view by
#'   chromosome. Default: \code{FALSE}.
#' @param regions_color_by Character. (heatmap only) Colour regions by
#'   \code{"Tm"} or \code{"chromosome"}. Default: \code{"Tm"}.

#' @importFrom plotly ggplotly plot_ly layout
#' @importFrom ggplot2 ggplot aes geom_rect scale_fill_viridis_c
#'   scale_y_continuous labs theme_bw theme element_text element_blank
#'   element_line geom_point geom_segment scale_color_gradientn
#'   scale_color_manual facet_wrap
#' @importFrom viridis viridis
#' @importFrom GenomicRanges seqnames start end mcols
#' @importFrom dplyr arrange mutate group_by ungroup n row_number
#' @importFrom rlang .data

#' @rdname plot_tm_interactive
#' @export
plot_tm_heatmap_interactive <- function(gr,
                                        genome_assembly = NULL,
                                        chromosome_to_plot = NULL,
                                        plot_type = c("karyogram", "faceted"),
                                        color_palette = c("viridis", "magma", "plasma", "inferno", "cividis"),
                                        title_name = NULL,
                                        zoom = NULL,
                                        x_axis = c("genome", "regions"),
                                        regions_facet = FALSE,
                                        regions_color_by = c("Tm", "chromosome")) {

  # Delegate entirely to plot_tm_heatmap which already returns a plotly object.
  p <- plot_tm_heatmap(
    gr                = gr,
    genome_assembly   = genome_assembly,
    chromosome_to_plot = chromosome_to_plot,
    plot_type         = plot_type,
    color_palette     = color_palette,
    title_name        = title_name,
    zoom              = zoom,
    x_axis            = x_axis,
    regions_facet     = regions_facet,
    regions_color_by  = regions_color_by
  )

  return(p)
}


#' Interactive karyotype scatter plot of Tm values
#'
#' @param gr A GRanges object containing the Tm values.
#' @param chromosomes A character vector specifying chromosomes to plot.
#' @param genome_assembly A string specifying the genome assembly.
#' @param colors A named vector of colours per chromosome.
#' @param shapes A named integer vector of pch values per chromosome.
#' @param plot_type Integer. karyoploteR plot type (ignored in regions mode).
#' @param point_cex Numeric. Point size (genome mode) or point size multiplier
#'   (regions mode). Default: \code{1.5}.
#' @param xaxis_cex Numeric. X-axis label size. Default: \code{0.7}.
#' @param yaxis_cex Numeric. Y-axis label size. Default: \code{0.8}.
#' @param chr_cex Numeric. Chromosome label size. Default: \code{1}.
#' @param tick_dist Numeric. X-axis tick spacing in bp. Default: \code{10000000}.
#' @param zoom A GRanges or string specifying a zoom region (genome mode only).
#' @param x_axis Character. \code{"genome"} for chromosomal coordinates;
#'   \code{"regions"} for region-indexed view.
#'
#' @return A plotly object.
#'
#' @rdname plot_tm_interactive
#' @export
plot_tm_karyotype_interactive <- function(gr,
                                          chromosomes = NULL,
                                          genome_assembly = NULL,
                                          colors = NULL,
                                          shapes = NULL,
                                          plot_type = 1,
                                          point_cex = 1.5,
                                          xaxis_cex = 0.7,
                                          yaxis_cex = 0.8,
                                          chr_cex = 1,
                                          tick_dist = 10000000,
                                          zoom = NULL,
                                          x_axis = c("genome", "regions")) {

  x_axis <- match.arg(x_axis)

  # -- REGIONS MODE -----------------------------------------------------------
  if (x_axis == "regions") {

    chrs_to_use <- if (!is.null(chromosomes)) chromosomes else
      sort(unique(as.character(seqnames(gr))))

    plot_data <- data.frame(
      chromosome   = as.character(seqnames(gr)),
      position     = (start(gr) + end(gr)) / 2,
      Tm           = mcols(gr)$Tm,
      start        = start(gr),
      end          = end(gr),
      stringsAsFactors = FALSE
    )
    plot_data <- plot_data[plot_data$chromosome %in% chrs_to_use, ]
    plot_data <- plot_data[order(plot_data$chromosome, plot_data$position), ]
    plot_data$region_index <- seq_len(nrow(plot_data))
    plot_data$label        <- paste0(plot_data$chromosome, ":",
                                      plot_data$start, "-", plot_data$end)

    if (is.null(colors)) {
      colors <- viridis::viridis(length(unique(plot_data$chromosome)))
      names(colors) <- sort(unique(plot_data$chromosome))
    }

    p <- plotly::plot_ly() %>%
      plotly::add_trace(
        data     = plot_data,
        x        = ~region_index,
        y        = ~Tm,
        type     = "scatter",
        mode     = "markers",
        color    = ~chromosome,
        colors   = colors,
        marker   = list(size = point_cex * 6),
        text     = ~paste0(
          "Region: ", label,
          "\nTm: ",   round(Tm, 2), "\u00B0C",
          "\nChr: ",  chromosome
        ),
        hoverinfo = "text"
      ) %>%
      plotly::layout(
        title = list(text = "Tm Values by Region Index", font = list(size = 16)),
        xaxis = list(
          title    = "Region Index (ordered by chr:position)",
          tickfont = list(size = xaxis_cex * 12)
        ),
        yaxis = list(
          title    = "Tm (\u00B0C)",
          tickfont = list(size = yaxis_cex * 12)
        ),
        showlegend = TRUE,
        legend     = list(title = list(text = "Chromosome"),
                          font  = list(size = chr_cex * 12))
      )

    return(p)
  }

  # -- GENOME MODE ------------------------------------------------------------
  plot_data <- data.frame(
    chromosome = as.character(seqnames(gr)),
    position   = (start(gr) + end(gr)) / 2,
    Tm         = mcols(gr)$Tm,
    start      = start(gr),
    end        = end(gr),
    stringsAsFactors = FALSE
  )

  if (!is.null(chromosomes)) {
    plot_data <- plot_data[plot_data$chromosome %in% chromosomes, ]
  }

  if (is.null(colors)) {
    colors <- viridis::viridis(length(unique(plot_data$chromosome)))
    names(colors) <- unique(plot_data$chromosome)
  }

  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      data  = plot_data,
      x     = ~position,
      y     = ~Tm,
      type  = "scatter",
      mode  = "markers",
      color = ~chromosome,
      colors = colors,
      marker = list(
        size   = point_cex * 10,
        symbol = if (is.null(shapes)) 16 else shapes[plot_data$chromosome]
      ),
      text = ~paste(
        "Chromosome:", chromosome,
        "\nPosition:", position,
        "\nTm:", round(Tm, 2), "(\u00B0C)",
        "\nRange:", start, "-", end
      ),
      hoverinfo = "text"
    ) %>%
    plotly::layout(
      title = list(text = "Tm Values by Chromosome", font = list(size = 16)),
      xaxis = list(
        title    = "Genomic Position",
        tickfont = list(size = xaxis_cex * 12)
      ),
      yaxis = list(
        title    = "Tm (\u00B0C)",
        tickfont = list(size = yaxis_cex * 12)
      ),
      showlegend = TRUE,
      legend     = list(title = list(text = "Chromosome"),
                        font  = list(size = chr_cex * 12))
    )

  return(p)
}


#' Interactive genome-tracks plot of Tm values
#'
#' @param gr A GRanges object containing the Tm values.
#' @param chromosome_to_plot A string specifying the chromosome to plot.
#' @param genome_assembly A string specifying the genome assembly.
#' @param tm_track_title A string specifying the track title.
#' @param color_palette A string specifying the viridis palette.
#' @param show_ideogram Logical. Show ideogram (genome mode only).
#' @param zoom A string specifying the zoom region (genome mode only).
#' @param x_axis Character. \code{"genome"} for chromosomal coordinates;
#'   \code{"regions"} for region-indexed view.
#' @param regions_show_bars Logical. Add lollipop bars in regions mode.
#'   Default: \code{TRUE}.
#'
#' @return A plotly object.
#'
#' @rdname plot_tm_interactive
#' @export
plot_tm_genome_tracks_interactive <- function(gr,
                                              chromosome_to_plot,
                                              genome_assembly = NULL,
                                              tm_track_title = "Melting Temperature (\u00B0C)",
                                              color_palette = c("viridis", "magma", "plasma", "inferno", "cividis"),
                                              show_ideogram = TRUE,
                                              zoom = NULL,
                                              x_axis = c("genome", "regions"),
                                              regions_show_bars = TRUE) {

  color_palette <- match.arg(color_palette)
  x_axis        <- match.arg(x_axis)

  # -- REGIONS MODE -----------------------------------------------------------
  if (x_axis == "regions") {

    gr_chr <- gr[seqnames(gr) == chromosome_to_plot]

    plot_data <- data.frame(
      start    = start(gr_chr),
      end      = end(gr_chr),
      midpoint = (start(gr_chr) + end(gr_chr)) / 2L,
      Tm       = mcols(gr_chr)$Tm,
      stringsAsFactors = FALSE
    )
    plot_data  <- plot_data[order(plot_data$start), ]
    plot_data$region_index <- seq_len(nrow(plot_data))
    plot_data$label        <- paste0(chromosome_to_plot, ":",
                                      plot_data$start, "-", plot_data$end)

    tm_range  <- range(plot_data$Tm, na.rm = TRUE)
    pal_cols  <- viridis::viridis(256, option = color_palette)

    hover_text <- paste0(
      "Region: ", plot_data$label,
      "\nTm: ",   round(plot_data$Tm, 2), "\u00B0C"
    )

    p <- plotly::plot_ly(data = plot_data)

    if (regions_show_bars) {
      # Add vertical segments (lollipop stems)
      for (i in seq_len(nrow(plot_data))) {
        p <- p %>% plotly::add_segments(
          x    = plot_data$region_index[i],
          xend = plot_data$region_index[i],
          y    = tm_range[1],
          yend = plot_data$Tm[i],
          line = list(color = "lightgray", width = 0.8),
          showlegend = FALSE,
          hoverinfo = "none"
        )
      }
    }

    p <- p %>%
      plotly::add_trace(
        x         = ~region_index,
        y         = ~Tm,
        type      = "scatter",
        mode      = "markers",
        marker    = list(
          color     = ~Tm,
          colorscale = color_palette,
          showscale  = TRUE,
          colorbar   = list(title = tm_track_title, len = 0.8, y = 0.5),
          size       = 8
        ),
        text      = hover_text,
        hoverinfo = "text"
      ) %>%
      plotly::layout(
        title = list(
          text = paste(tm_track_title, "\u2014", chromosome_to_plot),
          font = list(size = 16)
        ),
        xaxis = list(
          title    = paste0("Region Index (", chromosome_to_plot, ", ordered by position)"),
          showgrid = TRUE
        ),
        yaxis = list(
          title    = tm_track_title,
          showgrid = TRUE
        ),
        showlegend = FALSE
      )

    return(p)
  }

  # -- GENOME MODE ------------------------------------------------------------
  gr_filtered <- gr[seqnames(gr) == chromosome_to_plot]

  plot_data <- data.frame(
    position = (start(gr_filtered) + end(gr_filtered)) / 2,
    Tm       = mcols(gr_filtered)$Tm,
    start    = start(gr_filtered),
    end      = end(gr_filtered),
    stringsAsFactors = FALSE
  )

  colors <- viridis::viridis(100, option = color_palette)

  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      data  = plot_data,
      x     = ~position,
      y     = ~Tm,
      type  = "scatter",
      mode  = "markers",
      marker = list(
        color    = ~Tm,
        colorscale = color_palette,
        showscale  = TRUE,
        colorbar   = list(title = tm_track_title, len = 0.8, y = 0.5)
      ),
      text = ~paste(
        "Position:", position,
        "\nTm:", round(Tm, 2), "(\u00B0C)",
        "\nRange:", start, "-", end
      ),
      hoverinfo = "text"
    ) %>%
    plotly::layout(
      title = list(
        text = paste(tm_track_title, "-", chromosome_to_plot),
        font = list(size = 16)
      ),
      xaxis = list(title = "Genomic Position", showgrid = TRUE),
      yaxis = list(title = tm_track_title, showgrid = TRUE),
      showlegend = FALSE
    )

  return(p)
}
