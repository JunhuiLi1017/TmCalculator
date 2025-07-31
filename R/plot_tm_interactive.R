#' Convert Tm plots to interactive plotly versions
#'
#' These functions convert the standard Tm plots to interactive plotly versions
#' that can be used in Shiny applications or R Markdown documents.
#' @param gr A GRanges object containing the Tm values.
#' @param genome_assembly A string specifying the genome assembly.
#' @param chromosome_to_plot A string specifying the chromosome to plot.
#' @param plot_type A string specifying the plot type.
#' @param color_palette A string specifying the color palette.
#' @param title_name A string specifying the title name.
#' @param zoom A string specifying the zoom level.

#' @importFrom plotly ggplotly plot_ly layout
#' @importFrom ggplot2 ggplot aes geom_rect scale_fill_viridis_c scale_y_continuous labs theme_bw theme element_text element_blank element_line
#' @importFrom viridis viridis
#' @importFrom GenomicRanges seqnames start end mcols
#' @importFrom dplyr arrange mutate group_by ungroup n
#' @importFrom rlang .data

#' @rdname plot_tm_interactive
#' @export
plot_tm_heatmap_interactive <- function(gr,
                                      genome_assembly = NULL,
                                      chromosome_to_plot = NULL,
                                      plot_type = c("karyogram", "faceted"),
                                      color_palette = c("viridis", "magma", "plasma", "inferno", "cividis"),
                                      title_name = NULL,
                                      zoom = NULL) {
  # Create the base ggplot object using the original function
  p <- plot_tm_heatmap(
    gr = gr,
    genome_assembly = genome_assembly,
    chromosome_to_plot = chromosome_to_plot,
    plot_type = plot_type,
    color_palette = color_palette,
    title_name = title_name,
    zoom = zoom
  )
  
  # Extract the ggplot object from the ggbio object
  if (inherits(p, "GGbio")) {
    # Get the ggplot object from the ggbio object
    p_ggplot <- p@ggplot
  } else {
    p_ggplot <- p
  }
  
  # Convert to plotly
  p_interactive <- plotly::ggplotly(p_ggplot) %>%
    plotly::layout(
      hovermode = "closest",
      showlegend = TRUE,
      legend = list(
        title = list(text = "Tm (\u00B0C)"),
        orientation = "v",
        y = 0.5
      )
    )
  
  return(p_interactive)
}

#' Convert Tm karyotype plots to interactive plotly versions
#'
#' These functions convert the standard Tm karyotype plots to interactive plotly versions
#' that can be used in Shiny applications or R Markdown documents.
#' @param gr A GRanges object containing the Tm values.
#' @param chromosomes A vector of strings specifying the chromosomes to plot.
#' @param genome_assembly A string specifying the genome assembly.  
#' @param colors A vector of strings specifying the colors for the chromosomes.
#' @param shapes A vector of strings specifying the shapes for the chromosomes.
#' @param plot_type A string specifying the plot type.
#' @param point_cex A numeric value specifying the point size.
#' @param xaxis_cex A numeric value specifying the x-axis label size.
#' @param yaxis_cex A numeric value specifying the y-axis label size.
#' @param chr_cex A numeric value specifying the chromosome label size.
#' @param tick_dist A numeric value specifying the tick distance.
#' @param zoom A string specifying the zoom level.
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
                                        zoom = NULL) {
  # Create a data frame for plotly
  plot_data <- data.frame(
    chromosome = as.character(seqnames(gr)),
    position = (start(gr) + end(gr)) / 2,
    Tm = mcols(gr)$Tm,
    start = start(gr),
    end = end(gr)
  )
  
  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- viridis::viridis(length(unique(plot_data$chromosome)))
    names(colors) <- unique(plot_data$chromosome)
  }
  
  # Create the interactive plot
  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      data = plot_data,
      x = ~position,
      y = ~Tm,
      type = "scatter",
      mode = "markers",
      color = ~chromosome,
      colors = colors,
      marker = list(
        size = point_cex * 10,
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
      title = list(
        text = "Tm Values by Chromosome",
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Genomic Position",
        tickfont = list(size = xaxis_cex * 12)
      ),
      yaxis = list(
        title = "Tm (\u00B0C)",
        tickfont = list(size = yaxis_cex * 12)
      ),
      showlegend = TRUE,
      legend = list(
        title = list(text = "Chromosome"),
        font = list(size = chr_cex * 12)
      )
    )
  
  return(p)
}

#' Convert Tm genome tracks plots to interactive plotly versions
#'
#' These functions convert the standard Tm genome tracks plots to interactive plotly versions
#' that can be used in Shiny applications or R Markdown documents.
#' @param gr A GRanges object containing the Tm values.
#' @param chromosome_to_plot A string specifying the chromosome to plot.
#' @param genome_assembly A string specifying the genome assembly.
#' @param tm_track_title A string specifying the title name.
#' @param color_palette A string specifying the color palette.
#' @param show_ideogram A logical value specifying whether to show the ideogram.
#' @param zoom A string specifying the zoom level.
#' 
#' @rdname plot_tm_interactive
#' @export
plot_tm_genome_tracks_interactive <- function(gr,
                                            chromosome_to_plot,
                                            genome_assembly = NULL,
                                            tm_track_title = "Melting Temperature (\u00B0C)",
                                            color_palette = c("viridis", "magma", "plasma", "inferno", "cividis"),
                                            show_ideogram = TRUE,
                                            zoom = NULL) {
  # Filter data for the specified chromosome
  gr_filtered <- gr[seqnames(gr) == chromosome_to_plot]
  
  # Create a data frame for plotly
  plot_data <- data.frame(
    position = (start(gr_filtered) + end(gr_filtered)) / 2,
    Tm = mcols(gr_filtered)$Tm,
    start = start(gr_filtered),
    end = end(gr_filtered)
  )
  
  # Get color palette
  color_palette <- match.arg(color_palette)
  colors <- viridis::viridis(100, option = color_palette)
  
  # Create the interactive plot
  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      data = plot_data,
      x = ~position,
      y = ~Tm,
      type = "scatter",
      mode = "markers",
      marker = list(
        color = ~Tm,
        colorscale = color_palette,
        showscale = TRUE,
        colorbar = list(
          title = tm_track_title,
          len = 0.8,
          y = 0.5
        )
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
      xaxis = list(
        title = "Genomic Position",
        showgrid = TRUE
      ),
      yaxis = list(
        title = tm_track_title,
        showgrid = TRUE
      ),
      showlegend = FALSE
    )
  
  return(p)
} 