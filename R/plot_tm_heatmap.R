#' Plot Tm values as a heatmap using ggplot2
#'
#' This function generates a heatmap visualization of Tm values across chromosomes
#' using ggplot2. It supports both a genome-coordinate view (faceted by chromosome,
#' x-axis = genomic position) and a region-indexed view (x-axis = region index,
#' y-axis = Tm value) that makes sparse data legible at any genomic scale.
#'
#' @param gr A GRanges object. It MUST contain a metadata column named 'Tm'
#'           with numeric melting temperature values.
#' @param genome_assembly A character string indicating the genome assembly (e.g., "hg19", "mm10").
#'                        This is used for setting chromosome lengths when needed.
#' @param chromosome_to_plot A character vector specifying which chromosomes to visualize.
#'                          These chromosomes must exist in your GRanges object.
#' @param plot_type A character string specifying the type of plot to generate:
#'                  \itemize{
#'                    \item "karyogram" (default): Chromosomes arranged as rows,
#'                      each region shown as a colored rectangle at its genomic position.
#'                    \item "faceted": Same as "karyogram" but with more whitespace
#'                      between chromosome panels.
#'                  }
#'                  Ignored when \code{x_axis = "regions"}.
#' @param color_palette A character string specifying the viridis color palette to use.
#'                     Available options are:
#'                     \itemize{
#'                       \item "viridis" (default)
#'                       \item "magma"
#'                       \item "plasma"
#'                       \item "inferno"
#'                       \item "cividis"
#'                     }
#' @param title_name A character string for the plot title.
#' @param zoom A character string specifying the genomic region to zoom into.
#'            If NULL (default), the entire range of each chromosome will be shown.
#'            Example: c("chr1:1000000-2000000", "chr2:1000000-2000000").
#'            Only used when \code{x_axis = "genome"}.
#' @param x_axis Character. Controls how the x-axis is constructed:
#'   \itemize{
#'     \item \code{"genome"} (default): x-axis represents genomic position;
#'       regions are drawn as colored rectangles in chromosomal space.
#'     \item \code{"regions"}: x-axis represents a sequential region index
#'       (ordered by chromosome then start position). y-axis shows Tm value.
#'       This view makes it easy to compare Tm values across all input regions
#'       regardless of how far apart they are in the genome.
#'   }
#' @param regions_facet Logical. When \code{x_axis = "regions"}, facet the plot
#'   by chromosome so each panel has its own x index. Default: \code{FALSE}.
#' @param regions_color_by Character. When \code{x_axis = "regions"}, map colour
#'   to \code{"Tm"} (default, continuous gradient) or \code{"chromosome"}
#'   (discrete per-chromosome colour).
#'
#' @return A plotly object (interactive).
#'
#' @importFrom ggbio layout_karyogram
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame seqinfo seqnames
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_rect geom_point geom_segment scale_fill_viridis_c
#'   scale_color_viridis_c scale_color_manual scale_y_continuous labs
#'   theme_bw theme_minimal theme element_text element_blank element_line aes
#'   facet_grid facet_wrap expansion
#' @importFrom viridis viridis
#' @importFrom dplyr arrange mutate group_by ungroup n row_number
#' @importFrom GenomeInfoDb seqlengths seqlevels seqlevelsInUse genome
#' @importFrom rlang .data
#' @importFrom plotly ggplotly layout
#'
#' @examples
#' \dontrun{
#' # Create example GRanges object
#' gr_tm <- GenomicRanges::GRanges(
#'   seqnames = c("chr1", "chr2", "chr1", "chr2", "chr1"),
#'   ranges = IRanges::IRanges(
#'     start = c(100, 200, 300, 400, 150),
#'     end = c(150, 250, 350, 450, 200)
#'   ),
#'   Tm = c(65.5, 68.2, 70.1, 63.8, 72.0)
#' )
#'
#' # Genome-coordinate view (default)
#' plot_tm_heatmap(gr_tm, genome_assembly = "hg19", plot_type = "karyogram")
#'
#' # Region-indexed view: compare Tm values across all regions
#' plot_tm_heatmap(gr_tm, genome_assembly = "hg19", x_axis = "regions")
#'
#' # Region view, faceted and coloured by chromosome
#' plot_tm_heatmap(gr_tm, genome_assembly = "hg19", x_axis = "regions",
#'                 regions_facet = TRUE, regions_color_by = "chromosome")
#' }
#'
#' @importFrom GenomeInfoDb seqinfo genome seqlengths seqlevels seqlevelsInUse
#' @importFrom GenomicRanges makeGRangesFromDataFrame mcols
#' @importFrom ggplot2 ggplot geom_rect geom_point geom_segment
#'   scale_fill_viridis_c scale_color_viridis_c scale_color_manual
#'   scale_y_continuous labs theme_bw theme_minimal theme element_text
#'   element_blank element_line aes facet_grid facet_wrap
#' @importFrom viridis viridis
#' @importFrom dplyr arrange mutate group_by ungroup n row_number %>%
#' @importFrom rlang .data
#' @export

plot_tm_heatmap <- function(gr,
                            genome_assembly = NULL,
                            chromosome_to_plot = NULL,
                            plot_type = c("karyogram","faceted"),
                            color_palette = c("viridis", "magma", "plasma", "inferno", "cividis"),
                            title_name = NULL,
                            zoom = NULL,
                            x_axis = c("genome", "regions"),
                            regions_facet = FALSE,
                            regions_color_by = c("Tm", "chromosome")) {

  # Input validation
  if (!inherits(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  if (!"Tm" %in% names(mcols(gr))) {
    stop("GRanges object must contain a 'Tm' metadata column.")
  }

  plot_type        <- match.arg(plot_type)
  color_palette    <- match.arg(color_palette)
  x_axis           <- match.arg(x_axis)
  regions_color_by <- match.arg(regions_color_by)

  if (!is.null(chromosome_to_plot)) {
    gr_filtered <- gr[seqnames(gr) == chromosome_to_plot]
    if (length(gr_filtered) == 0) {
      stop(paste0("No data points found for chromosome: ", paste(chromosome_to_plot, collapse = ", ")))
    }
    gr <- gr_filtered
  }

  # -- REGIONS MODE -----------------------------------------------------------
  if (x_axis == "regions") {

    df <- data.frame(
      chromosome = as.character(seqnames(gr)),
      start      = start(gr),
      end        = end(gr),
      Tm         = mcols(gr)$Tm,
      stringsAsFactors = FALSE
    )
    df$label <- paste0(df$chromosome, ":", df$start, "-", df$end)
    df       <- df[order(df$chromosome, df$start), ]
    df$region_index <- seq_len(nrow(df))

    tm_range <- range(df$Tm, na.rm = TRUE)

    if (is.null(title_name))
      title_name <- "Tm Values across GRanges Regions"

    hover_text <- paste0("Region: ", df$label,
                          "\nTm: ",   round(df$Tm, 2), "\u00B0C",
                          "\nChr: ",  df$chromosome)
    df$hover <- hover_text

    # Build base aesthetics
    if (regions_color_by == "Tm") {
      aes_map <- ggplot2::aes(
        x     = .data$region_index,
        y     = .data$Tm,
        color = .data$Tm,
        text  = .data$hover
      )
    } else {
      aes_map <- ggplot2::aes(
        x     = .data$region_index,
        y     = .data$Tm,
        color = .data$chromosome,
        text  = .data$hover
      )
    }

    p <- ggplot2::ggplot(df, aes_map) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = .data$region_index, yend = tm_range[1]),
        linewidth = 0.3, alpha = 0.4
      ) +
      ggplot2::geom_point(size = 2)

    # Colour scale
    if (regions_color_by == "Tm") {
      p <- p + ggplot2::scale_color_viridis_c(
        option = color_palette, name = "Tm (\u00B0C)")
    } else {
      n_chrs     <- length(unique(df$chromosome))
      chr_colors <- stats::setNames(viridis::viridis(n_chrs, option = color_palette),
                              sort(unique(df$chromosome)))
      p <- p + ggplot2::scale_color_manual(values = chr_colors, name = "Chromosome")
    }

    if (regions_facet) {
      p <- p + ggplot2::facet_wrap(~ chromosome, scales = "free_x")
    }

    p <- p +
      ggplot2::labs(
        title = title_name,
        x     = "Region Index (ordered by chr:position)",
        y     = "Tm (\u00B0C)"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        strip.text      = ggplot2::element_text(face = "bold")
      )

    p_interactive <- plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        hovermode = "closest",
        showlegend = TRUE
      )
    return(p_interactive)
  }

  # -- GENOME MODE (original behaviour) --------------------------------------

  # Validate zoom parameter if provided
  if (!is.null(zoom)) {
    zoom_list <- list()
    for (i in seq_along(zoom)) {
      if (!is.character(zoom[i]) || !grepl("^chr[0-9]+:[0-9]+-[0-9]+$", zoom[i])) {
        stop("zoom must be a character string like 'chr1:1000000-2000000'")
      }
      zoom_range <- strsplit(zoom[i], ":")[[1]]
      chr_zoom   <- zoom_range[1]
      zoom_range_pos <- as.numeric(strsplit(zoom_range[2], "-")[[1]])
      zoom_start <- zoom_range_pos[1]
      zoom_end   <- zoom_range_pos[2]
      zoom_list[[i]] <- list(chr = chr_zoom, start = zoom_start, end = zoom_end)
    }
  } else {
    zoom_list <- NULL
  }

  if (is.null(genome_assembly)) {
    if (!is.null(genome(seqinfo(gr))) && all(!is.na(genome(seqinfo(gr))))) {
      if (length(unique(genome(seqinfo(gr)))) == 1) {
        actual_genome_assembly_name <- unique(genome(seqinfo(gr)))
        message("'genome_assembly' not provided, inferred from 'gr': ", actual_genome_assembly_name)
      } else {
        stop("'genome_assembly' not provided and 'gr' contains multiple genome assemblies.")
      }
    } else {
      stop("Please provide 'genome_assembly' as a character string, or ensure it's set in your 'gr' object.")
    }
  } else if (is.character(genome_assembly)) {
    actual_genome_assembly_name <- genome_assembly
  } else if (inherits(genome_assembly, "BSgenome")) {
    actual_genome_assembly_name <- genome(genome_assembly)
    message("'genome_assembly' provided as BSgenome object, extracted genome: ", actual_genome_assembly_name)
  } else if (inherits(genome_assembly, "Seqinfo")) {
    if (!is.null(genome(genome_assembly)) && all(!is.na(genome(genome_assembly)))) {
      if (length(unique(genome(genome_assembly))) == 1) {
        actual_genome_assembly_name <- unique(genome(genome_assembly))
        message("'genome_assembly' provided as Seqinfo object, extracted genome: ", actual_genome_assembly_name)
      } else {
        stop("'genome_assembly' provided as Seqinfo object with multiple genome assemblies.")
      }
    } else {
      stop("'genome_assembly' provided as Seqinfo object, but no genome assembly name is set within it.")
    }
  } else if (inherits(genome_assembly, "GRanges")) {
    if (!is.null(genome(seqinfo(genome_assembly))) && all(!is.na(genome(seqinfo(genome_assembly))))) {
      if (length(unique(genome(seqinfo(genome_assembly)))) == 1) {
        actual_genome_assembly_name <- unique(genome(seqinfo(genome_assembly)))
        message("'genome_assembly' provided as GRanges object, extracted genome: ", actual_genome_assembly_name)
      } else {
        stop("'genome_assembly' provided as GRanges object with multiple genome assemblies.")
      }
    } else {
      stop("'genome_assembly' provided as GRanges object, but no genome assembly name is set within its seqinfo.")
    }
  } else {
    stop("Invalid 'genome_assembly' type. Must be a character string, BSgenome, Seqinfo, or GRanges object.")
  }

  # Set default chromosome lengths
  .set_default_seqlengths <- function(gr, genome_assembly = NULL) {
    if (!is.null(seqlengths(gr)) && all(!is.na(seqlengths(gr)))) {
      message("seqlengths are already defined for the GRanges object. Skipping default setting.")
      return(gr)
    }

    if (is.null(genome_assembly)) {
      stop("Please provide a 'genome_assembly' (e.g., 'hg38', 'mm10') to set default seqlengths.")
    }

    tryCatch({
      fetched_lengths <- GenomeInfoDb::seqlengths(GenomeInfoDb::Seqinfo(genome = genome_assembly))

      chr_in_gr <- unique(as.character(GenomicRanges::seqnames(gr)))
      valid_chrs_to_set <- intersect(chr_in_gr, names(fetched_lengths))

      if (length(valid_chrs_to_set) == 0) {
        warning(paste0("No matching chromosome lengths found for '", genome_assembly,
                       "' and the chromosomes in your GRanges object (",
                       paste(chr_in_gr, collapse = ", "), "). seqlengths not set."))
        return(gr)
      }

      current_seqinfo <- GenomeInfoDb::seqinfo(gr)
      new_seqinfo <- GenomeInfoDb::Seqinfo(seqnames   = names(fetched_lengths),
                                           seqlengths = fetched_lengths,
                                           genome     = genome_assembly)
      GenomeInfoDb::seqinfo(gr) <- merge(current_seqinfo, new_seqinfo)

      GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevelsInUse(gr)
      GenomeInfoDb::seqlengths(gr)[valid_chrs_to_set] <- fetched_lengths[valid_chrs_to_set]

      message("Setting default chromosome lengths for '", genome_assembly, "' for: ",
              paste(valid_chrs_to_set, collapse = ", "))

    }, error = function(e) {
      stop(paste0("Failed to fetch sequence lengths for genome assembly '", genome_assembly,
                  "'. Please check the assembly name or ensure Bioconductor data packages are installed. Error: ", e$message))
    })

    return(gr)
  }

  # Process the GRanges object
  gr <- .set_default_seqlengths(gr, genome_assembly = genome_assembly)

  # Filter / zoom
  if (!is.null(zoom_list)) {
    gr_df_list_zoomed <- list()
    for (i in seq_along(zoom_list)) {
      gr_filtered <- gr[seqnames(gr) == zoom_list[[i]]$chr &
                         start(gr) >= zoom_list[[i]]$start &
                         end(gr)   <= zoom_list[[i]]$end]
      if (length(gr_filtered) > 0) {
        gr_df <- as.data.frame(gr_filtered) %>%
          dplyr::arrange(seqnames, start) %>%
          dplyr::group_by(seqnames) %>%
          dplyr::mutate(
            y_pos    = 1:dplyr::n(),
            seq_id   = paste0("seq_", dplyr::n()),
            is_zoomed = TRUE
          ) %>%
          dplyr::ungroup()
        gr_df_list_zoomed[[i]] <- gr_df
      }
    }

    gr_df <- do.call(rbind, gr_df_list_zoomed)
    if (is.null(gr_df)) {
      stop("No data points found in the specified zoom region.")
    }
  } else {
    gr_df <- as.data.frame(gr) %>%
      dplyr::arrange(seqnames, start) %>%
      dplyr::group_by(seqnames) %>%
      dplyr::mutate(
        y_pos    = 1:dplyr::n(),
        seq_id   = paste0("seq_", dplyr::n()),
        is_zoomed = FALSE
      ) %>%
      dplyr::ungroup()
  }

  # Convert back to GRanges for ggbio, carrying over the y_pos
  gr_plot <- GenomicRanges::makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)

  # Set default title_name if not provided
  if (is.null(title_name)) {
    title_name <- if (plot_type == "karyogram") {
      "Tm Values on Chromosomes"
    } else {
      "Tm Values Faceted by Chromosome"
    }
  }

  # Create plot based on type
  if (plot_type == "karyogram") {
    p <- ggplot2::ggplot(
      gr_df,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, fill = .data$Tm)
    ) +
      ggplot2::geom_rect(
        ggplot2::aes(ymin = .data$y_pos - 0.4, ymax = .data$y_pos + 0.4),
        color = "black", linewidth = 0.1
      ) +
      ggplot2::facet_grid(seqnames ~ ., scales = "free_x", space = "free_y") +
      ggplot2::scale_y_continuous(
        breaks = gr_df[["y_pos"]],
        labels = gr_df[["seq_id"]]
      ) +
      ggplot2::scale_fill_viridis_c(option = color_palette, name = "Tm (\u00B0C)") +
      ggplot2::labs(title = title_name, x = "Genomic Position", y = "Sequence ID") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        strip.text.y       = ggplot2::element_text(angle = 0),
        axis.ticks.y       = ggplot2::element_line(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position    = "right"
      )
  } else {  # faceted
    p <- ggplot2::ggplot(gr_df,
           ggplot2::aes(xmin = .data$start, xmax = .data$end, fill = .data$Tm)) +
      ggplot2::geom_rect(
        ggplot2::aes(ymin = .data$y_pos - 0.4, ymax = .data$y_pos + 0.4),
        color = "black", linewidth = 0.1
      ) +
      ggplot2::facet_grid(seqnames ~ ., scales = "free_x", space = "free_y") +
      ggplot2::scale_y_continuous(
        breaks = gr_df[["y_pos"]],
        labels = gr_df[["seq_id"]]
      ) +
      ggplot2::scale_fill_viridis_c(option = color_palette, name = "Tm (\u00B0C)") +
      ggplot2::labs(title = title_name, x = "Genomic Position", y = "Sequence ID") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        strip.text.y       = ggplot2::element_text(angle = 0),
        axis.ticks.y       = ggplot2::element_line(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position    = "right"
      )
  }

  p_interactive <- plotly::ggplotly(p) %>%
    plotly::layout(
      hovermode  = "closest",
      showlegend = TRUE,
      legend = list(
        title       = list(text = "Tm (\u00B0C)"),
        orientation = "v",
        y           = 0.5
      )
    )
  return(p_interactive)
}
