#' Plot Tm Values on Genomic Coordinates
#'
#' This function creates a plot of Tm values across genomic coordinates using
#' ggbio and ggplot2. It can display the data either as a single integrated plot
#' with ideograms or as faceted plots by chromosome.
#'
#' @param gr A GRanges object containing the genomic coordinates and Tm values.
#'           The Tm values should be in a metadata column named "Tm".
#' @param genome_assembly Character string specifying the genome assembly
#'                        (e.g., "hg19", "hg38", "mm10") for setting chromosome lengths.
#' @param chromosome_to_plot Character string specifying the chromosome to plot.
#'                          If NULL (default), all chromosomes will be plotted.
#'                          Example: "chr1" for plotting chr1 only.
#' @param plot_type Character string specifying the plot type: "karyogram" for
#'                  integrated ideogram view or "faceted" for separate chromosome panels.
#'                  Default is "karyogram".
#' @param color_palette Character string specifying the viridis color palette to use.
#'                      Default is "viridis". Available options are:
#'                      \itemize{
#'                        \item "viridis" (default): A perceptually uniform color map that works well for most people
#'                        \item "magma": A perceptually uniform color map with a dark purple to bright yellow range
#'                        \item "plasma": A perceptually uniform color map with a dark purple to bright yellow range
#'                        \item "inferno": A perceptually uniform color map with a dark purple to bright yellow range
#'                        \item "cividis": A perceptually uniform color map optimized for color vision deficiency
#'                      }
#'                      All palettes are colorblind-friendly and perceptually uniform.
#' @param title_name Character string for the plot title. If NULL, a default title will be used.
#' @param zoom A character string vector specifying the genomic region to zoom into.
#'            If NULL (default), the entire range of each chromosome will be shown.
#'            Example: c("chr1:1000000-2000000", "chr2:1000000-2000000") for zooming into chr1:1000000-2000000 and chr2:1000000-2000000
#'
#' @return A ggplot object displaying Tm values across genomic coordinates.
#'
#' @importFrom ggbio layout_karyogram
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame seqinfo seqnames
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_rect scale_fill_viridis_c scale_y_continuous labs theme_bw theme_minimal theme element_text element_blank element_line aes
#' @importFrom viridis viridis
#' @importFrom dplyr arrange mutate group_by ungroup n
#' @importFrom GenomeInfoDb seqlengths seqlevels seqlevelsInUse genome
#' @importFrom rlang .data
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
#' # Plot with ideograms
#' plot_tm_heatmap(gr_tm, genome_assembly = "hg19", plot_type = "karyogram")
#' 
#' # Faceted plot by chromosome
#' plot_tm_heatmap(gr_tm, genome_assembly = "hg19", plot_type = "faceted")
#' 
#' # Plot with zoom
#' plot_tm_heatmap(gr_tm, genome_assembly = "hg19", plot_type = "faceted", zoom = "chr1:100-200")
#' 
#' }
#' # calculate Tm values from a fasta file
#' fasta_file <- system.file("extdata", "BSgenome.Hsapiens.UCSC.hg38.fasta", package = "TmCalculator")
#' gr_tm <- tm_calculate(fasta_file)
#' 
#' # plot with zoom
#' plot_tm_heatmap(gr_tm$tm$Tm, genome_assembly = "hg38", chromosome_to_plot = c("chr1"),
#' plot_type = "faceted", zoom = "chr1:1000-2000000")
#' 
#' @importFrom GenomeInfoDb seqinfo genome seqlengths seqlevels seqlevelsInUse
#' @importFrom GenomicRanges makeGRangesFromDataFrame mcols
#' @importFrom ggplot2 ggplot geom_rect scale_fill_viridis_c scale_y_continuous labs theme_bw theme_minimal theme element_text element_blank element_line aes
#' @importFrom viridis viridis
#' @importFrom dplyr arrange mutate group_by ungroup n %>%
#' @importFrom rlang .data
#' @export

plot_tm_heatmap <- function(gr, 
                            genome_assembly = NULL,
                            chromosome_to_plot = NULL,
                            plot_type = c("karyogram","faceted"), 
                            color_palette = c("viridis", "magma", "plasma", "inferno", "cividis"),
                            title_name = NULL,
                            zoom = NULL) {
  
  # Input validation
  if (!inherits(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }
  
  if (!"Tm" %in% names(mcols(gr))) {
    stop("GRanges object must contain a 'Tm' metadata column.")
  }
  
  plot_type <- match.arg(plot_type)
  color_palette <- match.arg(color_palette)
  
  if (!is.null(chromosome_to_plot)) {
    gr_filtered <- gr[seqnames(gr) == chromosome_to_plot]
    if (length(gr_filtered) == 0) {
      stop(paste0("No data points found for chromosome: ", chromosome_to_plot))
    }
    gr <- gr_filtered
  }
  
  # Validate zoom parameter if provided
  if (!is.null(zoom)) {
    zoom_list <- list()
    for (i in seq_along(zoom)) {
      if (!is.character(zoom) || !grepl("^chr[0-9]+:[0-9]+-[0-9]+$", zoom)) {
        stop("zoom must be a character string like 'chr1:1000000-2000000'")
      }
      zoom_range <- strsplit(zoom, ":")[[1]]
      chr_zoom <- zoom_range[1]
      zoom_range_pos <- as.numeric(strsplit(zoom_range[2], "-")[[1]])
      zoom_start <- zoom_range_pos[1]
      zoom_end <- zoom_range_pos[2]
      zoom_list[[i]] <- list(chr = chr_zoom, start = zoom_start, end = zoom_end)
    }
  } else {
    zoom_list <- NULL
  }
  
  if (is.null(genome_assembly)) {
    # If genome_assembly is NULL, try to infer from gr
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
    # Extract genome name from BSgenome object
    actual_genome_assembly_name <- genome(genome_assembly)
    message("'genome_assembly' provided as BSgenome object, extracted genome: ", actual_genome_assembly_name)
  } else if (inherits(genome_assembly, "Seqinfo")) {
    # Extract genome name from Seqinfo object
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
    # Extract genome name from GRanges object (its seqinfo)
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
      new_seqinfo <- GenomeInfoDb::Seqinfo(seqnames = names(fetched_lengths), 
                                           seqlengths = fetched_lengths, 
                                           genome = genome_assembly)
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
  
  # Get unique sequences and assign a 'y_pos' for stacking and seq_id for labeling
  # If zoom is provided, filter the GRanges object to the zoom region for each chromosome
  if (!is.null(zoom_list)) {
    # Section 1: Extract regions that match zoom list criteria
    gr_df_list_zoomed <- list()
    for (i in seq_along(zoom_list)) {
      gr_filtered <- gr[seqnames(gr) == zoom_list[[i]]$chr & 
                         start(gr) >= zoom_list[[i]]$start & 
                         end(gr) <= zoom_list[[i]]$end]
      if (length(gr_filtered) > 0) {
        gr_df <- as.data.frame(gr_filtered) %>%
          dplyr::arrange(seqnames, start) %>%
          dplyr::group_by(seqnames) %>%
          dplyr::mutate(
            y_pos = 1:dplyr::n(), # Assign y-position per chromosome
            seq_id = paste0("seq_", y_pos), # Generate seq_id based on y_pos per chromosome
            is_zoomed = TRUE # Mark as zoomed region
          ) %>%
          dplyr::ungroup()
        gr_df_list_zoomed[[i]] <- gr_df
      }
    }
    
    # Section 2: Extract regions that are not in zoom list
    zoom_chrs <- sapply(zoom_list, function(x) x$chr)
    gr_non_zoomed <- gr[!seqnames(gr) %in% zoom_chrs]
    if (length(gr_non_zoomed) > 0) {
      gr_df_non_zoomed <- as.data.frame(gr_non_zoomed) %>%
        dplyr::arrange(seqnames, start) %>%
        dplyr::group_by(seqnames) %>%
        dplyr::mutate(
          y_pos = 1:dplyr::n(), # Assign y-position per chromosome
          seq_id = paste0("seq_", y_pos), # Generate seq_id based on y_pos per chromosome
          is_zoomed = FALSE # Mark as non-zoomed region
        ) %>%
        dplyr::ungroup()
    } else {
      gr_df_non_zoomed <- NULL
    }
    
    # Combine both sections
    gr_df_list <- c(gr_df_list_zoomed, list(gr_df_non_zoomed))
    gr_df_list <- gr_df_list[!sapply(gr_df_list, is.null)] # Remove NULL entries
    gr_df <- do.call(rbind, gr_df_list)
    if (is.null(gr_df)) {
      stop("No data points found in the specified zoom region.")
    }
  } else {
    gr_df <- as.data.frame(gr) %>%
      dplyr::arrange(seqnames, start) %>%
      dplyr::group_by(seqnames) %>%
      dplyr::mutate(
        y_pos = 1:dplyr::n(), # Assign y-position per chromosome
        seq_id = paste0("seq_", y_pos), # Generate seq_id based on y_pos per chromosome
        is_zoomed = FALSE # Mark as non-zoomed region
      ) %>%
      dplyr::ungroup()
  }
  
  # Convert back to GRanges for ggbio, carrying over the y_pos
  gr_plot <- GenomicRanges::makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)
  
  # Set default title_name if not provided
  if (is.null(title_name)) {
    title_name <- if (plot_type == "karyogram") {
      "Tm Values on Chromosomes (ggbio)"
    } else {
      "Tm Values Faceted by Chromosome (ggplot2/ggbio hybrid)"
    }
  }
  
  # Create plot based on type
  if (plot_type == "karyogram") {
    p <- ggplot2::ggplot(gr_plot) +
      ggbio::layout_karyogram() +
      ggplot2::geom_rect(
        data = gr_df,
        ggplot2::aes(xmin = .data$start, xmax = .data$end, 
                     ymin = .data$y_pos - 0.4, ymax = .data$y_pos + 0.4, 
                     fill = .data$Tm),
        color = "black",
        linewidth = 0.1
      ) +
      ggplot2::scale_fill_viridis_c(option = color_palette, name = "Tm (\u00B0C)") +
      ggplot2::scale_y_continuous(
        breaks = gr_df$y_pos,
        labels = gr_df$seq_id
      ) + 
      ggplot2::labs(
        title = title_name,
        y = "Sequence ID"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
  } else { # faceted plot
    p <- ggplot2::ggplot(gr_df, ggplot2::aes(xmin = .data$start, xmax = .data$end, fill = .data$Tm)) +
      ggplot2::geom_rect(ggplot2::aes(ymin = .data$y_pos - 0.4, ymax = .data$y_pos + 0.4), 
                         color = "black", linewidth = 0.1) +
      ggplot2::facet_grid(seqnames ~ ., scales = "free_x", space = "free_y") +
      ggplot2::scale_y_continuous(
        breaks = gr_df$y_pos,
        labels = gr_df$seq_id
      ) + 
      ggplot2::scale_fill_viridis_c(option = color_palette, name = "Tm (\u00B0C)") +
      ggplot2::labs(title = title_name, x = "Genomic Position", y = "Sequence ID") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(angle = 0),
        axis.ticks.y = ggplot2::element_line(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
  }
  
  return(p)
}