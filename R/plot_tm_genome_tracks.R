#' Plot Tm values as Genome Browser Tracks using Gviz
#'
#' This function generates Gviz plots displaying Tm values as DataTracks
#' alongside genome axes and ideograms for specified chromosomes.
#' Tm values are visualized using a heatmap-like color gradient.
#'
#' @param gr A GRanges object. It MUST contain a metadata column named 'Tm'
#'           with numeric melting temperature values.
#' @param chromosome_to_plot A character vector specifying the chromosomes to visualize.
#'                          These chromosomes must exist in your GRanges object.
#' @param genome_assembly A character string indicating the genome assembly (e.g., "hg19", "mm10").
#'                        This is used by IdeogramTrack for correct ideogram display.
#' @param tm_track_title A character string for the title of the Tm data track.
#' @param color_palette A character string specifying the viridis color palette to use.
#'                     Available options are:
#'                     \itemize{
#'                       \item "viridis" (default): A perceptually uniform color map that works well for most people
#'                       \item "magma": A perceptually uniform color map with a dark purple to bright yellow range
#'                       \item "plasma": A perceptually uniform color map with a dark purple to bright yellow range
#'                       \item "inferno": A perceptually uniform color map with a dark purple to bright yellow range
#'                       \item "cividis": A perceptually uniform color map optimized for color vision deficiency
#'                     }
#'                     All palettes are colorblind-friendly and perceptually uniform.
#' @param show_ideogram Logical, whether to display the chromosome ideogram tracks.
#' @param ncol Number of columns in the plot layout when multiple chromosomes are plotted.
#' @return Invisible NULL. The function generates a plot directly.
#' @examples
#' \dontrun{
#' # Example 1: Generate sample data with 150 sequences
#' set.seed(123)
#' 
#' # Generate 100 sequences for chr1
#' chr1_starts <- sort(sample(1:249250621, 100))  # chr1 length in hg19
#' chr1_lengths <- sample(50:200, 100, replace=TRUE)
#' chr1_ends <- chr1_starts + chr1_lengths
#' chr1_tms <- runif(100, min=60, max=80)
#' 
#' # Generate 50 sequences for chr2
#' chr2_starts <- sort(sample(1:243199373, 50))   # chr2 length in hg19
#' chr2_lengths <- sample(50:200, 50, replace=TRUE)
#' chr2_ends <- chr2_starts + chr2_lengths
#' chr2_tms <- runif(50, min=60, max=80)
#' 
#' # Create GRanges object
#' tm_results <- GRanges(
#'   seqnames = Rle(c(rep("chr1", 100), rep("chr2", 50))),
#'   ranges = IRanges(
#'     start = c(chr1_starts, chr2_starts),
#'     end = c(chr1_ends, chr2_ends)
#'   ),
#'   strand = Rle(sample(c("+", "-"), 150, replace=TRUE)),
#'   Tm = c(chr1_tms, chr2_tms)
#' )
#'
#' # Plot single chromosome
#' plot_tm_genome_tracks(
#'   gr = tm_results,
#'   chromosome_to_plot = "chr1",
#'   genome_assembly = "hg19",
#'   tm_track_title = "DNA Sequence Tm",
#'   show_ideogram = TRUE
#' )
#'
#' # Plot both chromosomes
#' plot_tm_genome_tracks(
#'   gr = tm_results,
#'   chromosome_to_plot = c("chr1", "chr2"),
#'   genome_assembly = "hg19",
#'   tm_track_title = "DNA Sequence Tm",
#'   show_ideogram = TRUE,
#'   ncol = 2
#' )
#'
#' # Example with custom color palette
#' plot_tm_genome_tracks(
#'   gr = tm_results,
#'   chromosome_to_plot = c("chr1"),
#'   genome_assembly = "hg19",
#'   color_palette = "plasma",
#'   show_ideogram = TRUE,
#'   ncol = 2
#' )
#' }
#'
#' @importFrom Gviz IdeogramTrack GenomeAxisTrack DataTrack plotTracks
#' @importFrom GenomicRanges GRanges seqnames start end strand mcols
#' @importFrom IRanges IRanges
#' @importFrom dplyr %>%
#' @importFrom viridis viridis
#' @importFrom grDevices colorRampPalette dev.cur dev.new
#' @importFrom graphics layout par text
#' @importFrom GenomeInfoDb seqlengths seqlevels seqlevelsInUse genome
#' @export

plot_tm_genome_tracks <- function(gr,
                                  chromosome_to_plot,
                                  genome_assembly = NULL,
                                  tm_track_title = "Melting Temperature (\u00B0C)",
                                  color_palette = c("viridis", "magma", "plasma", "inferno", "cividis"),
                                  show_ideogram = TRUE,
                                  ncol = 2) {
  
  # Helper function to get viridis colors for Gviz
  .get_viridis_gviz_colors <- function(palette = c("viridis", "magma", "plasma", "inferno", "cividis")) {
    # Get the first and last colors from the viridis palette
    if (palette == "viridis") {
      colors <- viridis::viridis(2, option = "viridis")
    } else if (palette == "magma") {
      colors <- viridis::magma(2, option = "magma")
    } else if (palette == "plasma") {
      colors <- viridis::plasma(2, option = "plasma")
    } else if (palette == "inferno") {
      colors <- viridis::inferno(2, option = "inferno")
    } else if (palette == "cividis") {
      colors <- viridis::cividis(2, option = "cividis")
    } else {
      stop("Invalid palette. Please choose from: viridis, magma, plasma, inferno, cividis")
    }
    return(colors)
  }

  # Helper function to set default seqlengths if missing
  .set_default_seqlengths <- function(gr, genome_assembly = NULL) {
    # If seqlengths are already fully defined, no need to do anything
    if (!is.null(GenomeInfoDb::seqlengths(gr)) && all(!is.na(GenomeInfoDb::seqlengths(gr)))) {
      message("seqlengths are already defined for the GRanges object. Skipping default setting.")
      return(gr)
    }
    
    if (is.null(genome_assembly)) {
      stop("Please provide a 'genome_assembly' (e.g., 'hg38', 'mm10') to set default seqlengths.")
    }
    
    # Fetch sequence lengths for the specified genome assembly
    # This uses GenomeInfoDb, which is a standard bioconductor package
    # It will download the necessary sequence information if not already cached.
    tryCatch({
      fetched_lengths <- GenomeInfoDb::seqlengths(GenomeInfoDb::Seqinfo(genome = genome_assembly))      
      # Get unique chromosomes present in the GRanges object
      chr_in_gr <- unique(as.character(GenomicRanges::seqnames(gr)))
      
      # Filter fetched lengths to only include chromosomes present in your GRanges object
      # and that are available in the fetched data.
      valid_chrs_to_set <- intersect(chr_in_gr, names(fetched_lengths))
      
      if (length(valid_chrs_to_set) == 0) {
        warning(paste0("No matching chromosome lengths found for '", genome_assembly,
                      "' and the chromosomes in your GRanges object (",
                      paste(chr_in_gr, collapse = ", "), "). seqlengths not set."))
        return(gr) # Return GRanges as is if no lengths can be set
      }
      
      # Set seqlengths for only the chromosomes present in the data and available
      GenomeInfoDb::seqlengths(gr) <- fetched_lengths[valid_chrs_to_set]
      
      message("Setting default chromosome lengths for '", genome_assembly, "' for: ",
              paste(valid_chrs_to_set, collapse = ", "))
      
    }, error = function(e) {
      stop(paste0("Failed to fetch sequence lengths for genome assembly '", genome_assembly,
                  "'. Please check the assembly name or ensure Bioconductor data packages are installed. Error: ", e$message))
    })
    
    return(gr)
  }

  # --- Input Validation ---
  color_palette <- match.arg(color_palette)
  
  if (!is.character(genome_assembly) || length(genome_assembly) != 1) {
    stop("genome_assembly must be a single character string.")
  }
  
  if (!inherits(gr, "GRanges")) {
    stop("Input 'gr' must be a GRanges object.")
  }
  if (!"Tm" %in% names(mcols(gr))) {
    stop("GRanges object must have a 'Tm' metadata column.")
  }
  if (!all(chromosome_to_plot %in% unique(as.character(GenomicRanges::seqnames(gr))))) {
    stop(paste0("One or more chromosomes not found in the GRanges object: ",
                paste(setdiff(chromosome_to_plot, unique(as.character(GenomicRanges::seqnames(gr)))), collapse = ", ")))
  }
  
  # Set default seqlengths if missing
  gr <- .set_default_seqlengths(gr, genome_assembly)
  
  # --- Prepare Data ---
  tm_range <- range(mcols(gr)$Tm, na.rm = TRUE)
  if (length(tm_range) < 2 || is.infinite(tm_range[1]) || is.infinite(tm_range[2])) {
    warning("Tm range is not well-defined or contains only one value. Using a default color range.")
    tm_range <- c(50, 90)
  }
  
  # --- Derive Gviz colors from the common color_palette ---
  gviz_colors <- .get_viridis_gviz_colors(color_palette)
  tm_color_low_derived <- gviz_colors[1]
  tm_color_high_derived <- gviz_colors[2]
  
  # --- Create plotting function for a single chromosome ---
  create_chr_plot <- function(chr) {
    gr_filtered <- gr[seqnames(gr) == chr]
    
    # Ensure all ranges have the same strand for DataTrack
    GenomicRanges::strand(gr_filtered) <- "*"
    
    # Create tracks
    iTrack <- NULL
    if (show_ideogram && !is.null(genome_assembly)) {
      iTrack <- IdeogramTrack(genome = genome_assembly, chromosome = chr)
    }
    
    gTrack <- GenomeAxisTrack()
    
    dtTrack <- DataTrack(
      range = gr_filtered,
      data = mcols(gr_filtered)$Tm,
      genome = genome_assembly,
      chromosome = chr,
      name = paste(tm_track_title, "-", chr),
      type = "gradient",
      gradient = c(tm_color_low_derived, tm_color_high_derived),
      cex.axis = 0.8,
      ylim = tm_range,
      background.panel = "#FFF",
      background.title = "gray",
      cex.title = 1,
      col.axis = "black",
      color.scheme = list(
        "gradient",
        limits = tm_range,
        color = colorRampPalette(c(tm_color_low_derived, tm_color_high_derived))(100)
      )
    )
    
    track_list <- list(gTrack, dtTrack)
    if (!is.null(iTrack)) {
      track_list <- c(iTrack, track_list)
    }
    
    Gviz::plotTracks(track_list,
               chromosome = chr,
               from = min(GenomicRanges::start(gr_filtered)),
               to = max(GenomicRanges::end(gr_filtered)),
               title.width = 1.5
    )
  }
  
  # --- Handle multiple chromosomes ---
  if (length(chromosome_to_plot) > 1) {
    # Calculate layout dimensions
    n_plots <- length(chromosome_to_plot)
    nrow <- ceiling(n_plots / ncol)
    
    # Set up the plotting layout
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # Create layout matrix
    layout_matrix <- matrix(1:(nrow * ncol), nrow = nrow, ncol = ncol, byrow = TRUE)
    
    # Create a new plot device if none exists
    if (dev.cur() == 1) {
      dev.new()
    }
    
    # Set up the layout
    layout(layout_matrix)
    
    # Plot each chromosome
    for (chr in chromosome_to_plot) {
      tryCatch({
        create_chr_plot(chr)
      }, error = function(e) {
        message(sprintf("Error plotting chromosome %s: %s", chr, e$message))
        # Create an empty plot with error message
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, sprintf("Error plotting %s", chr), col = "red")
      })
    }
  } else {
    # Single chromosome plot
    create_chr_plot(chromosome_to_plot)
  }
  
  invisible(NULL)
}
