#' Plot Tm values as Genome Browser Tracks using Gviz
#'
#' This function generates Gviz plots displaying Tm values as DataTracks
#' alongside genome axes and ideograms for specified chromosomes.
#' Tm values are visualized using a heatmap-like color gradient.
#'
#' @param gr A GRanges object. It MUST contain a metadata column named 'Tm'
#'           with numeric melting temperature values.
#' @param chromosome_to_plot A character string specifying the chromosome to visualize.
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
#' @param zoom A character string specifying the genomic region to zoom into.
#'            If NULL (default), the entire range of each chromosome will be shown.
#'            Example: "chr1:1000000-2000000" for zooming into chr1:1000000-2000000
#' @return Invisible NULL. The function generates a plot directly.
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' library(Gviz)
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
#' # Plot single chromosome with zoom
#' plot_tm_genome_tracks(
#'   gr = tm_results,
#'   chromosome_to_plot = "chr1",
#'   genome_assembly = "hg19",
#'   tm_track_title = "DNA Sequence Tm",
#'   zoom = "chr1:10062800-20000000"
#' )
#'
#' # Example with custom color palette and no zoom
#' plot_tm_genome_tracks(
#'   gr = tm_results,
#'   chromosome_to_plot = "chr2",
#'   genome_assembly = "hg19",
#'   color_palette = "plasma"
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
                                  zoom = NULL) {
  
  # Helper function to get viridis colors for Gviz
  .get_viridis_gviz_colors <- function(color_palette = c("viridis", "magma", "plasma", "inferno", "cividis")) {
    # Get the first and last colors from the viridis palette
    color_palette <- match.arg(color_palette)
    if (color_palette == "viridis") {
      colors <- viridis::viridis(2)
    } else if (color_palette == "magma") {
      colors <- viridis::magma(2)
    } else if (color_palette == "plasma") {
      colors <- viridis::plasma(2)
    } else if (color_palette == "inferno") {
      colors <- viridis::inferno(2)
    } else if (color_palette == "cividis") {
      colors <- viridis::cividis(2)
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
  if (length(chromosome_to_plot) > 1) {
    warning("Please provide a single chromosome to plot, only the first chromosome will be plotted. If you want to plot multiple chromosomes, use the 'chromosome_to_plot' argument multiple times.")
    chromosome_to_plot <- chromosome_to_plot[1]
  }

  if (!all(chromosome_to_plot %in% unique(as.character(GenomicRanges::seqnames(gr))))) {
    stop(paste0("Chromosomes not found in the GRanges object: ",
                paste(setdiff(chromosome_to_plot, unique(as.character(GenomicRanges::seqnames(gr)))), collapse = ", "),
                ". Please provide a valid chromosome name."))
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
        color = colorRampPalette(c(tm_color_low_derived, tm_color_high_derived))(nrow(gr_filtered))
      )
    )
    
    track_list <- list(gTrack, dtTrack)
    if (!is.null(iTrack)) {
      track_list <- c(iTrack, track_list)
    }
    
    # Determine plot range
    if (!is.null(zoom)) {
      # Validate zoom parameter
      if (!is.character(zoom) || !grepl("^chr[0-9]+:[0-9]+-[0-9]+$", zoom)) {
        stop("zoom must be a character string like 'chr1:1000000-2000000'")
      }
      if (length(zoom) != 1) {
        warning("zoom must be a character string like 'chr1:1000000-2000000', only the first chromosome will be plotted. If you want to plot multiple chromosomes, use the 'zoom' argument multiple times.")
        zoom <- zoom[1]
      }
      zoom_range <- strsplit(zoom, ":")[[1]]
      chr_zoom <- zoom_range[1]
      if (chr_zoom != chromosome_to_plot) {
        stop(paste0("zoom chromosome '", chr_zoom, "' does not match the chromosome to plot '", chromosome_to_plot, "'. Please provide a valid chromosome name."))
      }
      zoom_range_pos <- as.numeric(strsplit(zoom_range[2], "-")[[1]])
      from <- zoom_range_pos[1]
      to <- zoom_range_pos[2]

    } else {
      from <- min(GenomicRanges::start(gr_filtered))
      to <- max(GenomicRanges::end(gr_filtered))
    }
    
    Gviz::plotTracks(track_list,
               chromosome = chr,
               from = from,
               to = to,
               title.width = 1.5
    )
  }
  
  # --- Handle multiple chromosomes ---
  create_chr_plot(chromosome_to_plot)
  
  invisible(NULL)
}
