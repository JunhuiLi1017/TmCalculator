#' Plot Tm Values from GRanges with Per-Chromosome Colors and Shapes
#'
#' Creates a genome-wide plot of melting temperature (Tm) values from a GRanges object
#' using the karyoploteR package. The x-axis represents genomic positions across
#' chromosomes, and the y-axis represents Tm values. Points are plotted at the
#' midpoints of genomic ranges, with customizable colors and shapes per chromosome.
#'
#' @param gr A \code{GRanges} object with a \code{Tm} metadata column containing
#'   numeric melting temperature values.
#' @param chromosomes A character vector specifying chromosomes to plot. If \code{NULL},
#'   all unique chromosomes in \code{gr} are plotted. Defaults to \code{NULL}.
#' @param genome_assembly It can be either a UCSC style genome name (hg19, mm10, etc),
#'  a BSgenome, a Seqinfo object, a GRanges object with the chromosomes as ranges or 
#' in general any genome specification accepted by karyoploteR. If \code{NULL}, uses 
#' default or GRanges seqinfo. Defaults to \code{NULL}.
#' @param colors A named character vector specifying colors for each chromosome
#'   (e.g., \code{c(chr1 = "#FF0000", chr22 = "#00FF00")}). Names must match
#'   chromosomes in \code{gr} or \code{chromosomes}. If \code{NULL} or partially
#'   specified, unspecified chromosomes use the first viridis color. Defaults to \code{NULL}.
#' @param shapes A named integer vector specifying point shapes (pch values) for each
#'   chromosome (e.g., \code{c(chr1 = 16, chr22 = 17)}). Names must match
#'   chromosomes in \code{gr} or \code{chromosomes}. If \code{NULL} or partially
#'   specified, unspecified chromosomes use \code{pch = 16} (filled circles).
#'   Defaults to \code{NULL}.
#' @param plot_type An integer specifying the karyoploteR plot type (e.g., 1 for horizontal
#'   chromosomes, 4 or 7 for vertical or grid layouts). See \code{\link[karyoploteR]{plotKaryotype}}
#'   for details. Defaults to 1.
#' @param point_cex A numeric value for the size of plotted points. Defaults to 1.5.
#' @param xaxis_cex A numeric value for the text size of x-axis labels (base pair positions).
#'   Defaults to 0.7.
#' @param yaxis_cex A numeric value for the text size of y-axis labels (Tm values).
#'   Defaults to 0.8.
#' @param chr_cex A numeric value for the text size of chromosome names.
#'   Defaults to 1.
#' @param tick_dist A numeric value for the distance between tick marks on the x-axis.
#'   Defaults to 10000000.
#' @param zoom A \code{GRanges} object specifying a genomic region to zoom into.
#'   If \code{NULL}, the full chromosomes are plotted. Defaults to \code{NULL}.
#'
#' @return Invisibly returns \code{NULL}. The function generates a plot as a side effect.
#'
#' @details
#' The function validates that \code{gr} is a \code{GRanges} object with a \code{Tm} column.
#' It automatically sets the y-axis limits based on the range of \code{Tm} values, with
#' slight padding (floor and ceiling). The plot includes chromosome names (via karyoploteR's
#' default labeling), base pair positions, and a labeled y-axis. Points are colored and shaped
#' according to the \code{colors} and \code{shapes} parameters, with defaults applied for
#' unspecified chromosomes. The y-axis is placed on the left for \code{plot_type = 1} and on
#' the right for \code{plot_type = 4} or \code{7}, with the label positioned to the right of
#' the y-axis, vertically centered in the middle of the y-axis range. Text sizes for x-axis
#' labels, y-axis labels, and chromosome names can be customized using \code{xaxis_cex},
#' \code{yaxis_cex}, and \code{chr_cex}, respectively.
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' # Create a sample GRanges object
#' gr <- GRanges(
#'   seqnames = c("chr22", "chr1", "chr14", "chr22"),
#'   ranges = IRanges(
#'     start = c(13209021, 1, 13200, 13209150),
#'     end = c(13209099, 76, 13222, 13209200)
#'   ),
#'   strand = c("+", "*", "*", "-"),
#'   Tm = c(69.1147, 71.1160, 50.7169, 65.5000)
#' )
#' genome(seqinfo(gr)) <- "hg19"
#'
#' # Plot with default settings (plot_type=1)
#' plot_tm_karyotype(gr)
#'
#' # Plot with partial color and shape specification (plot_type=4)
#' plot_tm_karyotype(
#'   gr,
#'   genome_assembly="hg19",
#'   colors = c(chr1 = "#FF0000"),  # chr14, chr22 use default color
#'   shapes = c(chr1 = 16, chr14 = 17),  # chr22 uses pch=16
#'   plot_type = 4,
#'   xaxis_cex = 0.6,
#'   yaxis_cex = 0.9,
#'   chr_cex = 1.2
#' )
#'
#' # Plot with full color and shape specification (plot_type=7)
#' plot_tm_karyotype(
#'   gr,
#'   genome_assembly="hg38",
#'   colors = c(chr1 = "#FF0000", chr14 = "#00FF00", chr22 = "#0000FF"),
#'   shapes = c(chr1 = 16, chr14 = 17, chr22 = 16),
#'   plot_type = 5,
#'   xaxis_cex = 0.8,
#'   yaxis_cex = 0.7,
#'   chr_cex = 0.8
#' )
#'
#' # Plot with zoom into chr22
#' zoom_region <- GRanges("chr22:13200000-13220000")
#' plot_tm_karyotype(
#'   gr,
#'   genome_assembly="hg38",
#'   chromosomes = "chr22",
#'   zoom = zoom_region,
#'   xaxis_cex = 0.5,
#'   yaxis_cex = 1,
#'   chr_cex = 1.5
#' )
#' }
#'
#' @importFrom karyoploteR plotKaryotype kpDataBackground kpAxis kpPoints kpAddBaseNumbers
#' @importFrom GenomicRanges GRanges seqnames start end mcols
#' @importFrom viridis viridis
#' @export

plot_tm_karyotype <- function(gr, 
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
  # --- Input Validation ---
  if (!inherits(gr, "GRanges")) {
    stop("Input 'gr' must be a GRanges object.")
  }
  if (!"Tm" %in% names(mcols(gr))) {
    stop("GRanges object must have a 'Tm' metadata column.")
  }
# --- UPDATED genome_assembly HANDLING ---
  actual_genome_assembly_name <- NULL

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
  
  # Get chromosomes to plot
  if (is.null(chromosomes)) {
    chromosomes <- unique(as.character(seqnames(gr)))
  } else if (!all(chromosomes %in% unique(as.character(seqnames(gr))))) {
    stop(paste0("One or more chromosomes not found in the GRanges object: ",
                paste(setdiff(chromosomes, unique(as.character(seqnames(gr)))), collapse = ", ")))
  }
  
  # Set Tm range for y-axis
  tm_range <- range(mcols(gr)$Tm, na.rm = TRUE)
  ymin <- floor(min(tm_range))
  ymax <- ceiling(max(tm_range))
  if (ymin >= ymax) {
    warning("Invalid Tm range. Using default range 50-90.")
    ymin = 50
    ymax = 90
  }
  
  # Set default colors and shapes
  default_color <- "black"
  default_shape <- 16  # Filled circle
  color_map <- rep(default_color, length(chromosomes))
  shape_map <- rep(default_shape, length(chromosomes))
  names(color_map) <- chromosomes
  names(shape_map) <- chromosomes
  
  # Handle colors
  if (!is.null(colors)) {
    if (!is.character(colors) || is.null(names(colors))) {
      stop("'colors' must be a named character vector with chromosome names.")
    }
    valid_colors <- colors[names(colors) %in% chromosomes]
    if (length(valid_colors) == 0) {
      warning("No valid chromosome names in 'colors'. Using default color.")
    } else {
      color_map[names(valid_colors)] <- valid_colors
      if (length(valid_colors) < length(chromosomes)) {
        warning("Fewer colors than chromosomes. Using default color for unspecified chromosomes.")
        color_map[setdiff(chromosomes, names(valid_colors))] <- default_color
      }
    }
  }
  
  # Handle shapes
  if (!is.null(shapes)) {
    if (!is.numeric(shapes) || is.null(names(shapes))) {
      stop("'shapes' must be a named integer vector with chromosome names.")
    }
    valid_shapes <- shapes[names(shapes) %in% chromosomes]
    if (length(valid_shapes) == 0) {
      warning("No valid chromosome names in 'shapes'. Using default shape (pch=16).")
    } else {
      shape_map[names(valid_shapes)] <- valid_shapes
      if (length(valid_shapes) < length(chromosomes)) {
        warning("Fewer shapes than chromosomes. Using default shape (pch=16) for unspecified chromosomes.")
        shape_map[setdiff(chromosomes, names(valid_shapes))] <- default_shape
      }
    }
  }
  
  # --- Create Karyotype Plot ---
  kp <- plotKaryotype(
    chromosomes = chromosomes,
    plot.type = plot_type,
    genome = actual_genome_assembly_name,
    zoom = zoom,
    cex = chr_cex
  )
  
  # Add background panel
  kpDataBackground(kp, data.panel = 1, color = "#F5F5F5")
  
  # Add y-axis for Tm
  kpAxis(
    kp,
    ymin = ymin,
    ymax = ymax,
    numticks = 5,
    side = ifelse(plot_type %in% c(4, 7), 4, 2),  # Right for plot_type 4, 7; left for others
    cex = yaxis_cex
  )
  
  # Plot Tm data as points with per-chromosome colors and shapes
  for (chr in chromosomes) {
    gr_chr <- gr[seqnames(gr) == chr]
    if (length(gr_chr) > 0) {
      kpPoints(
        kp,
        data = gr_chr,
        x = (start(gr_chr) + end(gr_chr)) / 2,  # Midpoint of ranges
        y = mcols(gr_chr)$Tm,
        ymin = ymin,
        ymax = ymax,
        col = color_map[chr],
        pch = as.integer(shape_map[chr]),
        cex = point_cex
      )
    }
  }
  
  # Add base numbers
  kpAddBaseNumbers(kp, tick.dist = tick_dist, cex = xaxis_cex)
  
  invisible(NULL)
}