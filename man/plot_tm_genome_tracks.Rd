% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tm_genome_tracks.R
\name{plot_tm_genome_tracks}
\alias{plot_tm_genome_tracks}
\title{Plot Tm values as Genome Browser Tracks using Gviz}
\usage{
plot_tm_genome_tracks(
  gr,
  chromosome_to_plot,
  genome_assembly = NULL,
  tm_track_title = "Melting Temperature (°C)",
  color_palette = c("viridis", "magma", "plasma", "inferno", "cividis"),
  show_ideogram = TRUE,
  zoom = NULL
)
}
\arguments{
\item{gr}{A GRanges object. It MUST contain a metadata column named 'Tm'
with numeric melting temperature values.}

\item{chromosome_to_plot}{A character string specifying the chromosome to visualize.
These chromosomes must exist in your GRanges object.}

\item{genome_assembly}{A character string indicating the genome assembly (e.g., "hg19", "mm10").
This is used by IdeogramTrack for correct ideogram display.}

\item{tm_track_title}{A character string for the title of the Tm data track.}

\item{color_palette}{A character string specifying the viridis color palette to use.
Available options are:
\itemize{
  \item "viridis" (default): A perceptually uniform color map that works well for most people
  \item "magma": A perceptually uniform color map with a dark purple to bright yellow range
  \item "plasma": A perceptually uniform color map with a dark purple to bright yellow range
  \item "inferno": A perceptually uniform color map with a dark purple to bright yellow range
  \item "cividis": A perceptually uniform color map optimized for color vision deficiency
}
All palettes are colorblind-friendly and perceptually uniform.}

\item{show_ideogram}{Logical, whether to display the chromosome ideogram tracks.}

\item{zoom}{A character string specifying the genomic region to zoom into.
If NULL (default), the entire range of each chromosome will be shown.
Example: "chr1:1000000-2000000" for zooming into chr1:1000000-2000000}
}
\value{
Invisible NULL. The function generates a plot directly.
}
\description{
This function generates Gviz plots displaying Tm values as DataTracks
alongside genome axes and ideograms for specified chromosomes.
Tm values are visualized using a heatmap-like color gradient.
}
\examples{
\dontrun{
library(GenomicRanges)
library(Gviz)
# Example 1: Generate sample data with 150 sequences
set.seed(123)

# Generate 100 sequences for chr1
chr1_starts <- sort(sample(1:249250621, 100))  # chr1 length in hg19
chr1_lengths <- sample(50:200, 100, replace=TRUE)
chr1_ends <- chr1_starts + chr1_lengths
chr1_tms <- runif(100, min=60, max=80)

# Generate 50 sequences for chr2
chr2_starts <- sort(sample(1:243199373, 50))   # chr2 length in hg19
chr2_lengths <- sample(50:200, 50, replace=TRUE)
chr2_ends <- chr2_starts + chr2_lengths
chr2_tms <- runif(50, min=60, max=80)

# Create GRanges object
tm_results <- GRanges(
  seqnames = Rle(c(rep("chr1", 100), rep("chr2", 50))),
  ranges = IRanges(
    start = c(chr1_starts, chr2_starts),
    end = c(chr1_ends, chr2_ends)
  ),
  strand = Rle(sample(c("+", "-"), 150, replace=TRUE)),
  Tm = c(chr1_tms, chr2_tms)
)

# Plot single chromosome with zoom
plot_tm_genome_tracks(
  gr = tm_results,
  chromosome_to_plot = "chr1",
  genome_assembly = "hg19",
  tm_track_title = "DNA Sequence Tm",
  zoom = "chr1:10062800-20000000"
)

# Example with custom color palette and no zoom
plot_tm_genome_tracks(
  gr = tm_results,
  chromosome_to_plot = "chr2",
  genome_assembly = "hg19",
  color_palette = "plasma"
)
}

}
