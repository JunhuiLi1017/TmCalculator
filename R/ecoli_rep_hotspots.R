#' E. coli K-12 MG1655 replication-associated hotspot annotations
#'
#' A named list of genomic feature tables for *Escherichia coli* K-12 MG1655
#' (NCBI assembly GCF_000005845.2 / ASM584v2, chromosome \code{U00096.3}).
#' The data support the genome-wide Tm vignette and
#' \code{\link{plot_genome_track}} examples by providing independent multi-omics
#' layers that can be overlaid on the same genomic axis alongside Tm profiles
#' computed with \code{\link{tm_calculate}}.
#'
#' @format A named list with five elements:
#' \describe{
#'   \item{\code{all_peaks_IP_mutH}}{A data frame (38 rows) of MutL protein
#'     ChIP-seq peaks marking mismatch-repair-associated regions (MutL-AR).
#'     Columns: \code{chr}, \code{start}, \code{end}, \code{Sample}, \code{name},
#'     \code{Size}.}
#'   \item{\code{bins_rep}}{A data frame (4,642 rows) of tandem-repeat
#'     (microsatellite) counts in non-overlapping 1 kb bins. Columns:
#'     \code{chr}, \code{start}, \code{end}, \code{count}.}
#'   \item{\code{bins_cru}}{A data frame (4,642 rows) of cruciform-forming
#'     sequence counts in non-overlapping 1 kb bins. Columns: \code{chr},
#'     \code{start}, \code{end}, \code{count}.}
#'   \item{\code{ssdna}}{A data frame (2,636 rows) of single-stranded DNA
#'     regions. Columns: \code{chr}, \code{start}, \code{end},
#'     \code{Cells..strand.}, \code{Region}.}
#'   \item{\code{bins_gatc}}{A data frame (4,642 rows) of GATC methylation-site
#'     counts in non-overlapping 1 kb bins. Columns: \code{chr}, \code{start},
#'     \code{end}, \code{count}.}
#' }
#'
#' All coordinate-based tables use \code{chr = "U00096.3"} and are compatible
#' with \code{\link{plot_genome_track}} and \code{\link{compare_groups}}.
#'
#' @source Curated from published *E. coli* K-12 MG1655 multi-omics datasets
#'   used in the package vignette
#'   \code{vignette("genome_wide_tm_ecoli", package = "TmCalculator")}.
#'
#' @examples
#' data(ecoli_rep_hotspots)
#' names(ecoli_rep_hotspots)
#'
#' # MutL-AR peak coordinates
#' head(ecoli_rep_hotspots$all_peaks_IP_mutH)
#'
#' # Microsatellite density in 1 kb bins
#' summary(ecoli_rep_hotspots$bins_rep$count)
"ecoli_rep_hotspots"
