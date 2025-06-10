#' Launch TmCalculator Shiny Application
#'
#' Launches an interactive Shiny application for calculating melting temperatures (Tm) of DNA/RNA sequences.
#' The application provides a user-friendly interface for sequence analysis, Tm calculation, and visualization.
#'
#' @details The application provides comprehensive options for Tm calculations:
#' 
#' \strong{Input Methods:}
#' \itemize{
#'   \item Direct sequence input
#'   \item Genomic coordinate input (e.g., "chr1:100-200:+:hg38")
#'   \item FASTA file upload
#'   \item Complementary sequence input
#' }
#' 
#' \strong{Calculation Methods:}
#' \itemize{
#'   \item Nearest Neighbor thermodynamics (tm_nn):
#'     \itemize{
#'       \item Multiple parameter sets (DNA/DNA, RNA/RNA, RNA/DNA)
#'       \item Terminal and internal mismatch handling
#'       \item Dangling end corrections
#'     }
#'   \item GC content-based method (tm_gc):
#'     \itemize{
#'       \item Multiple variants (Primer3Plus, Chester1993, etc.)
#'       \item Custom coefficient sets
#'       \item Mismatch handling
#'     }
#'   \item Wallace rule (tm_wallace):
#'     \itemize{
#'       \item Simple 2°C per A/T, 4°C per G/C calculation
#'       \item Suitable for primers 14-20 nt in length
#'     }
#' }
#' 
#' \strong{Additional Features:}
#' \itemize{
#'   \item Salt concentration corrections
#'   \item Chemical modifications (DMSO, formamide)
#'   \item Ambiguous base handling
#'   \item Sequence alignment options
#'   \item Results export
#' }
#'
#' @return Launches the Shiny application in the user's default web browser.
#'
#' @examples
#' \dontrun{
#' # Launch the TmCalculator Shiny application
#' TmCalculatorShinyApp()
#' }
#'
#' @seealso \code{\link{tm_calculate}} for the core Tm calculation function
#' @seealso \code{\link{tm_nn}} for nearest neighbor thermodynamics
#' @seealso \code{\link{tm_gc}} for GC content-based calculations
#' @seealso \code{\link{tm_wallace}} for Wallace rule calculations
#'
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @import  shiny
#' @import  shinydashboard
#' @export
#'
TmCalculatorShinyApp <- function() {
  runApp(appDir = system.file('shiny', package = 'TmCalculator'))
}
