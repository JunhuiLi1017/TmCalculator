#' Launch the TmCalculatorShiny application
#'
#' Starts the interactive Shiny front end for the \pkg{TmCalculator} package.
#' The app provides four panels: \emph{Sequence Tm} (calculate melting
#' temperatures from pasted sequences, genomic coordinates, FASTA uploads, or
#' genome-wide BSgenome tiling), \emph{Visualize} (Tm plots and multi-omics
#' multi-track genome views), \emph{Compare groups} (group-wise statistical
#' comparison of stored results), and \emph{About}.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}
#'   (e.g. \code{port}, \code{host}, \code{launch.browser}).
#'
#' @return Invisibly returns \code{NULL}; called for the side effect of
#'   launching the application.
#'
#' @examples
#' \dontrun{
#' runTmCalculatorShiny()
#' }
#'
#' @export
runTmCalculatorShiny <- function(...) {
  if (!.tcs_has_pkg("shiny")) {
    stop("Package 'shiny' is required to run the application.")
  }
  app <- shiny::shinyApp(ui = app_ui(), server = app_server)
  shiny::runApp(app, ...)
}
