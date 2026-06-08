#' About panel UI
#' @noRd
mod_about_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header("About TmCalculatorShiny"),
    shiny::tags$div(
      shiny::tags$p(
        "TmCalculatorShiny is an interactive front end for the ",
        shiny::tags$b("TmCalculator"), " R/Bioconductor package, for ",
        "nucleic-acid melting-temperature analysis."
      ),
      shiny::tags$h5("Workflow"),
      shiny::tags$ol(
        shiny::tags$li(
          shiny::tags$b("Sequence Tm"), " - compute Tm from pasted sequences, ",
          "genomic coordinates, an uploaded FASTA, or genome-wide tiling of a ",
          "BSgenome package. Each calculation is saved to a shared store."
        ),
        shiny::tags$li(
          shiny::tags$b("Visualize"), " - plot stored Tm/GC results and, ",
          "optionally, integrate multi-omics feature ranges into multi-track ",
          "linear, circular, or karyotype genome views."
        ),
        shiny::tags$li(
          shiny::tags$b("Compare groups"), " - attach grouping metadata to a ",
          "stored result and test whether Tm/GC differ across groups."
        )
      ),
      shiny::tags$h5("Methods"),
      shiny::tags$p(
        "Tm is computed with TmCalculator's nearest-neighbour (tm_nn), ",
        "GC-based (tm_gc, including custom A/B/C/D coefficients), and ",
        "Wallace-rule (tm_wallace) engines."
      ),
      shiny::tags$p(
        shiny::tags$small(
          "Author: Junhui Li (ORCID 0000-0003-3973-1700). Co-author: ",
          "Lihua Julie Zhu."
        )
      )
    )
  )
}

#' About panel server
#' @noRd
mod_about_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    invisible(NULL)
  })
}
