#' Application server with a shared result store
#' @noRd
app_server <- function(input, output, session) {

  # Shared store of computed results. Each entry is a named list element:
  #   store$results[[label]] <- list(gr = <GRanges with Tm/GC>, options = <list>)
  store <- shiny::reactiveValues(results = list())

  # Register a new result and make it available to the other modules.
  add_result <- function(label, gr, options = NULL) {
    label <- as.character(label)
    if (is.null(label) || !nzchar(label)) {
      label <- paste0("result_", length(store$results) + 1L)
    }
    # Ensure uniqueness
    base_label <- label
    k <- 1L
    while (label %in% names(store$results)) {
      k <- k + 1L
      label <- paste0(base_label, " (", k, ")")
    }
    store$results[[label]] <- list(gr = gr, options = options)
    label
  }

  mod_sequence_server("seq", add_result)
  mod_visualize_server("viz", store)
  mod_compare_server("cmp", store)
  mod_about_server("about")
}
