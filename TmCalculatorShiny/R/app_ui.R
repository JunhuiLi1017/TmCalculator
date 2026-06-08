#' Assemble the application UI
#' @noRd
app_ui <- function() {
  bslib::page_navbar(
    title = "TmCalculatorShiny",
    theme = bslib::bs_theme(version = 5, preset = "cosmo"),
    id = "main_nav",
    bslib::nav_panel(
      title = "Sequence Tm",
      mod_sequence_ui("seq")
    ),
    bslib::nav_panel(
      title = "Visualize",
      mod_visualize_ui("viz")
    ),
    bslib::nav_panel(
      title = "Compare groups",
      mod_compare_ui("cmp")
    ),
    bslib::nav_panel(
      title = "About",
      mod_about_ui("about")
    )
  )
}
