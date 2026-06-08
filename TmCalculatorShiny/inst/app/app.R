# Entry point for `shiny::runApp(system.file("app", package = "TmCalculatorShiny"))`
library(TmCalculatorShiny)
shiny::shinyApp(ui = TmCalculatorShiny:::app_ui(),
                server = TmCalculatorShiny:::app_server)
