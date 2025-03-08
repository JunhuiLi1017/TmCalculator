# Load required libraries
library(shiny)
library(TmCalculator)

# Define UI
ui <- fluidPage(
  titlePanel("DNA Melting Temperature (Tm) Calculator"),
  sidebarLayout(
    sidebarPanel(
      textInput("sequence", "DNA Sequence", value = "ATGCATGCATGC"),
      numericInput("sodium", "Sodium Concentration (mM)", value = 50, min = 0),
      numericInput("mg", "Magnesium Concentration (mM)", value = 0, min = 0),
      numericInput("dntp", "dNTP Concentration (mM)", value = 0, min = 0),
      selectInput("method", "Tm Calculation Method", 
                  choices = c("Breslauer (1986)" = "breslauer",
                              "SantaLucia (1998)" = "santalucia",
                              "Owczarzy (2004)" = "owczarzy")),
      actionButton("calculate", "Calculate Tm")
    ),
    mainPanel(
      h3("Results"),
      verbatimTextOutput("tm_result")
    )
  )
)

# Define server logic
server <- function(input, output) {
  observeEvent(input$calculate, {
    # Extract inputs
    sequence <- toupper(input$sequence)
    sodium <- input$sodium
    mg <- input$mg
    dntp <- input$dntp
    method <- input$method
    
    # Validate sequence
    if (!grepl("^[ATGC]+$", sequence)) {
      output$tm_result <- renderText({
        "Invalid DNA sequence. Please enter a sequence containing only A, T, G, or C."
      })
      return()
    }
    
    # Calculate Tm
    tm <- TmCalculator::Tm(sequence, Na = sodium, Mg = mg, dNTPs = dntp, method = method)
    
    # Display result
    output$tm_result <- renderText({
      paste("Melting Temperature (Tm):", round(tm, 2), "°C")
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)