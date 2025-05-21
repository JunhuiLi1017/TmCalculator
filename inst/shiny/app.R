library(shiny)
library(TmCalculator)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "Tm Calculator"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Calculator", tabName = "calculator", icon = icon("calculator")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    
    tabItems(
      # Calculator Tab
      tabItem(tabName = "calculator",
        fluidRow(
          box(
            title = "Sequence Input",
            width = 6,
            selectInput("method", "Calculation Method",
                       choices = list(
                         "Nearest Neighbor (Tm_NN)" = "Tm_NN",
                         "GC Content (Tm_GC)" = "Tm_GC",
                         "Wallace Rule (Tm_wallace)" = "Tm_wallace"
                       ),
                       selected = "Tm_NN"),
            radioButtons("input_type", "Input Type",
                        choices = list(
                          "Direct Input" = "direct",
                          "FASTA File" = "fasta"
                        ),
                        selected = "direct"),
            conditionalPanel(
              condition = "input.input_type == 'direct'",
              textInput("primer_seq", "Primer Sequence (5' to 3')", 
                       placeholder = "Enter primer sequence..."),
              textInput("template_seq", "Template Sequence (5' to 3')", 
                       placeholder = "Enter template sequence...")
            ),
            conditionalPanel(
              condition = "input.input_type == 'fasta'",
              fileInput("fasta_file", "Upload FASTA File",
                       accept = c(".fasta", ".fa", ".txt")),
              selectInput("fasta_primer", "Select Primer Sequence",
                         choices = NULL),
              selectInput("fasta_template", "Select Template Sequence",
                         choices = NULL)
            ),
            checkboxInput("ambiguous", "Handle Ambiguous Bases", value = FALSE),
            actionButton("calculate", "Calculate Tm", class = "btn-primary")
          ),
          
          box(
            title = "Salt Conditions",
            width = 6,
            numericInput("Na", "Na+ (mM)", value = 50, min = 0),
            numericInput("K", "K+ (mM)", value = 0, min = 0),
            numericInput("Tris", "Tris (mM)", value = 0, min = 0),
            numericInput("Mg", "Mg2+ (mM)", value = 0, min = 0),
            numericInput("dNTPs", "dNTPs (mM)", value = 0, min = 0)
          )
        ),
        
        fluidRow(
          box(
            title = "Chemical Modifications",
            width = 6,
            numericInput("DMSO", "DMSO (%)", value = 0, min = 0, max = 100),
            numericInput("formamide", "Formamide", value = 0, min = 0),
            selectInput("formamide_unit", "Formamide Unit",
                       choices = list("Percent (%)" = "percent", "Molar (M)" = "molar"),
                       selected = "percent")
          ),
          
          box(
            title = "Results",
            width = 6,
            verbatimTextOutput("results"),
            downloadButton("download_results", "Download Results")
          )
        )
      ),
      
      # About Tab
      tabItem(tabName = "about",
        box(
          title = "About Tm Calculator",
          width = 12,
          includeMarkdown(system.file("shiny", "about.md", package = "TmCalculator"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Function to read FASTA file
  read_fasta <- function(file_path) {
    lines <- readLines(file_path)
    sequences <- list()
    current_seq <- NULL
    current_name <- NULL
    
    for (line in lines) {
      if (grepl("^>", line)) {
        if (!is.null(current_seq)) {
          sequences[[current_name]] <- current_seq
        }
        current_name <- sub("^>", "", line)
        current_seq <- ""
      } else {
        current_seq <- paste0(current_seq, line)
      }
    }
    if (!is.null(current_seq)) {
      sequences[[current_name]] <- current_seq
    }
    return(sequences)
  }
  
  # Update sequence selectors when FASTA file is uploaded
  observeEvent(input$fasta_file, {
    req(input$fasta_file)
    sequences <- read_fasta(input$fasta_file$datapath)
    updateSelectInput(session, "fasta_primer",
                     choices = names(sequences),
                     selected = names(sequences)[1])
    updateSelectInput(session, "fasta_template",
                     choices = names(sequences),
                     selected = if(length(sequences) > 1) names(sequences)[2] else names(sequences)[1])
  })
  
  # Calculate Tm when button is clicked
  observeEvent(input$calculate, {
    tryCatch({
      # Get sequences based on input type
      if (input$input_type == "direct") {
        primer <- input$primer_seq
        template <- input$template_seq
      } else {
        req(input$fasta_file)
        sequences <- read_fasta(input$fasta_file$datapath)
        primer <- sequences[[input$fasta_primer]]
        template <- sequences[[input$fasta_template]]
      }
      
      # Validate sequences
      if (primer == "" || template == "") {
        showNotification("Please enter both primer and template sequences", type = "error")
        return()
      }
      
      # Calculate Tm based on selected method
      result <- switch(input$method,
        "Tm_NN" = Tm_NN(
          input_seq = primer,
          complement_seq = template,
          Na = input$Na,
          K = input$K,
          Tris = input$Tris,
          Mg = input$Mg,
          dNTPs = input$dNTPs,
          DMSO = input$DMSO,
          formamide = list(value = input$formamide, unit = input$formamide_unit),
          ambiguous = input$ambiguous
        ),
        "Tm_GC" = Tm_GC(
          input_seq = primer,
          Na = input$Na,
          K = input$K,
          Tris = input$Tris,
          Mg = input$Mg,
          dNTPs = input$dNTPs,
          DMSO = input$DMSO,
          formamide = list(value = input$formamide, unit = input$formamide_unit),
          ambiguous = input$ambiguous
        ),
        "Tm_wallace" = Tm_wallace(
          input_seq = primer
        )
      )
      
      # Display results
      output$results <- renderPrint({
        cat("Method:", input$method, "\n")
        if (input$input_type == "fasta") {
          cat("Primer:", input$fasta_primer, "\n")
          cat("Template:", input$fasta_template, "\n")
        }
        print(result)
      })
      
      # Enable download button
      output$download_results <- downloadHandler(
        filename = function() {
          paste("tm_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt", sep = "")
        },
        content = function(file) {
          writeLines(capture.output(
            cat("Method:", input$method, "\n"),
            if (input$input_type == "fasta") {
              cat("Primer:", input$fasta_primer, "\n")
              cat("Template:", input$fasta_template, "\n")
            },
            print(result)
          ), file)
        }
      )
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server) 