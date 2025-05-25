library(shiny)
library(TmCalculator)
library(shinydashboard)

# Ensure the package is properly loaded
if (!requireNamespace("TmCalculator", quietly = TRUE)) {
  stop("TmCalculator package is required but not installed")
}

ui <- dashboardPage(
  dashboardHeader(
    title = tags$a(href = "https://github.com/JunhuiLi1017/TmCalculator", 
                   "TmCalculator", 
                   style = "color: black; font-family: Arial, Helvetica, sans-serif;"),
    titleWidth = 200
  ),
  
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
                # Left column - Input
                box(
                  title = "Sequence Input",
                  width = 6,
                  selectInput("method", "Calculation Method",
                              choices = list(
                                "Nearest Neighbor (tm_nn)" = "tm_nn",
                                "GC Content (tm_gc)" = "tm_gc",
                                "Wallace Rule (tm_wallace)" = "tm_wallace"
                              ),
                              selected = "tm_nn"),
                  radioButtons("input_type", "Input Type",
                               choices = list(
                                 "Direct Input" = "direct",
                                 "FASTA File" = "fasta"
                               ),
                               selected = "direct"),
                  conditionalPanel(
                    condition = "input.input_type == 'direct' && input.method == 'tm_nn'",
                    textInput("input_seq", "Primer Sequence (5' to 3')",
                              placeholder = "for example AACAGACT or AACAGACT, CGTGCATG"),
                    textInput("rev_input_seq", "Template Sequence (5' to 3')",
                              placeholder = "for example AGTCTGTT or AGTCTGTT, CATGCACG")
                  ),
                  conditionalPanel(
                    condition = "input.input_type == 'direct' && (input.method == 'tm_gc' || input.method == 'tm_wallace')",
                    textInput("input_seq", "Sequence (5' to 3')",
                              placeholder = "for example AACAGACT or AACAGACT, CGTGCATG")
                  ),
                  conditionalPanel(
                    condition = "input.input_type == 'fasta' && input.method == 'tm_nn'",
                    fileInput("fasta_file", "Upload FASTA File",
                              accept = c(".fasta", ".fa", ".txt")),
                    fileInput("fasta_file_complement", "Upload reverse complement FASTA File",
                              accept = c(".fasta", ".fa", ".txt"))
                  ),
                  conditionalPanel(
                    condition = "input.input_type == 'fasta' && (input.method == 'tm_gc' || input.method == 'tm_wallace')",
                    fileInput("fasta_file", "Upload FASTA File",
                              accept = c(".fasta", ".fa", ".txt"))
                  ),
                  actionButton("calculate", "Calculate Tm", class = "btn-primary")
                ),
                
                # Right column - Parameters
                box(
                  title = "Specific Parameters",
                  width = 6,
                  checkboxInput("ambiguous", "Ambiguous Bases", value = FALSE),
                  conditionalPanel(
                    condition = "input.method == 'tm_gc'",
                    selectInput("variant", "Variant",
                                choices = list(
                                  "Primer3Plus" = "Primer3Plus",
                                  "Chester1993" = "Chester1993",
                                  "QuikChange" = "QuikChange",
                                  "Schildkraut1965" = "Schildkraut1965",
                                  "Wetmur1991_MELTING" = "Wetmur1991_MELTING",
                                  "Wetmur1991_RNA" = "Wetmur1991_RNA",
                                  "Wetmur1991_RNA/DNA" = "Wetmur1991_RNA/DNA",
                                  "vonAhsen2001" = "vonAhsen2001"
                                ),
                                selected = "Primer3Plus")
                  ),
                  conditionalPanel(
                    condition = "input.method == 'tm_nn'",
                    numericInput("shift", "Alignment offset between primer and template sequences", value = 0),
                    numericInput("dnac_high", "High DNA concentration", value = 25, min = 0),
                    numericInput("dnac_low", "Low DNA concentration", value = 25, min = 0),
                    checkboxInput("self_comp", "Self-complementary", value = FALSE),
                    selectInput("nn_table", "Nearest Neighbor Table",
                                choices = list(
                                  "DNA_NN_SantaLucia_2004" = "DNA_NN_SantaLucia_2004",
                                  "DNA_NN_Breslauer_1986" = "DNA_NN_Breslauer_1986",
                                  "DNA_NN_Sugimoto_1996" = "DNA_NN_Sugimoto_1996",
                                  "DNA_NN_Allawi_1998" = "DNA_NN_Allawi_1998",
                                  "RNA_NN_Freier_1986" = "RNA_NN_Freier_1986",
                                  "RNA_NN_Xia_1998" = "RNA_NN_Xia_1998",
                                  "RNA_NN_Chen_2012" = "RNA_NN_Chen_2012",
                                  "RNA_DNA_NN_Sugimoto_1995" = "RNA_DNA_NN_Sugimoto_1995"
                                ),
                                selected = "DNA_NN_SantaLucia_2004"),
                    selectInput("tmm_table", "Terminal Mismatch Table",
                                choices = list(
                                  "DNA_TMM_Bommarito_2000" = "DNA_TMM_Bommarito_2000"
                                ),
                                selected = "DNA_TMM_Bommarito_2000"),
                    selectInput("imm_table", "Internal Mismatch Table",
                                choices = list(
                                  "DNA_IMM_Peyret_1999" = "DNA_IMM_Peyret_1999"
                                ),
                                selected = "DNA_IMM_Peyret_1999"),
                    selectInput("de_table", "Dangling Ends Table",
                                choices = list(
                                  "DNA_DE_Bommarito_2000" = "DNA_DE_Bommarito_2000",
                                  "RNA_DE_Turner_2010" = "RNA_DE_Turner_2010"
                                ),
                                selected = "DNA_DE_Bommarito_2000")
                  )
                )
              ),
              
              # Second row - Chemical and Salt Parameters
              fluidRow(
                conditionalPanel(
                  condition = "input.method == 'tm_nn' || input.method == 'tm_gc'",
                  box(
                    title = "Salt Corrections",
                    width = 6,
                    selectInput("salt_corr_method", "Salt Correction Method:",
                                choices = c("Schildkraut2010",
                                            "Wetmur1991",
                                            "SantaLucia1996",
                                            "SantaLucia1998-1",
                                            "Owczarzy2004",
                                            "Owczarzy2008"),
                                selected = "Schildkraut2010"),
                    numericInput("Na", "Na+ (mM)", value = 50, min = 0),
                    numericInput("K", "K+ (mM)", value = 0, min = 0),
                    numericInput("Tris", "Tris (mM)", value = 0, min = 0),
                    numericInput("Mg", "Mg2+ (mM)", value = 0, min = 0),
                    numericInput("dNTPs", "dNTPs (mM)", value = 0, min = 0)
                  ),
                  box(
                    title = "Chemical Corrections",
                    width = 6,
                    numericInput("DMSO", "DMSO (%)", value = 0, min = 0, max = 100),
                    numericInput("formamide", "Formamide", value = 0, min = 0),
                    selectInput("formamide_unit", "Formamide Unit",
                                choices = list("Percent (%)" = "percent", "Molar (M)" = "molar"),
                                selected = "percent"),
                    selectInput("formamide_factor", "Formamide Factor", 
                                choices = list("0.65" = "0.65", "0.6" = "0.6", "0.72" = "0.72"),
                                selected = "0.65"),
                    selectInput("DMSO_factor", "DMSO Factor",
                                choices = list("0.75" = "0.75", "0.5" = "0.5", "0.6" = "0.6", "0.65" = "0.65", "0.675" = "0.675"),
                                selected = "0.75")
                  )
                )
              ),
              
              # Third row - Results
              fluidRow(
                box(
                  title = "Results",
                  width = 12,
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
  # Update sequence selectors when FASTA file is uploaded
  
  # Calculate Tm when button is clicked
  observeEvent(input$calculate, {
    tryCatch({
      # Get sequences based on input type
      if (input$input_type == "direct") {
        primers <- unlist(strsplit(input$input_seq, ","))
        templates <- if(input$method == "tm_nn" && !is.null(input$rev_input_seq) && input$rev_input_seq != "") {
          unlist(strsplit(input$rev_input_seq, ","))
        } else {
          NULL
        }
      } else {
        req(input$fasta_file)
        primers <- input$fasta_file$datapath
        templates <- if(input$method == "tm_nn" && !is.null(input$fasta_file_complement)) {
          req(input$fasta_file_complement)
          input$fasta_file_complement$datapath
        } else {
          NULL
        }
      }
      
      # Calculate Tm using tm_calculate
      result <- tm_calculate(
        input_seq = primers,
        complement_seq = templates,
        method = input$method,
        ambiguous = input$ambiguous,
        shift = input$shift,
        nn_table = input$nn_table,
        tmm_table = input$tmm_table,
        imm_table = input$imm_table,
        de_table = input$de_table,
        dnac_high = input$dnac_high,
        dnac_low = input$dnac_low,
        self_comp = input$self_comp,
        variant = input$variant,
        Na = input$Na,
        K = input$K,
        Tris = input$Tris,
        Mg = input$Mg,
        dNTPs = input$dNTPs,
        salt_corr_method = input$salt_corr_method,
        DMSO = input$DMSO,
        formamide_value_unit = list(
          value = input$formamide,
          unit = input$formamide_unit
        ),
        dmso_factor = as.numeric(input$DMSO_factor),
        formamide_factor = as.numeric(input$formamide_factor)
      )
      
      
      # Display results
      output$results <- renderPrint({
        print(result)
      })
      
      # Enable download button
      output$download_results <- downloadHandler(
        filename = function() {
          paste("tm_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt", sep = "")
        },
        content = function(file) {
          writeLines(capture.output(
            cat("Method:", input$method, "\n\n"),
            for (i in seq_along(results)) {
              cat("Sequence:", names(results)[i], "\n")
              if (input$input_type == "fasta") {
                cat("Sequence: ", primers[[i]], "\n")
                if (!is.null(templates)) {
                  cat("Complement: ", templates[[i]], "\n")
                }
              }
              cat("\nResults:\n")
              print(results[[i]])
              cat("\n", paste(rep("-", 50), collapse = ""), "\n\n")
            }
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