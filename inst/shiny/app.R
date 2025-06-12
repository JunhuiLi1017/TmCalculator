library(shiny)
library(TmCalculator)
library(shinydashboard)
library(GenomicRanges) # Required for GRanges object manipulation

# Ensure the package is properly loaded
if (!requireNamespace("TmCalculator", quietly = TRUE)) {
  stop("TmCalculator package is required but not installed")
}
# Ensure karyoploteR is available if you're using plot_tm_karyotype
if (!requireNamespace("karyoploteR", quietly = TRUE)) {
  stop("karyoploteR package is required for karyotype plots but not installed")
}
# Ensure ggbio is available if you're using plot_tm_heatmap or plot_tm_genome_tracks
if (!requireNamespace("ggbio", quietly = TRUE)) {
  stop("ggbio package is required for heatmap and genome tracks plots but not installed")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 package is required for ggbio plots but not installed")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  stop("RColorBrewer package is required for ggbio plots but not installed")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
  stop("viridis package is required for color palettes but not installed")
}


# --- UI Definition ---
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
                                 "FASTA File" = "fasta",
                                 "Genomic Coordinates" = "genomic_coordinates"
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
                              placeholder = "for example AGTCTGTT or AGTCTGTT, CATGCACG")
                  ),
                  conditionalPanel(
                    condition = "input.input_type == 'fasta' && input.method == 'tm_nn'",
                    fileInput("fasta_file", "Upload FASTA File",
                              accept = c(".fasta", ".fa", ".txt")),
                    fileInput("fasta_file_complement", "Upload reverse complement FASTA File",
                              accept = c(".fasta", ".fa", ".txt"))
                  ),
                  conditionalPanel(
                    condition = "input.input_type == 'genomic_coordinates' && input.method == 'tm_nn'",
                    textInput("genomic_coordinates", "Genomic Coordinates",
                              placeholder = "for example chr1:100-200:+:hg38 or chr2:100-200:-:hg38"),
                    textInput("genomic_coordinates_complement", "Genomic Coordinates of reverse complement",
                              placeholder = "for example chr1:100-200:-:hg38 or chr2:100-200:+:hg38")
                  ),
                  conditionalPanel(
                    condition = "input.input_type == 'genomic_coordinates' && (input.method == 'tm_gc' || input.method == 'tm_wallace')",
                    textInput("genomic_coordinates", "Genomic Coordinates",
                              placeholder = "for example chr1:100-200:+:hg38 or chr2:100-200:-:hg38")
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
                ),
                box(
                  title = "Visualization",
                  width = 12,
                  fluidRow(
                    column(3,
                      selectInput("plot_type", "Plot Type",
                        choices = list(
                          "Karyotype Plot" = "karyotype",
                          "Heatmap Plot" = "heatmap",
                          "Genome Tracks" = "genome_tracks"
                        ),
                        selected = "karyotype"
                      )
                    ),
                    column(9, # Increased column width for plot options
                      # --- Karyotype Plot Options (Updated from previous turn) ---
                      conditionalPanel(
                        condition = "input.plot_type == 'karyotype'",
                        textInput("karyotype_genome_assembly", "Genome Assembly",
                          placeholder = "e.g., hg38, hg19, mm10, etc."
                        ),
                        # Chromosomes will typically be populated dynamically from your data on the server side
                        selectInput("karyotype_chromosomes", "Chromosomes to Plot",
                          choices = NULL, # Choices will be set by the server
                          multiple = TRUE
                        ),
                        textInput("karyotype_colors", "Chromosome Colors (named, e.g., chr1=red,chr2=blue)",
                          placeholder = "e.g., chr1=#FF0000,chr2=blue,chrX=#00FF00"
                        ),
                        textInput("karyotype_shapes", "Chromosome Shapes (named, e.g., chr1=16,chr2=17)",
                          placeholder = "e.g., chr1=16,chr2=17,chrY=18 (see ggplot2 pch values 0-25, 32-127)"
                        ),
                        selectInput("karyotype_plot_type", "Karyotype Layout (plot_type)",
                          choices = as.list(setNames(1:7, 1:7)), # Options 1 through 7 for plotKaryotype layout
                          selected = "1"
                        ),
                        numericInput("karyotype_point_cex", "Point Size (cex)", value = 1.5, min = 0.1, max = 3, step = 0.1),
                        numericInput("karyotype_xaxis_cex", "X-axis Text Size (cex)", value = 0.7, min = 0.1, max = 1.5, step = 0.1),
                        numericInput("karyotype_yaxis_cex", "Y-axis Text Size (cex)", value = 0.8, min = 0.1, max = 1.5, step = 0.1),
                        numericInput("karyotype_chr_cex", "Chromosome Name Text Size (chr_cex)", value = 1, min = 0.1, max = 2, step = 0.1),
                        numericInput("karyotype_tick_dist", "X-axis Tick Distance (bp)", value = 10000000, min = 100000, max = 100000000, step = 1000000)
                      ),
                      # --- Heatmap Plot Options ---
                      conditionalPanel(
                        condition = "input.plot_type == 'heatmap'",
                        textInput("heatmap_genome_assembly", "Genome Assembly",
                          placeholder = "e.g., hg38, hg19, mm10, etc."
                        ),
                        selectInput("heatmap_chromosomes", "Chromosome to Plot",
                          choices = NULL, # Choices will be set by the server
                          multiple = TRUE  # Multiple chromosome selection
                        ),
                        textInput("heatmap_zoom", "Zoom Region",
                          placeholder = "e.g., chr1:1000000-2000000, chr2:1000000-2000000"
                        ),
                        selectInput("heatmap_plot_type", "Heatmap Type",
                          choices = list(
                            "Karyogram" = "karyogram",
                            "Faceted" = "faceted"
                          ),
                          selected = "karyogram"
                        ),
                        selectInput("heatmap_color_palette", "Color Palette",
                          choices = list(
                            "Viridis" = "viridis",
                            "Magma" = "magma",
                            "Plasma" = "plasma",
                            "Inferno" = "inferno",
                            "Cividis" = "cividis"
                          ),
                          selected = "viridis"
                        ),
                        textInput("heatmap_title", "Plot Title", value = "Tm Values on Chromosomes (ggbio)",
                                  placeholder = "e.g., My Tm Heatmap")
                      ),
                      # --- Genome Tracks Plot Options ---
                      conditionalPanel(
                        condition = "input.plot_type == 'genome_tracks'",
                        textInput("genome_tracks_genome_assembly", "Genome Assembly",
                          placeholder = "e.g., hg38, hg19, mm10, etc."
                        ),
                        selectInput("genome_tracks_chromosomes", "Chromosome to Plot",
                          choices = NULL, # Choices will be set by the server
                          multiple = FALSE  # Single chromosome selection
                        ),
                        textInput("genome_tracks_zoom", "Zoom Region",
                          placeholder = "e.g., chr1:1000000-2000000"
                        ),
                        selectInput("genome_tracks_color_palette", "Color Palette",
                          choices = list(
                            "Viridis" = "viridis",
                            "Magma" = "magma",
                            "Plasma" = "plasma",
                            "Inferno" = "inferno",
                            "Cividis" = "cividis"
                          ),
                          selected = "viridis"
                        ),
                        checkboxInput("genome_tracks_show_ideogram", "Show Ideogram", value = TRUE)
                      )
                    )
                  ),
                  plotOutput("plot_output"),
                  downloadButton("download_plot", "Download Plot")
                )
              )
      ),

      # About Tab
      tabItem(tabName = "about",
              box(
                title = "About Tm Calculator",
                width = 12,
                # Make sure about.md exists in inst/shiny/about.md within your package
                includeMarkdown(system.file("shiny", "about.md", package = "TmCalculator"))
              )
      )
    )
  )
)

# --- Server Definition ---
server <- function(input, output, session) {

  # Initialize result to NULL globally within the server session
  # This makes 'result' accessible to observers outside of observeEvent
  result <- reactiveVal(NULL)

  # Helper function to parse named string inputs (e.g., "chr1=red,chr2=blue")
  parse_named_vector <- function(text_input_string, type = "character") {
    if (is.null(text_input_string) || text_input_string == "") {
      return(NULL)
    }

    pairs <- strsplit(text_input_string, ",")[[1]]
    names_vec <- character(length(pairs))
    values_vec <- character(length(pairs))

    for (i in seq_along(pairs)) {
      parts <- strsplit(pairs[i], "=")[[1]]
      if (length(parts) == 2) {
        names_vec[i] <- trimws(parts[1])
        values_vec[i] <- trimws(parts[2])
      } else {
        # Warning for malformed input
        warning(paste("Invalid format for pair:", pairs[i], ". Skipping this entry."))
      }
    }

    # Filter out any pairs that were malformed (i.e., names_vec[i] is empty)
    valid_indices <- nchar(names_vec) > 0
    names_vec <- names_vec[valid_indices]
    values_vec <- values_vec[valid_indices]

    if (type == "numeric") {
      values_vec <- as.numeric(values_vec)
    }

    if (length(names_vec) == 0) return(NULL) # Return NULL if no valid pairs were found
    setNames(values_vec, names_vec)
  }

  # Update chromosome choices for karyotype, genome tracks, and heatmap plots
  observe({
    req(result()) # Ensure result is not NULL
    current_gr <- result()$tm$Tm
    if (!is.null(current_gr) && inherits(current_gr, "GRanges")) {
      chromosomes <- unique(as.character(GenomicRanges::seqnames(current_gr)))
      # Update for karyotype
      updateSelectInput(session, "karyotype_chromosomes",
        choices = chromosomes,
        selected = chromosomes # Select all chromosome by default
      )
      # Update for genome tracks
      updateSelectInput(session, "genome_tracks_chromosomes",
        choices = chromosomes,
        selected = chromosomes[1] # Select first chromosome by default
      )
      # Update for heatmap
      updateSelectInput(session, "heatmap_chromosomes",
        choices = chromosomes,
        selected = chromosomes # Select all chromosomes by default
      )
    }
  })

  # Calculate Tm when button is clicked
  observeEvent(input$calculate, {
    # Create progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Calculating Tm values", value = 0)
    on.exit(progress$close())

    tryCatch({
      # Get sequences based on input type
      primers <- NULL
      templates <- NULL

      if (input$input_type == "direct") {
        primers <- unlist(strsplit(input$input_seq, ","))
        if(input$method == "tm_nn" && !is.null(input$rev_input_seq) && input$rev_input_seq != "") {
          templates <- unlist(strsplit(input$rev_input_seq, ","))
        }
      } else if (input$input_type == "genomic_coordinates") {
        primers <- input$genomic_coordinates
        if(input$method == "tm_nn" && !is.null(input$genomic_coordinates_complement) && input$genomic_coordinates_complement != "") {
          templates <- input$genomic_coordinates_complement
        }
      } else if (input$input_type == "fasta") {
        req(input$fasta_file)
        primers <- input$fasta_file$datapath
        if(input$method == "tm_nn" && !is.null(input$fasta_file_complement)) {
          req(input$fasta_file_complement)
          templates <- input$fasta_file_complement$datapath
        }
      }

      # Update progress
      progress$set(message = "Processing sequences", value = 0.2)

      # Calculate Tm using tm_calculate
      calculated_result <- tm_calculate(
        input_seq = primers,
        complement_seq = templates,
        method = input$method,
        ambiguous = input$ambiguous,
        # Only pass NN-specific parameters if method is tm_nn
        shift = if(input$method == "tm_nn") input$shift else NULL,
        nn_table = if(input$method == "tm_nn") input$nn_table else NULL,
        tmm_table = if(input$method == "tm_nn") input$tmm_table else NULL,
        imm_table = if(input$method == "tm_nn") input$imm_table else NULL,
        de_table = if(input$method == "tm_nn") input$de_table else NULL,
        dnac_high = if(input$method == "tm_nn") input$dnac_high else NULL,
        dnac_low = if(input$method == "tm_nn") input$dnac_low else NULL,
        self_comp = if(input$method == "tm_nn") input$self_comp else FALSE,
        # Only pass GC-specific parameters if method is tm_gc
        variant = if(input$method == "tm_gc") input$variant else NULL,
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

      # Update progress
      progress$set(message = "Preparing results", value = 0.8)

      # Store result in reactiveVal
      result(calculated_result)

      # Display results
      output$results <- renderPrint({
        print(result())
      })

      # Create plot based on selected type
      output$plot_output <- renderPlot({
        req(result()) # Ensure result is available before plotting
        current_gr <- result()$tm$Tm
        req(current_gr) # Ensure the GRanges object is not NULL

        # Create progress object for plotting
        plot_progress <- shiny::Progress$new()
        plot_progress$set(message = "Generating plot", value = 0)
        on.exit(plot_progress$close())

        # Parse colors and shapes for karyotype plot
        karyotype_colors_parsed <- parse_named_vector(input$karyotype_colors, type = "character")
        karyotype_shapes_parsed <- parse_named_vector(input$karyotype_shapes, type = "numeric")

        # Determine genome assembly for plotting functions from input
        plot_genome_assembly_karyotype <- if (input$karyotype_genome_assembly != "") input$karyotype_genome_assembly else NULL
        plot_genome_assembly_heatmap <- if (input$heatmap_genome_assembly != "") input$heatmap_genome_assembly else NULL
        plot_genome_assembly_genome_tracks <- if (input$genome_tracks_genome_assembly != "") input$genome_tracks_genome_assembly else NULL

        # Update plot progress
        plot_progress$set(message = "Rendering plot", value = 0.5)

        p <- switch(input$plot_type,
          "karyotype" = {
            plot_tm_karyotype(
              gr = current_gr,
              chromosomes = input$karyotype_chromosomes,
              genome_assembly = plot_genome_assembly_karyotype,
              colors = karyotype_colors_parsed,
              shapes = karyotype_shapes_parsed,
              plot_type = as.numeric(input$karyotype_plot_type),
              point_cex = input$karyotype_point_cex,
              xaxis_cex = input$karyotype_xaxis_cex,
              yaxis_cex = input$karyotype_yaxis_cex,
              chr_cex = input$karyotype_chr_cex,
              tick_dist = input$karyotype_tick_dist
            )
          },
          "heatmap" = {
            plot_tm_heatmap(
              gr = current_gr,
              genome_assembly = plot_genome_assembly_heatmap,
              chromosome_to_plot = input$heatmap_chromosomes,
              plot_type = input$heatmap_plot_type,
              color_palette = input$heatmap_color_palette,
              title_name = input$heatmap_title,
              zoom = if (input$heatmap_zoom != "") input$heatmap_zoom else NULL
            )
          },
          "genome_tracks" = {
            plot_tm_genome_tracks(
              gr = current_gr,
              chromosome_to_plot = input$genome_tracks_chromosomes,
              genome_assembly = plot_genome_assembly_genome_tracks,
              color_palette = input$genome_tracks_color_palette,
              show_ideogram = input$genome_tracks_show_ideogram,
              zoom = if (input$genome_tracks_zoom != "") input$genome_tracks_zoom else NULL
            )
          }
        )

        # Update plot progress
        plot_progress$set(message = "Finalizing plot", value = 0.9)

        return(p)
      })

      # Update progress
      progress$set(message = "Calculation complete", value = 1)

    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
