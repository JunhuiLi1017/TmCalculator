#' Sequence Tm module UI
#'
#' Four input modes are merged into a single panel:
#'   1. paste     - raw sequences (one per line, optional >name headers)
#'   2. coords    - genomic coordinate strings
#'   3. fasta     - uploaded FASTA file
#'   4. bsgenome  - genome-wide tiling of an installed BSgenome package
#'
#' @noRd
mod_sequence_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      title = "Input & parameters",

      ## ---- Input mode ------------------------------------------------------
      shiny::radioButtons(
        ns("input_mode"), "Input type",
        choices = c(
          "Paste sequences"     = "paste",
          "Genomic coordinates" = "coords",
          "Upload FASTA"        = "fasta",
          "BSgenome (genome-wide)" = "bsgenome"
        ),
        selected = "paste"
      ),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'paste'", ns("input_mode")),
        shiny::textAreaInput(
          ns("paste_seq"), "Sequences",
          rows = 5,
          placeholder = ">probe1\nATCGTGCGTAGCAGTACGATCAGTAG\nGCGCATATGCGC"
        )
      ),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'coords'", ns("input_mode")),
        shiny::textAreaInput(
          ns("coord_seq"), "Coordinate strings",
          rows = 4,
          placeholder = "chr1:1000-1200:+:BSgenome.Hsapiens.UCSC.hg38:region1"
        ),
        shiny::helpText(
          "One per line, fields separated by ':' in the order ",
          shiny::tags$code("chr:start-end:strand:pkg_name:region_id"),
          ". region_id is optional. The BSgenome package named in pkg_name ",
          "must be installed."
        )
      ),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'fasta'", ns("input_mode")),
        shiny::fileInput(
          ns("fasta_file"), "FASTA file",
          accept = c(".fa", ".fasta", ".fna", ".txt")
        )
      ),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'bsgenome'", ns("input_mode")),
        shiny::uiOutput(ns("genome_picker")),
        shiny::textInput(
          ns("chromosomes"), "Chromosome(s)",
          placeholder = "chr1  (blank = all; comma-separated for several)"
        ),
        shiny::numericInput(ns("window"), "Window (bp)", value = 200, min = 1),
        shiny::numericInput(ns("slide"), "Slide / step (bp)", value = 50, min = 1),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput(ns("start"), "Start (opt.)",
                                               value = NA)),
          shiny::column(6, shiny::numericInput(ns("end"), "End (opt.)",
                                               value = NA))
        ),
        shiny::selectInput(ns("strand"), "Strand",
                           choices = c("+", "-"), selected = "+"),
        shiny::selectInput(
          ns("trim_N"), "N handling",
          choices = c("Trim ends" = "ends", "Filter windows" = "filter",
                      "Keep (none)" = "none"),
          selected = "ends"
        ),
        shiny::numericInput(ns("max_N_frac"), "Max N fraction",
                            value = 0.10, min = 0, max = 1, step = 0.05),
        shiny::helpText(
          "Genome-wide tiling can produce many windows; restrict the ",
          "chromosome or start/end range to keep runtimes reasonable."
        )
      ),

      shiny::hr(),

      ## ---- Method ----------------------------------------------------------
      shiny::selectInput(
        ns("method"), "Method",
        choices = c("Nearest neighbour" = "tm_nn",
                    "GC content"        = "tm_gc",
                    "Wallace rule"      = "tm_wallace"),
        selected = "tm_nn"
      ),

      ## NN options
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'tm_nn'", ns("method")),
        shiny::selectInput(ns("nn_table"), "NN table",
                           choices = .tcs_nn_tables(),
                           selected = "DNA_NN_SantaLucia_2004"),
        shiny::numericInput(ns("dnac_high"), "dnac_high (nM)", value = 25, min = 0),
        shiny::numericInput(ns("dnac_low"), "dnac_low (nM)", value = 25, min = 0),
        shiny::checkboxInput(ns("self_comp"), "Self-complementary", FALSE)
      ),

      ## GC options (incl. custom A/B/C/D)
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'tm_gc'", ns("method")),
        shiny::checkboxInput(ns("gc_custom"),
                             "Custom coefficients (A, B, C, D)", FALSE),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == false", ns("gc_custom")),
          shiny::selectInput(ns("gc_variant"), "Variant",
                             choices = .tcs_gc_variants(),
                             selected = "Primer3Plus")
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("gc_custom")),
          shiny::helpText(
            "Tm = A + B * %GC - C / N - D * %mismatch"
          ),
          shiny::fluidRow(
            shiny::column(6, shiny::numericInput(ns("gc_A"), "A (intercept)",
                                                 value = 81.5)),
            shiny::column(6, shiny::numericInput(ns("gc_B"), "B (x %GC)",
                                                 value = 0.41))
          ),
          shiny::fluidRow(
            shiny::column(6, shiny::numericInput(ns("gc_C"), "C (/ N)",
                                                 value = 600)),
            shiny::column(6, shiny::numericInput(ns("gc_D"), "D (mismatch)",
                                                 value = 1))
          ),
          shiny::selectInput(ns("gc_salt"), "Salt correction",
                             choices = .tcs_salt_methods(),
                             selected = "Schildkraut2010")
        ),
        shiny::checkboxInput(ns("gc_mismatch"), "Count 'X' as mismatch", TRUE)
      ),

      ## Ion / salt conditions (shared by nn & gc)
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] != 'tm_wallace'", ns("method")),
        shiny::hr(),
        shiny::tags$b("Ion concentrations (mM)"),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput(ns("Na"), "Na+", value = 50, min = 0)),
          shiny::column(6, shiny::numericInput(ns("K"), "K+", value = 0, min = 0))
        ),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput(ns("Tris"), "Tris", value = 0, min = 0)),
          shiny::column(6, shiny::numericInput(ns("Mg"), "Mg2+", value = 0, min = 0))
        ),
        shiny::numericInput(ns("dNTPs"), "dNTPs", value = 0, min = 0),
        shiny::checkboxInput(ns("ambiguous"), "Account for ambiguous bases", FALSE)
      ),

      shiny::hr(),
      shiny::textInput(ns("label"), "Result label",
                       placeholder = "auto"),
      shiny::actionButton(ns("run"), "Calculate Tm",
                          class = "btn-primary", width = "100%")
    ),

    ## ---- Main output -------------------------------------------------------
    bslib::card(
      bslib::card_header("Results"),
      shiny::uiOutput(ns("status")),
      DT::DTOutput(ns("table"))
    )
  )
}

#' Sequence Tm module server
#'
#' @param id module id.
#' @param add_result function(label, gr, options) registering a result in the
#'   shared store; returns the (possibly de-duplicated) label used.
#' @noRd
mod_sequence_server <- function(id, add_result) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$genome_picker <- shiny::renderUI({
      genomes <- .tcs_installed_genomes()
      if (length(genomes) == 0) {
        shiny::helpText(
          "No BSgenome.* packages are installed. Install one, e.g. ",
          shiny::tags$code(
            "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
          ), "."
        )
      } else {
        shiny::selectInput(ns("genome"), "BSgenome package",
                           choices = genomes, selected = genomes[1])
      }
    })

    rv <- shiny::reactiveValues(df = NULL, msg = NULL, ok = FALSE)

    # ---- Build the input the engine expects --------------------------------
    build_input <- function() {
      mode <- input$input_mode

      if (mode == "paste") {
        seqs <- .tcs_parse_sequence_box(input$paste_seq)
        if (is.null(seqs)) stop("Please paste at least one sequence.")
        return(list(kind = "vector", value = unname(seqs),
                    region_ids = names(seqs)))
      }

      if (mode == "coords") {
        txt <- input$coord_seq
        if (is.null(txt) || !nzchar(trimws(txt)))
          stop("Please enter at least one coordinate string.")
        coords <- trimws(strsplit(txt, "\r?\n")[[1]])
        coords <- coords[nzchar(coords)]
        return(list(kind = "vector", value = coords, region_ids = NULL))
      }

      if (mode == "fasta") {
        f <- input$fasta_file
        if (is.null(f)) stop("Please upload a FASTA file.")
        return(list(kind = "file", value = f$datapath, region_ids = NULL))
      }

      if (mode == "bsgenome") {
        if (is.null(input$genome) || !nzchar(input$genome))
          stop("Please select an installed BSgenome package.")
        if (!.tcs_has_pkg(input$genome))
          stop(sprintf("Package '%s' is not installed.", input$genome))
        bsg <- getExportedValue(input$genome, input$genome)

        chrom <- NULL
        if (!is.null(input$chromosomes) && nzchar(trimws(input$chromosomes))) {
          chrom <- trimws(strsplit(input$chromosomes, ",")[[1]])
          chrom <- chrom[nzchar(chrom)]
        }
        start <- if (is.null(input$start) || is.na(input$start)) NULL else input$start
        end <- if (is.null(input$end) || is.na(input$end)) NULL else input$end

        coords <- TmCalculator::make_genomiccoord(
          bsgenome    = bsg,
          chromosomes = chrom,
          window      = as.integer(input$window),
          slide       = as.integer(input$slide),
          start       = start,
          end         = end,
          strand      = input$strand,
          trim_N      = input$trim_N,
          max_N_frac  = input$max_N_frac,
          genome_pkg_name = input$genome,
          as_vector   = TRUE,
          verbose     = FALSE
        )
        if (length(coords) == 0)
          stop("No windows produced. Relax the N filter or widen the range.")
        return(list(kind = "vector", value = coords, region_ids = NULL))
      }

      stop("Unknown input mode.")
    }

    # ---- Assemble the tm_calculate argument list ---------------------------
    build_args <- function(input_value) {
      args <- list(input_seq = input_value,
                   method = input$method,
                   ambiguous = isTRUE(input$ambiguous))

      if (input$method == "tm_nn") {
        args$nn_table <- input$nn_table
        args$dnac_high <- input$dnac_high
        args$dnac_low <- input$dnac_low
        args$self_comp <- isTRUE(input$self_comp)
      }

      if (input$method == "tm_gc") {
        args$mismatch <- isTRUE(input$gc_mismatch)
        if (isTRUE(input$gc_custom)) {
          args$userset <- c(input$gc_A, input$gc_B, input$gc_C, input$gc_D)
          args$salt_method <- input$gc_salt
        } else {
          args$variant <- input$gc_variant
        }
      }

      if (input$method != "tm_wallace") {
        args$Na <- input$Na
        args$K <- input$K
        args$Tris <- input$Tris
        args$Mg <- input$Mg
        args$dNTPs <- input$dNTPs
      }
      args
    }

    shiny::observeEvent(input$run, {
      rv$ok <- FALSE
      res <- tryCatch({
        spec <- build_input()
        args <- build_args(spec$value)
        out <- do.call(TmCalculator::tm_calculate, args)
        gr <- out$gr
        if (is.null(gr)) gr <- TmCalculator::to_genomic_ranges(spec$value)
        # Attach region ids when we have them (paste mode with names)
        if (!is.null(spec$region_ids) &&
            length(spec$region_ids) == length(gr) &&
            !"region_id" %in% names(GenomicRanges::mcols(gr))) {
          GenomicRanges::mcols(gr)$region_id <- spec$region_ids
        }
        opts <- out$options
        list(gr = gr, options = opts)
      }, error = function(e) {
        rv$msg <- conditionMessage(e)
        NULL
      })

      if (is.null(res)) {
        rv$df <- NULL
        return(invisible(NULL))
      }

      label <- if (nzchar(trimws(input$label %||% ""))) trimws(input$label) else
        paste0(input$method, "_", format(Sys.time(), "%H%M%S"))
      used_label <- add_result(label, res$gr, res$options)

      rv$df <- as.data.frame(res$gr)
      rv$msg <- sprintf("Saved %d region(s) as \"%s\".",
                        length(res$gr), used_label)
      rv$ok <- TRUE
    })

    output$status <- shiny::renderUI({
      if (is.null(rv$msg)) return(NULL)
      cls <- if (isTRUE(rv$ok)) "alert alert-success" else "alert alert-danger"
      shiny::tags$div(class = cls, rv$msg)
    })

    output$table <- DT::renderDT({
      shiny::req(rv$df)
      DT::datatable(rv$df, options = list(pageLength = 10, scrollX = TRUE),
                    rownames = FALSE)
    })
  })
}
