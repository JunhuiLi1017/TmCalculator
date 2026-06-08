#' Visualize module UI
#'
#' Two views:
#'   - "Tm plots": linear / heatmap / karyotype views of a stored result.
#'   - "Multi-omics tracks": integrate uploaded feature ranges with a stored
#'     Tm result and render multi-track linear / circular / karyotype genome
#'     plots.
#' @noRd
mod_visualize_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      title = "Plot settings",

      shiny::selectInput(ns("result"), "Stored Tm result", choices = NULL),

      shiny::radioButtons(
        ns("view"), "View",
        choices = c("Tm plots" = "tm",
                    "Multi-omics tracks" = "multi"),
        selected = "tm"
      ),

      ## ---- Tm plots --------------------------------------------------------
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'tm'", ns("view")),
        shiny::selectInput(
          ns("tm_plot"), "Plot type",
          choices = c("Linear (Tm vs position)" = "linear",
                      "Heatmap (interactive)"   = "heatmap",
                      "Karyotype (interactive)" = "karyo")
        ),
        shiny::selectInput(ns("palette"), "Palette",
                           choices = c("viridis", "magma", "plasma",
                                       "inferno", "cividis"),
                           selected = "viridis"),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'linear'", ns("tm_plot")),
          shiny::selectInput(ns("lin_x"), "X axis",
                             choices = c("index", "label", "position"),
                             selected = "index"),
          shiny::checkboxInput(ns("lin_line"), "Connect points with line", FALSE),
          shiny::checkboxInput(ns("lin_facet"), "Facet by chromosome", FALSE)
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] != 'linear'", ns("tm_plot")),
          shiny::textInput(ns("assembly"), "Genome assembly (e.g. hg38)",
                           value = "")
        )
      ),

      ## ---- Multi-omics tracks ---------------------------------------------
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'multi'", ns("view")),
        shiny::fileInput(ns("feat_file"), "Feature ranges (BED / CSV)",
                         accept = c(".bed", ".csv", ".tsv", ".txt")),
        shiny::selectInput(ns("feat_col"), "Feature value column", choices = NULL),
        shiny::selectInput(
          ns("engine"), "Track engine",
          choices = c("Linear genome"  = "linear",
                      "Circular (circos)" = "circos",
                      "Karyotype"       = "karyo")
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] != 'circos'", ns("engine")),
          shiny::textInput(ns("m_assembly"), "Genome assembly (e.g. hg38)",
                           value = ""),
          shiny::textInput(ns("m_chroms"), "Chromosomes (comma-separated)",
                           value = "")
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'circos'", ns("engine")),
          shiny::textInput(ns("m_genome_name"), "Genome name", value = "genome"),
          shiny::numericInput(ns("m_genome_size"),
                              "Genome size (bp, single chromosome)",
                              value = NA)
        ),
        shiny::hr(),
        shiny::tags$b("Annotate Tm ranges (optional)"),
        shiny::selectInput(ns("strategy"), "Integration strategy",
                           choices = c("overlap", "nearest", "window", "bin"),
                           selected = "overlap"),
        shiny::actionButton(ns("integrate"), "Integrate & store result",
                            class = "btn-outline-secondary btn-sm",
                            width = "100%"),
        shiny::helpText(
          "Integration runs integrate_granges() and saves an annotated result ",
          "that also becomes available in Compare groups."
        )
      ),

      shiny::hr(),
      shiny::actionButton(ns("draw"), "Draw", class = "btn-primary",
                          width = "100%")
    ),

    bslib::card(
      bslib::card_header("Plot"),
      shiny::uiOutput(ns("status")),
      shiny::uiOutput(ns("plot_area"))
    )
  )
}

#' Visualize module server
#' @param id module id.
#' @param store shared reactiveValues with `$results`. `store` must also expose
#'   nothing else; integration registers results via the same list.
#' @noRd
mod_visualize_server <- function(id, store) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    result_names <- shiny::reactive(names(store$results))
    shiny::observe({
      shiny::updateSelectInput(session, "result", choices = result_names())
    })

    current_gr <- shiny::reactive({
      shiny::req(input$result)
      .tcs_result_gr(store$results[[input$result]])
    })

    # Uploaded feature ranges
    features <- shiny::reactive({
      f <- input$feat_file
      shiny::req(f)
      .tcs_read_features(f$datapath)
    })

    shiny::observeEvent(features(), {
      gr <- features()
      cols <- names(GenomicRanges::mcols(gr))
      num_cols <- cols[vapply(cols, function(c)
        is.numeric(GenomicRanges::mcols(gr)[[c]]), logical(1))]
      shiny::updateSelectInput(session, "feat_col",
                               choices = if (length(num_cols)) num_cols else cols,
                               selected = if (length(num_cols)) num_cols[1]
                               else cols[1])
    })

    rv <- shiny::reactiveValues(kind = NULL, obj = NULL, msg = NULL, ok = FALSE)

    .chrom_vec <- function(txt) {
      if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
      v <- trimws(strsplit(txt, ",")[[1]]); v[nzchar(v)]
    }

    # Build a track_list from a Tm GRanges + (optional) feature GRanges
    build_tracks <- function(gr_tm, gr_feat, feat_col, type_tm = "points") {
      tl <- list(
        list(type = type_tm, data = gr_tm, value_col = "Tm",
             col = "#e41a1c", name = "Tm")
      )
      if ("GC" %in% names(GenomicRanges::mcols(gr_tm))) {
        tl[[length(tl) + 1L]] <- list(type = "line", data = gr_tm,
                                      value_col = "GC", col = "#377eb8",
                                      name = "GC")
      }
      if (!is.null(gr_feat) && !is.null(feat_col) &&
          feat_col %in% names(GenomicRanges::mcols(gr_feat))) {
        tl[[length(tl) + 1L]] <- list(type = "line", data = gr_feat,
                                      value_col = feat_col, col = "#4daf4a",
                                      name = feat_col)
      }
      tl
    }

    # ---- Integrate & store -------------------------------------------------
    shiny::observeEvent(input$integrate, {
      out <- tryCatch({
        gr_tm <- current_gr()
        if (is.null(gr_tm)) stop("Select a stored Tm result first.")
        gr_feat <- features()
        ann <- TmCalculator::integrate_granges(gr_tm, gr_feat,
                                               strategy = input$strategy)
        ann
      }, error = function(e) { rv$msg <- conditionMessage(e); rv$ok <- FALSE; NULL })

      if (!is.null(out)) {
        label <- paste0(input$result, "+features(", input$strategy, ")")
        store$results[[label]] <- list(gr = out, options = NULL)
        rv$msg <- sprintf("Stored integrated result \"%s\" (%d ranges).",
                          label, length(out))
        rv$ok <- TRUE
      }
    })

    # ---- Draw --------------------------------------------------------------
    shiny::observeEvent(input$draw, {
      rv$ok <- FALSE
      res <- tryCatch({
        gr <- current_gr()
        if (is.null(gr)) stop("Select a stored result first.")

        if (input$view == "tm") {
          if (input$tm_plot == "linear") {
            p <- TmCalculator::plot_tm_linear(
              gr, x_axis = input$lin_x, color_palette = input$palette,
              add_line = isTRUE(input$lin_line),
              facet_by_chr = isTRUE(input$lin_facet))
            list(kind = "ggplot", obj = p)
          } else if (input$tm_plot == "heatmap") {
            p <- TmCalculator::plot_tm_heatmap_interactive(
              gr,
              genome_assembly = if (nzchar(input$assembly)) input$assembly else NULL,
              color_palette = input$palette)
            list(kind = "plotly", obj = p)
          } else {
            p <- TmCalculator::plot_tm_karyotype_interactive(
              gr,
              genome_assembly = if (nzchar(input$assembly)) input$assembly else NULL)
            list(kind = "plotly", obj = p)
          }
        } else {
          gr_feat <- features()
          tl <- build_tracks(gr, gr_feat, input$feat_col)
          if (input$engine == "linear") {
            list(kind = "base", obj = function() {
              TmCalculator::plot_linear_genome(
                genome_name = if (nzchar(input$m_assembly)) input$m_assembly else "genome",
                track_list = tl,
                chromosomes = .chrom_vec(input$m_chroms))
            })
          } else if (input$engine == "circos") {
            gsize <- input$m_genome_size
            if (is.null(gsize) || is.na(gsize))
              gsize <- max(GenomicRanges::end(gr))
            list(kind = "base", obj = function() {
              TmCalculator::plot_circos_genome(
                genome_name = input$m_genome_name %||% "genome",
                genome_size = gsize,
                track_list = tl)
            })
          } else {
            list(kind = "base", obj = function() {
              TmCalculator::plot_karyotype_genome(
                genome_assembly = if (nzchar(input$m_assembly)) input$m_assembly else "hg38",
                track_list = tl,
                chromosomes = .chrom_vec(input$m_chroms) %||% "auto")
            })
          }
        }
      }, error = function(e) { rv$msg <- conditionMessage(e); rv$ok <- FALSE; NULL })

      if (is.null(res)) { rv$kind <- NULL; rv$obj <- NULL; return(invisible(NULL)) }
      rv$kind <- res$kind
      rv$obj <- res$obj
      rv$msg <- NULL
      rv$ok <- TRUE
    })

    output$status <- shiny::renderUI({
      if (length(result_names()) == 0)
        return(shiny::tags$div(
          class = "alert alert-info",
          "No stored results yet. Calculate Tm on the Sequence Tm panel first."))
      if (is.null(rv$msg)) return(NULL)
      cls <- if (isTRUE(rv$ok)) "alert alert-success" else "alert alert-danger"
      shiny::tags$div(class = cls, rv$msg)
    })

    output$plot_area <- shiny::renderUI({
      kind <- rv$kind
      if (is.null(kind)) return(NULL)
      if (kind == "plotly") plotly::plotlyOutput(ns("plotly_plot"), height = "600px")
      else shiny::plotOutput(ns("static_plot"), height = "600px")
    })

    output$plotly_plot <- plotly::renderPlotly({
      shiny::req(rv$kind == "plotly")
      rv$obj
    })

    output$static_plot <- shiny::renderPlot({
      shiny::req(rv$kind %in% c("ggplot", "base"))
      if (rv$kind == "ggplot") print(rv$obj) else rv$obj()
    })
  })
}
