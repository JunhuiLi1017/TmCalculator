#' Compare-groups module UI
#'
#' Consumes a stored Sequence-Tm result, lets the user attach grouping
#' metadata (quantile bins of a numeric column, a threshold split, or pasted
#' labels), then runs \code{TmCalculator::compare_groups}.
#' @noRd
mod_compare_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      title = "Group definition & test",

      shiny::selectInput(ns("result"), "Stored result", choices = NULL),

      shiny::selectInput(ns("target"), "Value(s) to compare",
                         choices = NULL, multiple = TRUE),

      shiny::hr(),
      shiny::tags$b("Add grouping metadata"),
      shiny::radioButtons(
        ns("group_mode"), NULL,
        choices = c("Quantile bins of a column" = "quantile",
                    "Threshold split of a column" = "threshold",
                    "Use existing metadata column" = "existing",
                    "Paste group labels (one per row)" = "paste"),
        selected = "quantile"
      ),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'quantile'", ns("group_mode")),
        shiny::selectInput(ns("q_col"), "Numeric column", choices = NULL),
        shiny::numericInput(ns("q_n"), "Number of bins", value = 3,
                            min = 2, max = 10, step = 1)
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'threshold'", ns("group_mode")),
        shiny::selectInput(ns("t_col"), "Numeric column", choices = NULL),
        shiny::numericInput(ns("t_value"), "Threshold", value = NA)
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'existing'", ns("group_mode")),
        shiny::selectInput(ns("e_col"), "Grouping column", choices = NULL)
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'paste'", ns("group_mode")),
        shiny::textAreaInput(ns("p_labels"), "Labels", rows = 4,
                             placeholder = "high\nhigh\nlow\nlow"),
        shiny::helpText("Provide exactly one label per region, in table order.")
      ),

      shiny::hr(),
      shiny::selectInput(ns("method"), "Test",
                         choices = c("Wilcoxon / Kruskal-Wallis" = "wilcoxon",
                                     "t-test / ANOVA" = "t.test"),
                         selected = "wilcoxon"),
      shiny::selectInput(ns("alternative"), "Alternative (2 groups)",
                         choices = c("two.sided", "less", "greater"),
                         selected = "two.sided"),
      shiny::checkboxInput(ns("posthoc"), "Post-hoc pairwise (>= 3 groups)", TRUE),
      shiny::selectInput(ns("padj"), "p-adjust",
                         choices = c("holm", "BH", "bonferroni", "hochberg",
                                     "hommel", "BY", "fdr", "none"),
                         selected = "holm"),
      shiny::actionButton(ns("run"), "Compare", class = "btn-primary",
                          width = "100%")
    ),

    bslib::card(
      bslib::card_header("Comparison"),
      shiny::uiOutput(ns("status")),
      bslib::navset_card_tab(
        bslib::nav_panel("Omnibus", DT::DTOutput(ns("results"))),
        bslib::nav_panel("Group summary", DT::DTOutput(ns("summary"))),
        bslib::nav_panel("Pairwise", DT::DTOutput(ns("pairwise")))
      )
    )
  )
}

#' Compare-groups module server
#' @param id module id.
#' @param store shared reactiveValues with `$results`.
#' @noRd
mod_compare_server <- function(id, store) {
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

    # Keep column pickers in sync with the chosen result
    shiny::observeEvent(current_gr(), {
      gr <- current_gr()
      shiny::req(gr)
      cols <- names(GenomicRanges::mcols(gr))
      num_cols <- cols[vapply(cols, function(c)
        is.numeric(GenomicRanges::mcols(gr)[[c]]), logical(1))]
      tgt_default <- intersect(c("Tm", "GC"), num_cols)
      shiny::updateSelectInput(session, "target", choices = num_cols,
                               selected = if (length(tgt_default)) tgt_default
                               else num_cols[1])
      shiny::updateSelectInput(session, "q_col", choices = num_cols,
                               selected = num_cols[1])
      shiny::updateSelectInput(session, "t_col", choices = num_cols,
                               selected = num_cols[1])
      shiny::updateSelectInput(session, "e_col", choices = cols,
                               selected = cols[1])
    })

    rv <- shiny::reactiveValues(results = NULL, summary = NULL,
                                pairwise = NULL, msg = NULL, ok = FALSE)

    # Build a group factor of length(gr) according to the chosen mode
    build_group <- function(gr) {
      n <- length(gr)
      mode <- input$group_mode

      if (mode == "quantile") {
        shiny::req(input$q_col)
        x <- GenomicRanges::mcols(gr)[[input$q_col]]
        nb <- max(2L, as.integer(input$q_n))
        brk <- stats::quantile(x, probs = seq(0, 1, length.out = nb + 1L),
                               na.rm = TRUE)
        brk <- unique(brk)
        if (length(brk) < 3L)
          stop("Column has too few distinct values for that many bins.")
        g <- cut(x, breaks = brk, include.lowest = TRUE,
                 labels = paste0("Q", seq_len(length(brk) - 1L)))
        return(as.character(g))
      }

      if (mode == "threshold") {
        shiny::req(input$t_col)
        if (is.null(input$t_value) || is.na(input$t_value))
          stop("Please enter a numeric threshold.")
        x <- GenomicRanges::mcols(gr)[[input$t_col]]
        return(ifelse(x >= input$t_value,
                      paste0(">=", input$t_value),
                      paste0("<", input$t_value)))
      }

      if (mode == "existing") {
        shiny::req(input$e_col)
        return(as.character(GenomicRanges::mcols(gr)[[input$e_col]]))
      }

      if (mode == "paste") {
        labs <- trimws(strsplit(input$p_labels %||% "", "\r?\n")[[1]])
        labs <- labs[nzchar(labs)]
        if (length(labs) != n)
          stop(sprintf("Pasted %d label(s) but the result has %d region(s).",
                       length(labs), n))
        return(labs)
      }

      stop("Unknown grouping mode.")
    }

    shiny::observeEvent(input$run, {
      rv$ok <- FALSE
      res <- tryCatch({
        gr <- current_gr()
        if (is.null(gr)) stop("Select a stored result first.")
        if (length(input$target) == 0) stop("Select at least one value to compare.")
        g <- build_group(gr)
        GenomicRanges::mcols(gr)[["group"]] <- g
        cg <- TmCalculator::compare_groups(
          gr,
          target = input$target,
          method = input$method,
          group = "group",
          alternative = input$alternative,
          posthoc = isTRUE(input$posthoc),
          p.adjust.method = input$padj
        )
        cg
      }, error = function(e) {
        rv$msg <- conditionMessage(e); rv$ok <- FALSE; NULL
      })

      if (is.null(res)) {
        rv$results <- NULL; rv$summary <- NULL; rv$pairwise <- NULL
        return(invisible(NULL))
      }
      rv$results <- res$results
      rv$summary <- res$summary
      rv$pairwise <- res$pairwise
      rv$msg <- "Comparison complete."
      rv$ok <- TRUE
    })

    output$status <- shiny::renderUI({
      if (is.null(rv$msg)) {
        if (length(result_names()) == 0)
          return(shiny::tags$div(
            class = "alert alert-info",
            "No stored results yet. Calculate Tm on the Sequence Tm panel first."
          ))
        return(NULL)
      }
      cls <- if (isTRUE(rv$ok)) "alert alert-success" else "alert alert-danger"
      shiny::tags$div(class = cls, rv$msg)
    })

    output$results <- DT::renderDT({
      shiny::req(rv$results)
      DT::datatable(rv$results, options = list(scrollX = TRUE), rownames = FALSE)
    })
    output$summary <- DT::renderDT({
      shiny::req(rv$summary)
      DT::datatable(rv$summary, options = list(scrollX = TRUE), rownames = FALSE)
    })
    output$pairwise <- DT::renderDT({
      shiny::validate(shiny::need(!is.null(rv$pairwise),
                                  "No pairwise table (need >= 3 groups with post-hoc on)."))
      DT::datatable(rv$pairwise, options = list(scrollX = TRUE), rownames = FALSE)
    })
  })
}
