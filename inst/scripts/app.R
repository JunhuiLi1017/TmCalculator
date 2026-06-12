library(shiny)
library(shinydashboard)
library(plotly)
library(TmCalculator)
library(GenomicRanges)   # GRanges manipulation
library(IRanges)

## ---- dependency guards -----------------------------------------------------
for (pkg in c("TmCalculator", "GenomicRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' is required but not installed.", pkg))
}
# Soft guards: only needed for some plot types; warn rather than stop.
.has <- function(p) requireNamespace(p, quietly = TRUE)

# ===========================================================================
# UI
# ===========================================================================
ui <- dashboardPage(
  dashboardHeader(
    title = tags$a(href = "https://github.com/JunhuiLi1017/TmCalculator",
                   "TmCalculator",
                   style = "color: black; font-family: Arial, Helvetica, sans-serif;"),
    titleWidth = 220
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Calculator", tabName = "calculator", icon = icon("calculator")),
      menuItem("Compare groups", tabName = "compare", icon = icon("layer-group")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),

  dashboardBody(
    tabItems(

      # ===================== CALCULATOR TAB =================================
      tabItem(
        tabName = "calculator",
        fluidRow(
          # ---- Input -------------------------------------------------------
          box(
            title = "Sequence Input", width = 6,
            selectInput("method", "Calculation Method",
                        choices = c("Nearest Neighbor (tm_nn)" = "tm_nn",
                                    "GC Content (tm_gc)"        = "tm_gc",
                                    "Wallace Rule (tm_wallace)" = "tm_wallace"),
                        selected = "tm_nn"),
            radioButtons("input_type", "Input Type",
                         choices = c("Direct Input"         = "direct",
                                     "FASTA File"           = "fasta",
                                     "Genomic Coordinates"  = "genomic_coordinates",
                                     "BSgenome (genome-wide)" = "bsgenome"),
                         selected = "direct"),

            ## Direct input
            conditionalPanel(
              condition = "input.input_type == 'direct' && input.method == 'tm_nn'",
              textInput("input_seq", "Primer Sequence (5' to 3')",
                        placeholder = "e.g. AACAGACT or AACAGACT, CGTGCATG"),
              textInput("rev_input_seq", "Template Sequence (5' to 3')",
                        placeholder = "e.g. AGTCTGTT or AGTCTGTT, CATGCACG")
            ),
            conditionalPanel(
              condition = "input.input_type == 'direct' && (input.method == 'tm_gc' || input.method == 'tm_wallace')",
              textInput("input_seq_gc", "Sequence (5' to 3')",
                        placeholder = "e.g. AGTCTGTT or AGTCTGTT, CATGCACG")
            ),

            ## FASTA
            conditionalPanel(
              condition = "input.input_type == 'fasta' && input.method == 'tm_nn'",
              fileInput("fasta_file", "Upload FASTA File",
                        accept = c(".fasta", ".fa", ".txt")),
              fileInput("fasta_file_complement", "Upload reverse-complement FASTA",
                        accept = c(".fasta", ".fa", ".txt"))
            ),
            conditionalPanel(
              condition = "input.input_type == 'fasta' && (input.method == 'tm_gc' || input.method == 'tm_wallace')",
              fileInput("fasta_file_gc", "Upload FASTA File",
                        accept = c(".fasta", ".fa", ".txt"))
            ),

            ## Genomic coordinates
            conditionalPanel(
              condition = "input.input_type == 'genomic_coordinates' && input.method == 'tm_nn'",
              textInput("genomic_coordinates", "Genomic Coordinates",
                        placeholder = "chr1:1000-1200:+:BSgenome.Hsapiens.UCSC.hg38:region1"),
              textInput("genomic_coordinates_complement",
                        "Coordinates of reverse complement",
                        placeholder = "chr1:1000-1200:-:BSgenome.Hsapiens.UCSC.hg38")
            ),
            conditionalPanel(
              condition = "input.input_type == 'genomic_coordinates' && (input.method == 'tm_gc' || input.method == 'tm_wallace')",
              textInput("genomic_coordinates_gc", "Genomic Coordinates",
                        placeholder = "chr1:1000-1200:+:BSgenome.Hsapiens.UCSC.hg38:region1")
            ),
            conditionalPanel(
              condition = "input.input_type == 'genomic_coordinates'",
              helpText("Field order: chr:start-end:strand:pkg_name:region_id ",
                       "(region_id optional). The named BSgenome package must be installed.")
            ),

            ## BSgenome genome-wide tiling
            conditionalPanel(
              condition = "input.input_type == 'bsgenome'",
              uiOutput("genome_picker"),
              textInput("bs_chroms", "Chromosome(s)",
                        placeholder = "chr1 (blank = all; comma-separated)"),
              fluidRow(
                column(6, numericInput("bs_window", "Window (bp)", 200, min = 1)),
                column(6, numericInput("bs_slide", "Slide (bp)", 50, min = 1))
              ),
              fluidRow(
                column(6, numericInput("bs_start", "Start (opt.)", NA)),
                column(6, numericInput("bs_end", "End (opt.)", NA))
              ),
              selectInput("bs_strand", "Strand", c("+", "-"), "+"),
              selectInput("bs_trimN", "N handling",
                          c("Trim ends" = "ends", "Filter windows" = "filter",
                            "Keep (none)" = "none"), "ends"),
              numericInput("bs_maxN", "Max N fraction", 0.10, 0, 1, 0.05),
              helpText("Genome-wide tiling can create many windows; restrict the ",
                       "chromosome or start/end range to keep runtimes sensible.")
            ),

            actionButton("calculate", "Calculate Tm", class = "btn-primary")
          ),

          # ---- Method parameters ------------------------------------------
          box(
            title = "Specific Parameters", width = 6,
            checkboxInput("ambiguous", "Ambiguous Bases", FALSE),

            ## GC parameters incl. custom A/B/C/D
            conditionalPanel(
              condition = "input.method == 'tm_gc'",
              checkboxInput("gc_custom", "Custom coefficients (A, B, C, D)", FALSE),
              conditionalPanel(
                condition = "input.gc_custom == false",
                selectInput("variant", "Variant",
                            choices = c("Primer3Plus", "Chester1993", "QuikChange",
                                        "Schildkraut1965", "Wetmur1991_MELTING",
                                        "Wetmur1991_RNA", "Wetmur1991_RNA/DNA",
                                        "vonAhsen2001"),
                            selected = "Primer3Plus")
              ),
              conditionalPanel(
                condition = "input.gc_custom == true",
                helpText("Tm = A + B * %GC - C / N - D * %mismatch"),
                fluidRow(
                  column(6, numericInput("gc_A", "A (intercept)", 81.5)),
                  column(6, numericInput("gc_B", "B (x %GC)", 0.41))
                ),
                fluidRow(
                  column(6, numericInput("gc_C", "C (/ N)", 600)),
                  column(6, numericInput("gc_D", "D (mismatch)", 1))
                )
              ),
              checkboxInput("gc_mismatch", "Count 'X' as mismatch", TRUE)
            ),

            ## NN parameters
            conditionalPanel(
              condition = "input.method == 'tm_nn'",
              numericInput("shift", "Alignment offset (primer vs template)", 0),
              numericInput("dnac_high", "High DNA concentration (nM)", 25, min = 0),
              numericInput("dnac_low", "Low DNA concentration (nM)", 25, min = 0),
              checkboxInput("self_comp", "Self-complementary", FALSE),
              selectInput("nn_table", "Nearest Neighbor Table",
                          choices = c("DNA_NN_SantaLucia_2004", "DNA_NN_Breslauer_1986",
                                      "DNA_NN_Sugimoto_1996", "DNA_NN_Allawi_1998",
                                      "RNA_NN_Freier_1986", "RNA_NN_Xia_1998",
                                      "RNA_NN_Chen_2012", "RNA_DNA_NN_Sugimoto_1995"),
                          selected = "DNA_NN_SantaLucia_2004"),
              selectInput("tmm_table", "Terminal Mismatch Table",
                          choices = c("DNA_TMM_Bommarito_2000"),
                          selected = "DNA_TMM_Bommarito_2000"),
              selectInput("imm_table", "Internal Mismatch Table",
                          choices = c("DNA_IMM_Peyret_1999"),
                          selected = "DNA_IMM_Peyret_1999"),
              selectInput("de_table", "Dangling Ends Table",
                          choices = c("DNA_DE_Bommarito_2000", "RNA_DE_Turner_2010"),
                          selected = "DNA_DE_Bommarito_2000")
            )
          )
        ),

        # ---- Salt + chemical corrections ----------------------------------
        fluidRow(
          conditionalPanel(
            condition = "input.method == 'tm_nn' || input.method == 'tm_gc'",
            box(
              title = "Salt Corrections", width = 6,
              selectInput("salt_method", "Salt Correction Method",
                          choices = c("Schildkraut2010", "Wetmur1991",
                                      "SantaLucia1996", "SantaLucia1998-1",
                                      "Owczarzy2004", "Owczarzy2008"),
                          selected = "Schildkraut2010"),
              numericInput("Na", "Na+ (mM)", 50, min = 0),
              numericInput("K", "K+ (mM)", 0, min = 0),
              numericInput("Tris", "Tris (mM)", 0, min = 0),
              numericInput("Mg", "Mg2+ (mM)", 0, min = 0),
              numericInput("dNTPs", "dNTPs (mM)", 0, min = 0)
            ),
            box(
              title = "Chemical Corrections", width = 6,
              numericInput("DMSO", "DMSO (%)", 0, min = 0, max = 100),
              numericInput("formamide", "Formamide", 0, min = 0),
              selectInput("formamide_unit", "Formamide Unit",
                          choices = c("Percent (%)" = "percent", "Molar (M)" = "molar"),
                          selected = "percent"),
              selectInput("formamide_factor", "Formamide Factor",
                          choices = c("0.65", "0.6", "0.72"), selected = "0.65"),
              selectInput("DMSO_factor", "DMSO Factor",
                          choices = c("0.75", "0.5", "0.6", "0.65", "0.675"),
                          selected = "0.75")
            )
          )
        ),

        # ---- Results ------------------------------------------------------
        fluidRow(
          box(
            title = "Results", width = 12,
            verbatimTextOutput("results"),
            downloadButton("download_results", "Download Results (CSV)")
          ),

          box(
            title = "Visualization", width = 12,
            fluidRow(
              column(3,
                selectInput("plot_type", "Plot Type",
                            choices = c("Linear (Tm vs position)" = "linear",
                                        "Karyotype (interactive)" = "karyotype",
                                        "Heatmap (interactive)"   = "heatmap",
                                        "Genome Tracks"           = "genome_tracks",
                                        "Multi-omics tracks"      = "multiomics"),
                            selected = "linear")
              ),
              column(9,
                ## Linear
                conditionalPanel(
                  condition = "input.plot_type == 'linear'",
                  selectInput("lin_x", "X axis",
                              c("index", "label", "position"), "index"),
                  selectInput("lin_palette", "Palette",
                              c("viridis", "magma", "plasma", "inferno", "cividis"),
                              "viridis"),
                  checkboxInput("lin_line", "Connect points", FALSE),
                  checkboxInput("lin_facet", "Facet by chromosome", FALSE)
                ),

                ## Karyotype
                conditionalPanel(
                  condition = "input.plot_type == 'karyotype'",
                  textInput("karyotype_genome_assembly", "Genome Assembly",
                            placeholder = "e.g. hg38, hg19, mm10"),
                  selectInput("karyotype_chromosomes", "Chromosomes to Plot",
                              choices = NULL, multiple = TRUE),
                  textInput("karyotype_colors", "Chromosome Colors",
                            placeholder = "chr1=#FF0000,chr2=blue"),
                  textInput("karyotype_shapes", "Chromosome Shapes",
                            placeholder = "chr1=16,chr2=17"),
                  selectInput("karyotype_plot_type", "Layout (plot_type)",
                              choices = as.character(1:7), selected = "1"),
                  numericInput("karyotype_point_cex", "Point Size", 1.5, 0.1, 3, 0.1),
                  numericInput("karyotype_tick_dist", "Tick Distance (bp)",
                               10000000, 100000, 100000000, 1000000)
                ),

                ## Heatmap
                conditionalPanel(
                  condition = "input.plot_type == 'heatmap'",
                  textInput("heatmap_genome_assembly", "Genome Assembly",
                            placeholder = "e.g. hg38"),
                  selectInput("heatmap_chromosomes", "Chromosome(s) to Plot",
                              choices = NULL, multiple = TRUE),
                  textInput("heatmap_zoom", "Zoom Region",
                            placeholder = "chr1:1000000-2000000"),
                  selectInput("heatmap_plot_type", "Heatmap Type",
                              c("Karyogram" = "karyogram", "Faceted" = "faceted"),
                              "karyogram"),
                  selectInput("heatmap_color_palette", "Color Palette",
                              c("viridis", "magma", "plasma", "inferno", "cividis"),
                              "viridis"),
                  textInput("heatmap_title", "Plot Title", "Tm Values")
                ),

                ## Genome tracks
                conditionalPanel(
                  condition = "input.plot_type == 'genome_tracks'",
                  textInput("genome_tracks_genome_assembly", "Genome Assembly",
                            placeholder = "e.g. hg38"),
                  selectInput("genome_tracks_chromosomes", "Chromosome to Plot",
                              choices = NULL, multiple = FALSE),
                  textInput("genome_tracks_zoom", "Zoom Region",
                            placeholder = "chr1:1000000-2000000"),
                  selectInput("genome_tracks_track_type", "Track Type",
                              c("gradient", "points", "line", "bars"), "gradient"),
                  selectInput("genome_tracks_color_palette", "Color Palette",
                              c("viridis", "magma", "plasma", "inferno", "cividis"),
                              "viridis"),
                  checkboxInput("genome_tracks_show_ideogram", "Show Ideogram", TRUE)
                ),

                ## Multi-omics
                conditionalPanel(
                  condition = "input.plot_type == 'multiomics'",
                  fileInput("multi_feat_file", "Feature ranges (BED / CSV)",
                            accept = c(".bed", ".csv", ".tsv", ".txt")),
                  selectInput("multi_feat_col", "Feature value column", choices = NULL),
                  selectInput("multi_engine", "Track engine",
                              c("Linear genome" = "linear",
                                "Circular (circos)" = "circos",
                                "Karyotype" = "karyo"), "linear"),
                  conditionalPanel(
                    condition = "input.multi_engine != 'circos'",
                    textInput("multi_assembly", "Genome assembly", placeholder = "hg38"),
                    textInput("multi_chroms", "Chromosomes (comma-separated)")
                  ),
                  conditionalPanel(
                    condition = "input.multi_engine == 'circos'",
                    textInput("multi_genome_name", "Genome name", "genome"),
                    numericInput("multi_genome_size",
                                 "Genome size (bp, single chromosome)", NA)
                  ),
                  selectInput("multi_strategy", "Integration strategy",
                              c("overlap", "nearest", "window", "bin"), "overlap"),
                  actionButton("multi_integrate", "Integrate (annotate result)",
                               class = "btn-outline-secondary")
                ),
                actionButton("draw", "Draw plot", class = "btn-primary")
              )
            ),
            conditionalPanel(
              condition = "input.plot_type == 'heatmap' || input.plot_type == 'karyotype'",
              plotlyOutput("plotly_output", height = "600px")
            ),
            conditionalPanel(
              condition = "input.plot_type != 'heatmap' && input.plot_type != 'karyotype'",
              plotOutput("plot_output", height = "600px")
            ),
            downloadButton("download_plot", "Download Plot (PNG)")
          )
        )
      ),

      # ===================== COMPARE TAB ===================================
      tabItem(
        tabName = "compare",
        fluidRow(
          box(
            title = "Group definition & test", width = 4,
            helpText("Operates on the most recent Tm calculation."),
            selectInput("cmp_target", "Value(s) to compare",
                        choices = NULL, multiple = TRUE),
            radioButtons("cmp_group_mode", "Add grouping metadata",
                         c("Quantile bins of a column" = "quantile",
                           "Threshold split of a column" = "threshold",
                           "Existing metadata column" = "existing",
                           "Paste group labels" = "paste"),
                         "quantile"),
            conditionalPanel(
              condition = "input.cmp_group_mode == 'quantile'",
              selectInput("cmp_qcol", "Numeric column", choices = NULL),
              numericInput("cmp_qn", "Number of bins", 3, 2, 10, 1)
            ),
            conditionalPanel(
              condition = "input.cmp_group_mode == 'threshold'",
              selectInput("cmp_tcol", "Numeric column", choices = NULL),
              numericInput("cmp_tval", "Threshold", NA)
            ),
            conditionalPanel(
              condition = "input.cmp_group_mode == 'existing'",
              selectInput("cmp_ecol", "Grouping column", choices = NULL)
            ),
            conditionalPanel(
              condition = "input.cmp_group_mode == 'paste'",
              textAreaInput("cmp_plabels", "Labels (one per region)", rows = 4,
                            placeholder = "high\nhigh\nlow\nlow")
            ),
            selectInput("cmp_method", "Test",
                        c("Wilcoxon / Kruskal-Wallis" = "wilcoxon",
                          "t-test / ANOVA" = "t.test"), "wilcoxon"),
            selectInput("cmp_alt", "Alternative (2 groups)",
                        c("two.sided", "less", "greater"), "two.sided"),
            checkboxInput("cmp_posthoc", "Post-hoc pairwise (>= 3 groups)", TRUE),
            selectInput("cmp_padj", "p-adjust",
                        c("holm", "BH", "bonferroni", "hochberg", "hommel",
                          "BY", "fdr", "none"), "holm"),
            actionButton("cmp_run", "Compare", class = "btn-primary")
          ),
          box(
            title = "Omnibus", width = 8, tableOutput("cmp_results")
          )
        ),
        fluidRow(
          box(title = "Group summary", width = 6, tableOutput("cmp_summary")),
          box(title = "Pairwise", width = 6, tableOutput("cmp_pairwise"))
        )
      ),

      # ===================== ABOUT TAB =====================================
      tabItem(
        tabName = "about",
        box(
          title = "About TmCalculator", width = 12,
          tags$p("Interactive front end for the ", tags$b("TmCalculator"),
                 " package for nucleic-acid melting-temperature analysis."),
          tags$h4("Calculator"),
          tags$p("Compute Tm from direct sequences, a FASTA upload, genomic ",
                 "coordinate strings, or genome-wide tiling of an installed ",
                 "BSgenome package. Methods: tm_nn, tm_gc (incl. custom ",
                 "A/B/C/D coefficients), and tm_wallace."),
          tags$h4("Visualization"),
          tags$p("Linear, karyotype, heatmap, and genome-track views of the ",
                 "result, plus multi-omics integration (integrate_granges) with ",
                 "multi-track linear / circular / karyotype genome plots."),
          tags$h4("Compare groups"),
          tags$p("Attach grouping metadata to the result and test whether Tm/GC ",
                 "differ across groups via compare_groups."),
          tags$p(tags$small("Author: Junhui Li (ORCID 0000-0003-3973-1700). ",
                            "Co-author: Lihua Julie Zhu."))
        )
      )
    )
  )
)

# ===========================================================================
# SERVER
# ===========================================================================
server <- function(input, output, session) {

  result <- reactiveVal(NULL)     # latest tm_calculate() output (has $gr)

  ## ---- helpers ------------------------------------------------------------
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

  parse_named_vector <- function(s, type = "character") {
    if (is.null(s) || s == "") return(NULL)
    pairs <- strsplit(s, ",")[[1]]
    nm <- character(0); val <- character(0)
    for (p in pairs) {
      parts <- strsplit(p, "=")[[1]]
      if (length(parts) == 2) { nm <- c(nm, trimws(parts[1])); val <- c(val, trimws(parts[2])) }
    }
    if (length(nm) == 0) return(NULL)
    if (type == "numeric") val <- as.numeric(val)
    setNames(val, nm)
  }

  chrom_vec <- function(txt) {
    if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
    v <- trimws(strsplit(txt, ",")[[1]]); v[nzchar(v)]
  }

  read_features <- function(path) {
    first <- readLines(path, n = 1L, warn = FALSE)
    sep <- if (grepl("\t", first)) "\t" else ","
    df <- utils::read.csv(path, sep = sep, header = TRUE,
                          stringsAsFactors = FALSE, check.names = TRUE)
    nms <- tolower(names(df))
    ci <- which(nms %in% c("seqnames", "chrom", "chr", "chromosome"))[1]
    si <- which(nms %in% c("start", "chromstart", "begin"))[1]
    ei <- which(nms %in% c("end", "chromend", "stop"))[1]
    if (is.na(ci) || is.na(si) || is.na(ei)) {
      df <- utils::read.csv(path, sep = sep, header = FALSE, stringsAsFactors = FALSE)
      gr <- GenomicRanges::GRanges(as.character(df[[1]]),
              IRanges::IRanges(as.integer(df[[2]]) + 1L, as.integer(df[[3]])))
      if (ncol(df) > 3) for (j in 4:ncol(df)) {
        v <- df[[j]]; num <- suppressWarnings(as.numeric(v))
        GenomicRanges::mcols(gr)[[paste0("V", j)]] <-
          if (all(is.na(num) == is.na(v))) num else v
      }
      return(gr)
    }
    gr <- GenomicRanges::GRanges(as.character(df[[ci]]),
            IRanges::IRanges(as.integer(df[[si]]), as.integer(df[[ei]])))
    meta <- setdiff(seq_len(ncol(df)), c(ci, si, ei))
    for (j in meta) {
      v <- df[[j]]; num <- suppressWarnings(as.numeric(v))
      GenomicRanges::mcols(gr)[[names(df)[j]]] <-
        if (all(is.na(num) == is.na(v))) num else v
    }
    gr
  }

  ## ---- installed BSgenome picker -----------------------------------------
  output$genome_picker <- renderUI({
    pkgs <- tryCatch(rownames(utils::installed.packages()),
                     error = function(e) character(0))
    genomes <- sort(grep("^BSgenome\\.", pkgs, value = TRUE))
    if (length(genomes) == 0)
      helpText("No BSgenome.* packages installed. e.g. ",
               tags$code("BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"))
    else
      selectInput("bs_genome", "BSgenome package", genomes, genomes[1])
  })

  ## ---- build tm_calculate argument list ----------------------------------
  build_input_seq <- function() {
    it <- input$input_type; m <- input$method
    if (it == "direct") {
      raw <- if (m == "tm_nn") input$input_seq else input$input_seq_gc
      if (is.null(raw) || !nzchar(trimws(raw))) stop("Enter a sequence.")
      seqs <- trimws(strsplit(raw, ",")[[1]]); seqs <- seqs[nzchar(seqs)]
      comp <- NULL
      if (m == "tm_nn" && !is.null(input$rev_input_seq) && nzchar(trimws(input$rev_input_seq)))
        comp <- trimws(strsplit(input$rev_input_seq, ",")[[1]])
      return(list(seq = toupper(seqs), comp = comp))
    }
    if (it == "fasta") {
      f <- if (m == "tm_nn") input$fasta_file else input$fasta_file_gc
      req(f)
      comp <- NULL
      if (m == "tm_nn" && !is.null(input$fasta_file_complement))
        comp <- input$fasta_file_complement$datapath
      return(list(seq = f$datapath, comp = comp))
    }
    if (it == "genomic_coordinates") {
      raw <- if (m == "tm_nn") input$genomic_coordinates else input$genomic_coordinates_gc
      if (is.null(raw) || !nzchar(trimws(raw))) stop("Enter coordinate string(s).")
      coords <- trimws(strsplit(raw, "\n")[[1]]); coords <- coords[nzchar(coords)]
      comp <- NULL
      if (m == "tm_nn" && !is.null(input$genomic_coordinates_complement) &&
          nzchar(trimws(input$genomic_coordinates_complement)))
        comp <- trimws(strsplit(input$genomic_coordinates_complement, "\n")[[1]])
      return(list(seq = coords, comp = comp))
    }
    if (it == "bsgenome") {
      req(input$bs_genome)
      if (!requireNamespace(input$bs_genome, quietly = TRUE))
        stop(sprintf("Package '%s' is not installed.", input$bs_genome))
      bsg <- getExportedValue(input$bs_genome, input$bs_genome)
      st <- if (is.na(input$bs_start)) NULL else input$bs_start
      en <- if (is.na(input$bs_end)) NULL else input$bs_end
      coords <- TmCalculator::make_genomiccoord(
        bsgenome = bsg, chromosomes = chrom_vec(input$bs_chroms),
        window = as.integer(input$bs_window), slide = as.integer(input$bs_slide),
        start = st, end = en, strand = input$bs_strand,
        trim_N = input$bs_trimN, max_N_frac = input$bs_maxN,
        genome_pkg_name = input$bs_genome, as_vector = TRUE, verbose = FALSE)
      if (length(coords) == 0) stop("No windows produced; relax filters/range.")
      return(list(seq = coords, comp = NULL))
    }
    stop("Unknown input type.")
  }

  build_tm_args <- function(spec) {
    a <- list(input_seq = spec$seq, method = input$method,
              ambiguous = isTRUE(input$ambiguous))
    if (!is.null(spec$comp)) a$complement_seq <- spec$comp

    if (input$method == "tm_nn") {
      a$shift <- input$shift
      a$nn_table <- input$nn_table
      a$tmm_table <- input$tmm_table
      a$imm_table <- input$imm_table
      a$de_table <- input$de_table
      a$dnac_high <- input$dnac_high
      a$dnac_low <- input$dnac_low
      a$self_comp <- isTRUE(input$self_comp)
    }
    if (input$method == "tm_gc") {
      a$mismatch <- isTRUE(input$gc_mismatch)
      if (isTRUE(input$gc_custom))
        a$userset <- c(input$gc_A, input$gc_B, input$gc_C, input$gc_D)
      else
        a$variant <- input$variant
    }
    if (input$method != "tm_wallace") {
      a$Na <- input$Na; a$K <- input$K; a$Tris <- input$Tris
      a$Mg <- input$Mg; a$dNTPs <- input$dNTPs
      a$salt_method <- input$salt_method
      a$DMSO <- input$DMSO
      a$formamide_unit <- list(value = input$formamide, unit = input$formamide_unit)
      a$dmso_factor <- as.numeric(input$DMSO_factor)
      a$formamide_factor <- as.numeric(input$formamide_factor)
    }
    a
  }

  ## ---- run calculation ----------------------------------------------------
  observeEvent(input$calculate, {
    progress <- shiny::Progress$new(); on.exit(progress$close())
    progress$set(message = "Calculating Tm", value = 0.2)
    tryCatch({
      spec <- build_input_seq()
      out <- do.call(TmCalculator::tm_calculate, build_tm_args(spec))
      progress$set(message = "Done", value = 1)
      result(out)
    }, error = function(e) showNotification(paste("Error:", e$message),
                                            type = "error", duration = 8))
  })

  output$results <- renderPrint({
    req(result()); print(result())
  })

  output$download_results <- downloadHandler(
    filename = function() paste0("tm_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(result())
      df <- result()$df %||% as.data.frame(result()$gr)
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  ## ---- sync chromosome / column pickers ----------------------------------
  observe({
    req(result())
    gr <- result()$gr; req(inherits(gr, "GRanges"))
    chrs <- unique(as.character(GenomicRanges::seqnames(gr)))
    updateSelectInput(session, "karyotype_chromosomes", choices = chrs, selected = chrs)
    updateSelectInput(session, "heatmap_chromosomes", choices = chrs, selected = chrs)
    updateSelectInput(session, "genome_tracks_chromosomes", choices = chrs, selected = chrs[1])

    cols <- names(GenomicRanges::mcols(gr))
    numc <- cols[vapply(cols, function(c) is.numeric(GenomicRanges::mcols(gr)[[c]]), logical(1))]
    tgt <- intersect(c("Tm", "GC"), numc)
    updateSelectInput(session, "cmp_target", choices = numc,
                      selected = if (length(tgt)) tgt else numc[1])
    updateSelectInput(session, "cmp_qcol", choices = numc, selected = numc[1])
    updateSelectInput(session, "cmp_tcol", choices = numc, selected = numc[1])
    updateSelectInput(session, "cmp_ecol", choices = cols, selected = cols[1])
  })

  ## ---- multi-omics feature upload ----------------------------------------
  features <- reactive({ req(input$multi_feat_file); read_features(input$multi_feat_file$datapath) })
  observeEvent(features(), {
    gr <- features(); cols <- names(GenomicRanges::mcols(gr))
    numc <- cols[vapply(cols, function(c) is.numeric(GenomicRanges::mcols(gr)[[c]]), logical(1))]
    updateSelectInput(session, "multi_feat_col",
                      choices = if (length(numc)) numc else cols,
                      selected = if (length(numc)) numc[1] else cols[1])
  })

  observeEvent(input$multi_integrate, {
    tryCatch({
      gr <- result()$gr; req(gr)
      ann <- TmCalculator::integrate_granges(gr, features(), strategy = input$multi_strategy)
      cur <- result(); cur$gr <- ann; cur$df <- as.data.frame(ann); result(cur)
      showNotification(sprintf("Integrated %d feature column(s) into the result.",
                               length(names(GenomicRanges::mcols(ann))) -
                                 length(names(GenomicRanges::mcols(gr)))),
                       type = "message")
    }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
  })

  build_tracks <- function(gr_tm, gr_feat, feat_col) {
    tl <- list(list(type = "points", data = gr_tm, value_col = "Tm",
                    col = "#e41a1c", name = "Tm"))
    if ("GC" %in% names(GenomicRanges::mcols(gr_tm)))
      tl[[length(tl) + 1L]] <- list(type = "line", data = gr_tm, value_col = "GC",
                                    col = "#377eb8", name = "GC")
    if (!is.null(gr_feat) && !is.null(feat_col) &&
        feat_col %in% names(GenomicRanges::mcols(gr_feat)))
      tl[[length(tl) + 1L]] <- list(type = "line", data = gr_feat,
                                    value_col = feat_col, col = "#4daf4a", name = feat_col)
    tl
  }

  ## ---- drawing ------------------------------------------------------------
  draw_obj <- reactiveVal(NULL)   # holds a function() for base plots / ggplot

  observeEvent(input$draw, {
    tryCatch({
      gr <- result()$gr; req(gr)
      pt <- input$plot_type

      if (pt == "linear") {
        p <- TmCalculator::plot_tm_linear(gr, x_axis = input$lin_x,
               color_palette = input$lin_palette, add_line = isTRUE(input$lin_line),
               facet_by_chr = isTRUE(input$lin_facet))
        draw_obj(function() print(p))
      } else if (pt == "genome_tracks") {
        p <- TmCalculator::plot_tm_genome_tracks(gr,
               chromosome_to_plot = input$genome_tracks_chromosomes,
               genome_assembly = if (nzchar(input$genome_tracks_genome_assembly))
                 input$genome_tracks_genome_assembly else NULL,
               color_palette = input$genome_tracks_color_palette,
               track_type = input$genome_tracks_track_type,
               show_ideogram = isTRUE(input$genome_tracks_show_ideogram),
               zoom = if (nzchar(input$genome_tracks_zoom)) input$genome_tracks_zoom else NULL)
        draw_obj(function() print(p))
      } else if (pt == "multiomics") {
        tl <- build_tracks(gr, features(), input$multi_feat_col)
        if (input$multi_engine == "linear") {
          draw_obj(function() TmCalculator::plot_genome_track(
            genome_name = if (nzchar(input$multi_assembly)) input$multi_assembly else "genome",
            track_list = tl, chromosomes = chrom_vec(input$multi_chroms)))
        } else if (input$multi_engine == "circos") {
          gs <- if (is.na(input$multi_genome_size)) max(GenomicRanges::end(gr)) else input$multi_genome_size
          draw_obj(function() TmCalculator::plot_circos_genome(
            genome_name = input$multi_genome_name %||% "genome",
            genome_size = gs, track_list = tl))
        } else {
          draw_obj(function() TmCalculator::plot_karyotype_genome(
            genome_assembly = if (nzchar(input$multi_assembly)) input$multi_assembly else "hg38",
            track_list = tl, chromosomes = chrom_vec(input$multi_chroms) %||% "auto"))
        }
      }
      # heatmap / karyotype are interactive -> handled in renderPlotly
    }, error = function(e) showNotification(paste("Error:", e$message),
                                            type = "error", duration = 8))
  })

  output$plot_output <- renderPlot({
    req(input$plot_type %in% c("linear", "genome_tracks", "multiomics"))
    f <- draw_obj(); req(f); f()
  })

  output$plotly_output <- renderPlotly({
    req(result()); gr <- result()$gr; req(gr)
    if (input$plot_type == "heatmap") {
      TmCalculator::plot_tm_heatmap_interactive(gr,
        genome_assembly = if (nzchar(input$heatmap_genome_assembly)) input$heatmap_genome_assembly else NULL,
        chromosome_to_plot = input$heatmap_chromosomes,
        plot_type = input$heatmap_plot_type,
        color_palette = input$heatmap_color_palette,
        title_name = input$heatmap_title,
        zoom = if (nzchar(input$heatmap_zoom)) input$heatmap_zoom else NULL)
    } else if (input$plot_type == "karyotype") {
      TmCalculator::plot_tm_karyotype_interactive(gr,
        chromosomes = input$karyotype_chromosomes,
        genome_assembly = if (nzchar(input$karyotype_genome_assembly)) input$karyotype_genome_assembly else NULL,
        colors = parse_named_vector(input$karyotype_colors, "character"),
        shapes = parse_named_vector(input$karyotype_shapes, "numeric"),
        plot_type = as.numeric(input$karyotype_plot_type),
        point_cex = input$karyotype_point_cex,
        tick_dist = input$karyotype_tick_dist)
    }
  })

  output$download_plot <- downloadHandler(
    filename = function() paste0("tm_plot_", Sys.Date(), ".png"),
    content = function(file) {
      grDevices::png(file, width = 9, height = 7, units = "in", res = 300)
      f <- draw_obj(); if (!is.null(f)) f()
      grDevices::dev.off()
    }
  )

  ## ===================== COMPARE ==========================================
  build_group <- function(gr) {
    n <- length(gr); mode <- input$cmp_group_mode
    if (mode == "quantile") {
      req(input$cmp_qcol)
      x <- GenomicRanges::mcols(gr)[[input$cmp_qcol]]
      nb <- max(2L, as.integer(input$cmp_qn))
      brk <- unique(stats::quantile(x, probs = seq(0, 1, length.out = nb + 1L), na.rm = TRUE))
      if (length(brk) < 3L) stop("Too few distinct values for that many bins.")
      return(as.character(cut(x, brk, include.lowest = TRUE,
                              labels = paste0("Q", seq_len(length(brk) - 1L)))))
    }
    if (mode == "threshold") {
      req(input$cmp_tcol)
      if (is.na(input$cmp_tval)) stop("Enter a numeric threshold.")
      x <- GenomicRanges::mcols(gr)[[input$cmp_tcol]]
      return(ifelse(x >= input$cmp_tval, paste0(">=", input$cmp_tval), paste0("<", input$cmp_tval)))
    }
    if (mode == "existing") {
      req(input$cmp_ecol); return(as.character(GenomicRanges::mcols(gr)[[input$cmp_ecol]]))
    }
    labs <- trimws(strsplit(input$cmp_plabels %||% "", "\n")[[1]]); labs <- labs[nzchar(labs)]
    if (length(labs) != n) stop(sprintf("Pasted %d label(s) but result has %d region(s).", length(labs), n))
    labs
  }

  cmp_out <- reactiveVal(NULL)
  observeEvent(input$cmp_run, {
    tryCatch({
      gr <- result()$gr
      if (is.null(gr)) stop("Run a Tm calculation first.")
      if (length(input$cmp_target) == 0) stop("Select at least one value to compare.")
      GenomicRanges::mcols(gr)[["group"]] <- build_group(gr)
      cmp_out(TmCalculator::compare_groups(gr, target = input$cmp_target,
              method = input$cmp_method, group = "group", alternative = input$cmp_alt,
              posthoc = isTRUE(input$cmp_posthoc), p.adjust.method = input$cmp_padj))
    }, error = function(e) { cmp_out(NULL); showNotification(paste("Error:", e$message),
                                                             type = "error", duration = 8) })
  })

  output$cmp_results <- renderTable({ req(cmp_out()); cmp_out()$results })
  output$cmp_summary <- renderTable({ req(cmp_out()); cmp_out()$summary })
  output$cmp_pairwise <- renderTable({
    validate(need(!is.null(cmp_out()) && !is.null(cmp_out()$pairwise),
                  "No pairwise table (needs >= 3 groups with post-hoc on)."))
    cmp_out()$pairwise
  })
}

shinyApp(ui = ui, server = server)
