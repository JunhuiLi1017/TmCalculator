## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo      = TRUE,
  message   = FALSE,
  warning   = FALSE,
  fig.align = "center",
  # html_vignette inherits fig.retina = 2 from rmarkdown, which renders every
  # figure at twice the nominal resolution and then scales it down in the
  # browser. The pixel count, and therefore the size of the base64 blob
  # embedded in the self-contained HTML, is four times larger for no visible
  # gain in a vignette. Setting it to 1 is what keeps the installed size of
  # doc/ within the limit R CMD check reports on.
  fig.retina = 1,
  dpi        = 72,
  fig.width = 7,
  fig.height = 6
)

## ----packages-----------------------------------------------------------------
library(TmCalculator)
library(BSgenome)
library(GenomicRanges)

## ----forge, eval=FALSE--------------------------------------------------------
# library(BSgenomeForge)
# 
# forgeBSgenomeDataPkgFromNCBI(
#   assembly_accession = "GCF_000005845.2",
#   pkg_maintainer     = "Junhui Li <ljh.biostat@gmail.com>",
#   destdir            = "."
# )
# 
# install.packages(
#   "./BSgenome.Ecoli.NCBI.ASM584v2",
#   repos = NULL,
#   type  = "source"
# )

## ----setup-ecoli-bsgenome, message=FALSE, warning=FALSE-----------------------
ecoli_pkg   <- "BSgenome.Ecoli.NCBI.ASM584v2"
genome_obj  <- "Ecoli"   # BSgenomeObjname in DESCRIPTION; not the package name

.ecoli_genome_ready <- function() {
  if (!requireNamespace(ecoli_pkg, quietly = TRUE)) return(FALSE)
  exists(genome_obj, envir = asNamespace(ecoli_pkg), inherits = FALSE)
}

if (!requireNamespace(ecoli_pkg, quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github(
    "JunhuiLi1017/BSgenome.Ecoli.NCBI.ASM584v2",
    upgrade = "never",
    quiet = TRUE
  )
}

if (!.ecoli_genome_ready()) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("BSgenomeForge", quietly = TRUE)) {
    BiocManager::install("BSgenomeForge", ask = FALSE, update = FALSE)
  }
  pkgdir <- BSgenomeForge::forgeBSgenomeDataPkgFromNCBI(
    assembly_accession = "GCF_000005845.2",
    pkg_maintainer     = "Junhui Li <ljh.biostat@gmail.com>",
    destdir            = tempdir()
  )
  utils::install.packages(pkgdir, repos = NULL, type = "source", quiet = TRUE)
  if (ecoli_pkg %in% loadedNamespaces()) {
    unloadNamespace(ecoli_pkg)
  }
}

if (!.ecoli_genome_ready()) {
  stop(
    "Could not load genome object '", genome_obj, "' from package '", ecoli_pkg, "'.\n",
    "Run the manual 'forge' chunk above, or forgeBSgenomeDataPkgFromNCBI() locally.",
    call. = FALSE
  )
}

## ----load-genome--------------------------------------------------------------
ecoli_pkg  <- "BSgenome.Ecoli.NCBI.ASM584v2"
genome_obj <- "Ecoli"

suppressPackageStartupMessages(library(ecoli_pkg, character.only = TRUE))
genome      <- base::get(genome_obj, envir = asNamespace(ecoli_pkg))
genome_name <- ecoli_pkg
chr_name    <- "U00096.3"
chr_length  <- length(genome[[chr_name]])

cat("Chromosome:", chr_name, "\n")
cat("Length:    ", format(chr_length, big.mark = ","), "bp\n")

## ----make-coord---------------------------------------------------------------
runtime0 <- system.time({
  bins_gc <- make_genomiccoord(
    bsgenome    = genome_name,
    chromosomes = chr_name,
    window      = 200L,
    slide       = 200L,
    start       = 1,
    end         = chr_length,
    strand      = "+"
  )
})

cat("Total windows:", length(bins_gc), "\n")
cat(sprintf("Window generation: %.2f s (elapsed)\n", runtime0["elapsed"]))

## ----coor-to-gr---------------------------------------------------------------
input_new <- list(pkg_name = genome_name, seq = bins_gc)
runtime1 <- system.time({
  gr_batch <- to_genomic_ranges_fast(input_new)
})

cat(sprintf(
  "Coordinate resolution: %.2f s (elapsed)\n",
  runtime1["elapsed"]
))

## ----tm-calculate-------------------------------------------------------------
runtime2 <- system.time({
  tm_ASM584v2 <- tm_calculate(
    gr_batch,
    method   = "tm_nn",
    nn_table = "DNA_NN_Breslauer_1986",
    Na       = 50            # mM; standard PCR-like conditions
  )
})

cat(sprintf(
  "Tm calculation: %.2f s (elapsed) for %s windows\n",
  runtime2["elapsed"],
  format(length(bins_gc), big.mark = ",")
))

Tm <- as.data.frame(tm_ASM584v2$gr[, c("Tm", "GC")])
summary(Tm[, c("Tm", "GC")])

## ----track-list---------------------------------------------------------------
# Reference labels: replication origin (ori) and terminus (dif)
label <- data.frame(
  seqnames = genome_name,
  start    = c(3925804, 1590777),
  end      = c(3925804, 1590777),
  label    = c("ori", "dif")
)
data(ecoli_rep_hotspots)
tracks <- list(
  # Ideogram: MutL-AR peaks shown inside the chromosome bar
  list(type = "rect", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
       col = "#2C3E50", bg.col = "grey", name = "MutL-AR",
       legend_font_col = "#2C3E50", ideogram = TRUE, height = 0.5),

  # Sequence thermodynamics
  list(type = "line", data = Tm, value_col = "GC",
       name = "GC content", col = "#4A90E2",
       legend_font_col = "#4A90E2"),
  list(type = "line", data = Tm, value_col = "Tm",
       name = "Melting temp", col = "#E06666",
       legend_font_col = "#E06666", height = 2),

  # Repeat / structural features
  list(type = "line", data = ecoli_rep_hotspots$bins_rep,
       value_col = "count", name = "Microsatellites", col = "#2ECC71",
       legend_font_col = "#2ECC71"),
  list(type = "line", data = ecoli_rep_hotspots$bins_cru,
       value_col = "count", name = "Cruciform", col = "#3B3E6B",
       legend_font_col = "#3B3E6B"),

  # ssDNA regions
  list(data = ecoli_rep_hotspots$ssdna, name = "ssDNA",
       col = "#8E44AD", legend_font_col = "#8E44AD"),

  # GATC methylation sites
  list(type = "line", data = ecoli_rep_hotspots$bins_gatc,
       value_col = "count", name = "GATC sites", col = "#D35400",
       legend_font_col = "#D35400"),

  # Global highlight: translucent bands at MutL-AR peaks across all tracks
  list(type = "highlight", data = ecoli_rep_hotspots$all_peaks_IP_mutH,
       col = "#F1C40F", alpha = 0.18)
)

## ----circular-full, fig.width=8, fig.height=8, fig.cap="Circular genome map of E. coli K-12 MG1655. Concentric rings from outside in: MutL-AR peaks (grey ideogram), GC content, melting temperature, microsatellite density, cruciform sequences, ssDNA regions, and GATC site density. Yellow highlight bands mark MutL-AR peak regions."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks,
  circular    = TRUE,
  label       = label
)

## ----linear-full, fig.width=10, fig.height=5, fig.cap="Linear genome view of E. coli K-12 MG1655. MutL-AR peaks are drawn inside the chromosome ideogram bar."----
plot_genome_track(
  genome_name = genome_name,
  genome_size = chr_length,
  track_list  = tracks
)

