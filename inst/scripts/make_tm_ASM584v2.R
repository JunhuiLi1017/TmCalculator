#!/usr/bin/env Rscript
# make_tm_ASM584v2.R
#
# Generate the tm_ASM584v2 dataset (genome-wide Tm for E. coli K-12 MG1655)
# and add it to ecoli_rep_hotspots as ecoli_rep_hotspots$tm_ASM584v2.
# Requires BSgenome.Ecoli.NCBI.ASM584v2 to be installed.
#
# Run from the package root:
#   Rscript inst/scripts/make_tm_ASM584v2.R

library(TmCalculator)
library(BSgenome)

ecoli_pkg  <- "BSgenome.Ecoli.NCBI.ASM584v2"
genome_obj <- "Ecoli"

suppressPackageStartupMessages(library(ecoli_pkg, character.only = TRUE))
genome      <- base::get(genome_obj, envir = asNamespace(ecoli_pkg))
chr_name    <- "U00096.3"
chr_length  <- length(genome[[chr_name]])

cat("Chromosome:", chr_name, "\n")
cat("Length:    ", format(chr_length, big.mark = ","), "bp\n")

# Step 1 — genomic windows (200 bp, non-overlapping)
bins_gc <- make_genomiccoord(
  bsgenome    = genome,
  chromosomes = chr_name,
  window      = 200L,
  slide       = 200L,
  start       = 1,
  end         = chr_length,
  strand      = "+"
)
cat("Total windows:", length(bins_gc), "\n")

# Step 2 — resolve coordinates
input_new <- list(pkg_name = ecoli_pkg, seq = bins_gc)
gr_batch  <- to_genomic_ranges_fast(input_new)

# Step 3 — compute Tm
tm_ASM584v2 <- tm_calculate(
  gr_batch,
  method   = "tm_nn",
  nn_table = "DNA_NN_SantaLucia_2004",
  Na       = 50
)

cat("Tm range: ", paste(range(tm_ASM584v2$gr$Tm), collapse = " - "), "°C\n")
cat("GC range: ", paste(range(tm_ASM584v2$gr$GC), collapse = " - "), "\n")

# Step 4 — merge into ecoli_rep_hotspots
load("data/ecoli_rep_hotspots.rda")
ecoli_rep_hotspots$tm_ASM584v2 <- tm_ASM584v2
save(ecoli_rep_hotspots, file = "data/ecoli_rep_hotspots.rda", compress = "xz")
cat("Saved ecoli_rep_hotspots$tm_ASM584v2 into data/ecoli_rep_hotspots.rda\n")

# Remove standalone file if it exists
if (file.exists("data/tm_ASM584v2.rda")) {
  file.remove("data/tm_ASM584v2.rda")
  cat("Removed standalone data/tm_ASM584v2.rda\n")
}
