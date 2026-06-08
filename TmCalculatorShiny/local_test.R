# Local smoke test for TmCalculatorShiny -------------------------------------
# Run this on your own machine (the build sandbox cannot install R/Bioconductor).
#
#   setwd("path/to/TmCalculator/TmCalculatorShiny")
#   source("local_test.R")
#
# It loads the package with devtools, runs the test suite, and launches the app.

if (!requireNamespace("devtools", quietly = TRUE))
  stop("install.packages('devtools') first")

devtools::document()        # regenerate man/ + NAMESPACE from roxygen
devtools::load_all(".")     # attach the package from source

# 1. Unit tests -------------------------------------------------------------
testthat::test_local(".")

# 2. Quick programmatic checks ----------------------------------------------
res <- TmCalculator::tm_calculate(
  c("ATCGTGCGTAGCAGTACGATCAGTAG", "GGGCCCATATATGCGC"),
  method = "tm_nn")
print(as.data.frame(res$gr))

gc_custom <- TmCalculator::tm_calculate(
  "ATCGTGCGTAGCAGTACGATCAGTAG", method = "tm_gc",
  userset = c(81.5, 0.41, 600, 1))     # A, B, C, D
print(as.data.frame(gc_custom$gr))

# 3. Launch the app ---------------------------------------------------------
message("Launching app... close the browser tab / press Esc to stop.")
runTmCalculatorShiny()
