#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# bench_worker.R -- run ONE tool on ONE input size, in a clean process.
#
# Invoked by benchmark_tools.R; not meant to be run by hand. Writing the Tm
# values to disk lets the driver check cross-tool agreement afterwards.
#
#   Rscript bench_worker.R <tool> <seqfile> <n> <outcsv>
#     tool    "tmcalculator" | "rmelting" | "rmelting_nn"
#     seqfile one sequence per line, plain text
#     n       number of sequences to process (first n lines)
#     outcsv  where to write the Tm values
#
# Shared conditions across all tools (see benchmark_tools.R for why):
#   SantaLucia 2004 NN parameters, 50 mM Na+, 25 nM of each strand,
#   Schildkraut/Lifson salt correction, no other cosolutes.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("usage: bench_worker.R <tool> <seqfile> <n> <outcsv>")
tool <- args[1]; seqfile <- args[2]; n <- as.integer(args[3]); outcsv <- args[4]

seqs <- readLines(seqfile, n = n)

## Load the tool BEFORE timing starts: package/JVM startup is not part of
## the per-dataset cost being compared, and Biopython's import is likewise
## excluded on the Python side.
if (tool %in% c("tmcalculator", "tmcalculator_serial")) {
  suppressPackageStartupMessages(library(TmCalculator))
} else if (tool %in% c("rmelting", "rmelting_nn")) {
  suppressPackageStartupMessages(library(rmelting))
} else {
  stop("unknown tool: ", tool)
}

t_prep <- 0

if (tool %in% c("tmcalculator", "tmcalculator_serial")) {
  ## Input preparation (building the GRanges container) is reported
  ## separately: the other two tools take bare strings and have no
  ## equivalent step, so folding it into the compute time would not be a
  ## like-for-like comparison.
  t_prep <- system.time(gr <- to_genomic_ranges(seqs))[["elapsed"]]
  t0 <- proc.time()
  res <- suppressWarnings(tm_calculate(
    gr, method = "tm_nn",
    nn_table = "DNA_NN_SantaLucia_2004",
    salt_method = "Schildkraut2010",
    Na = 50, dnac_high = 25, dnac_low = 25
  ))
  tm <- as.numeric(res$gr$Tm)
  elapsed <- (proc.time() - t0)[["elapsed"]]

} else {
  ## rmelting wraps the MELTING 5 Java program (melting5jars) through rJava
  ## and exposes no vectorised entry point, so one JVM call per sequence is
  ## the only interface available.
  ##
  ## Two configurations are benchmarked, because MELTING's behaviour on
  ## 200 bp windows depends on how it is called:
  ##   "rmelting"    - no method.nn given, so MELTING selects a model by
  ##                   sequence length. Above `size.threshold` (default 60)
  ##                   that is an approximative GC-based formula, i.e. what a
  ##                   user of the defaults actually gets on genomic windows.
  ##   "rmelting_nn" - method.nn = "san04" forces the nearest-neighbour model
  ##                   regardless of length (verified: the returned enthalpy
  ##                   of ~-1.6e6 cal for a 200 bp duplex is a per-stack sum,
  ##                   which the approximative formula does not produce),
  ##                   giving a like-for-like comparison with TmCalculator.
  ##
  ## Strand-concentration convention: MELTING divides the supplied total
  ## nucleic acid concentration by 4 for non-self-complementary duplexes
  ## (see the "Correction factor" field of the returned Environment), whereas
  ## TmCalculator and Biopython both use dnac1 - dnac2/2. Passing 50 nM here
  ## gives 50/4 = 12.5 nM, matching the 25 - 25/2 = 12.5 nM used by the other
  ## two tools, so the comparison is not confounded by this convention.
  force_nn <- identical(tool, "rmelting_nn")
  t0 <- proc.time()
  tm <- vapply(seqs, function(s) {
    out <- tryCatch({
      a <- list(sequence = s,
                nucleic.acid.conc = 50e-9,    # mol/L; /4 -> 12.5 nM effective
                hybridisation.type = "dnadna",
                Na.conc = 0.05)               # mol/L
      if (force_nn) a$method.nn <- "san04"    # SantaLucia 2004
      do.call(melting, a)
    }, error = function(e) NULL)
    if (is.null(out)) return(NA_real_)
    r <- out$Results                          # a named list, not a data.frame
    j <- grep("melting temperature", tolower(names(r)))
    if (!length(j)) return(NA_real_)
    v <- suppressWarnings(as.numeric(r[[j[1]]]))
    if (length(v) == 0L || is.na(v[1])) NA_real_ else v[1]
  }, numeric(1), USE.NAMES = FALSE)
  elapsed <- (proc.time() - t0)[["elapsed"]]
}

utils::write.csv(data.frame(i = seq_along(tm), Tm = tm), outcsv, row.names = FALSE)
cat(sprintf("BENCH_ELAPSED %.3f\n", elapsed))
cat(sprintf("BENCH_PREP %.3f\n", t_prep))
cat(sprintf("BENCH_N %d\n", length(tm)))
cat(sprintf("BENCH_NA %d\n", sum(is.na(tm))))
