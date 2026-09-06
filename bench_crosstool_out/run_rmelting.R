
suppressPackageStartupMessages(library(rmelting))
a    <- commandArgs(trailingOnly = TRUE)
seqs <- readLines(a[1])
t0   <- proc.time()[["elapsed"]]
tm <- vapply(seqs, function(s) {
  r <- rmelting::melting(sequence           = s,
                         nucleic.acid.conc  = 1.25e-08,
                         hybridisation.type = "dnadna",
                         Na.conc            = 0.05,
                         method.nn          = "san04",
                         correction.ion     = "san96")
  as.numeric(r$Results[["Melting temperature (C)"]])
}, numeric(1), USE.NAMES = FALSE)
t1 <- proc.time()[["elapsed"]]
writeLines(format(tm, digits = 10), a[2])
cat(sprintf("COMPUTE_SECONDS %.4f\n", t1 - t0))

