
tA   <- proc.time()[["elapsed"]]
suppressPackageStartupMessages(library(rmelting))
tB   <- proc.time()[["elapsed"]]
a    <- commandArgs(trailingOnly = TRUE)
seqs <- readLines(a[1])
tC   <- proc.time()[["elapsed"]]
tm <- vapply(seqs, function(s) {
  r <- rmelting::melting(sequence           = s,
                         nucleic.acid.conc  = 1.25e-08,
                         hybridisation.type = "dnadna",
                         Na.conc            = 0.05,
                         method.nn          = "san04",
                         correction.ion     = "san96",
                         size.threshold     = 201)
  as.numeric(r$Results[["Melting temperature (C)"]])
}, numeric(1), USE.NAMES = FALSE)
tD   <- proc.time()[["elapsed"]]
writeLines(format(tm, digits = 10), a[2])
tE   <- proc.time()[["elapsed"]]
tim <- sprintf(c("LOAD_SECONDS %.4f", "READ_SECONDS %.4f",
                 "COMPUTE_SECONDS %.4f", "WRITE_SECONDS %.4f"),
               c(tB - tA, tC - tB, tD - tC, tE - tD))
if (length(a) >= 3L) writeLines(tim, a[3])
cat(tim, sep = "\n")

