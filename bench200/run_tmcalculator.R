
## Four timestamps, not two. Reading the input and writing the results are
## artefacts of running each tool as a separate process, and they scale with
## the input, so folding them into "start-up" makes that quantity grow with n
## and describe nothing.
tA   <- proc.time()[["elapsed"]]
suppressPackageStartupMessages(library(TmCalculator))
tB   <- proc.time()[["elapsed"]]
a    <- commandArgs(trailingOnly = TRUE)
seqs <- readLines(a[1])
tC   <- proc.time()[["elapsed"]]
res  <- tm_calculate(seqs, method = "tm_nn",
                     nn_table    = "DNA_NN_SantaLucia_2004",
                     salt_method = "SantaLucia1996",
                     Na = 50, dnac_high = 25, dnac_low = 25,
                     self_comp = FALSE)
tD   <- proc.time()[["elapsed"]]
writeLines(format(res$gr$Tm, digits = 10), a[2])
tE   <- proc.time()[["elapsed"]]
## timing goes to a file, not stdout: rmelting drives the Java engine through
## rJava, which redirects stdout, so a printed value can be swallowed
tim <- sprintf(c("LOAD_SECONDS %.4f", "READ_SECONDS %.4f",
                 "COMPUTE_SECONDS %.4f", "WRITE_SECONDS %.4f"),
               c(tB - tA, tC - tB, tD - tC, tE - tD))
if (length(a) >= 3L) writeLines(tim, a[3])
cat(tim, sep = "\n")

