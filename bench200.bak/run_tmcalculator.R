
suppressPackageStartupMessages(library(TmCalculator))
a    <- commandArgs(trailingOnly = TRUE)
seqs <- readLines(a[1])
t0   <- proc.time()[["elapsed"]]
res  <- tm_calculate(seqs, method = "tm_nn",
                     nn_table    = "DNA_NN_SantaLucia_2004",
                     salt_method = "SantaLucia1996",
                     Na = 50, dnac_high = 25, dnac_low = 25,
                     self_comp = FALSE,
                     BPPARAM = BiocParallel::SerialParam())
t1 <- proc.time()[["elapsed"]]
writeLines(format(res$gr$Tm, digits = 10), a[2])
## timing goes to a file, not stdout: rmelting drives the Java engine through
## rJava, which redirects stdout, so a printed value can be swallowed
if (length(a) >= 3L) writeLines(sprintf("COMPUTE_SECONDS %.4f", t1 - t0), a[3])
cat(sprintf("COMPUTE_SECONDS %.4f\n", t1 - t0))

