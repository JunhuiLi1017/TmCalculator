#!/usr/bin/env Rscript
# ===========================================================================
# find_isoGC_pairs.R -- windows with identical GC content but different Tm
#
# Supports the additional Figure 2 panel requested in review: a case in which
# two genomic windows of the same length and the same GC content melt at
# different temperatures, which is the effect a nearest-neighbor model
# captures and a GC-content formula cannot.
#
# Run after the E. coli case study, with `tm_ASM584v2` in the session:
#   source(system.file("scripts", "find_isoGC_pairs.R", package = "TmCalculator"))
# ===========================================================================

suppressPackageStartupMessages({
  library(TmCalculator)
  library(GenomicRanges)
})

stopifnot(exists("tm_ASM584v2"))

gr <- tm_ASM584v2$gr
d  <- data.frame(
  idx   = seq_along(gr),
  chr   = as.character(seqnames(gr)),
  start = start(gr),
  end   = end(gr),
  Tm    = gr$Tm,
  GC    = gr$GC,
  seq   = as.character(gr$sequence),
  stringsAsFactors = FALSE
)
d <- d[is.finite(d$Tm) & is.finite(d$GC), ]
d <- d[nchar(d$seq) == max(nchar(d$seq)), ]   # keep full-width windows only

## Within each exact GC value, the extreme pair is the one to show: same
## composition, maximal difference in stacking energy.
best <- do.call(rbind, lapply(split(d, round(d$GC, 6)), function(g) {
  if (nrow(g) < 2) return(NULL)
  lo <- g[which.min(g$Tm), ]
  hi <- g[which.max(g$Tm), ]
  data.frame(
    GC        = lo$GC,
    n_windows = nrow(g),
    dTm       = hi$Tm - lo$Tm,
    lo_locus  = sprintf("%s:%d-%d", lo$chr, lo$start, lo$end),
    lo_Tm     = lo$Tm,
    hi_locus  = sprintf("%s:%d-%d", hi$chr, hi$start, hi$end),
    hi_Tm     = hi$Tm,
    lo_seq    = lo$seq,
    hi_seq    = hi$seq,
    stringsAsFactors = FALSE
  )
}))
best <- best[order(-best$dTm), ]

## Restrict to GC values that are well represented, so the pair is typical
## of the genome rather than an isolated outlier.
common <- best[best$n_windows >= 100, ]

cat("\nTop iso-GC pairs (GC identical, Tm maximally different)\n")
cat("-------------------------------------------------------\n")
print(format(head(common[, c("GC", "n_windows", "dTm",
                             "lo_locus", "lo_Tm", "hi_locus", "hi_Tm")], 10),
             digits = 4), row.names = FALSE)

cat("\nAt the median GC of the genome:\n")
med <- common[which.min(abs(common$GC - median(d$GC, na.rm = TRUE))), ]
cat(sprintf("  GC = %.3f, %d windows, dTm = %.2f C\n",
            med$GC, med$n_windows, med$dTm))
cat(sprintf("  low  %s  Tm = %.2f C\n", med$lo_locus, med$lo_Tm))
cat(sprintf("       %s\n", med$lo_seq))
cat(sprintf("  high %s  Tm = %.2f C\n", med$hi_locus, med$hi_Tm))
cat(sprintf("       %s\n", med$hi_seq))

## Dinucleotide composition explains the difference; report the steps that
## differ most, since that is what the panel should annotate.
steps <- function(x) {
  k <- substring(x, 1:(nchar(x) - 1L), 2:nchar(x))
  table(factor(k, levels = c(outer(c("A","C","G","T"), c("A","C","G","T"), paste0))))
}
diff_steps <- steps(med$hi_seq) - steps(med$lo_seq)
cat("\nDinucleotide steps enriched in the high-Tm window (top 6):\n")
print(head(sort(diff_steps, decreasing = TRUE), 6))
cat("Depleted (top 6):\n")
print(head(sort(diff_steps), 6))

utils::write.csv(common, "isoGC_pairs.csv", row.names = FALSE)
cat("\nFull table written to isoGC_pairs.csv\n")

invisible(list(all = best, common = common, chosen = med))

getwd()

png("Figure2.png", width = 180, height = 180, units = "mm", res = 600)   # 或
tiff("Figure2.tif", width = 180, height = 180, units = "mm", res = 600, compression = "lzw")




dev.off()

