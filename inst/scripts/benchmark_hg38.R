# ==========================================================================
# benchmark_hg38.R -- single controlled run producing every number cited in
# vignettes/hg38_performance_parallel.Rmd
#
# HOW TO RUN (important):
#   1. Reboot or make sure swap is near zero and no stray R workers exist
#      (macOS: Activity Monitor -> Memory; `ps aux | grep exec/R`).
#   2. Start a FRESH R session. Do NOT use devtools::load_all().
#   3. source(system.file("scripts", "benchmark_hg38.R",
#                         package = "TmCalculator"))
#      or source the file from the repo after devtools::install().
#
# Output: a data.frame `bench` printed as markdown, and benchmark_hg38.csv
# in the working directory. Total runtime ~12-15 min with run_genome=TRUE,
# ~6 min without.
# ==========================================================================

run_genome <- FALSE          # TRUE = also run the whole-genome parallel case
n_workers  <- 3L             # match your RAM: min(physical cores - 1, (GB-4)/3)

## -- 0. Environment guards -------------------------------------------------
stopifnot(
  "Run in a fresh session, not load_all()" =
    !"devtools_shims" %in% search(),
  "Session restored a saved workspace; start R with --vanilla or disable .RData restore" =
    length(ls(envir = .GlobalEnv)) <= 2L
)
suppressPackageStartupMessages({
  library(TmCalculator)
  library(BiocParallel)
})
pkg <- "BSgenome.Hsapiens.UCSC.hg38"
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
genome <- get(pkg, envir = asNamespace(pkg))
sl     <- GenomeInfoDb::seqlengths(genome)

has_ps <- requireNamespace("ps", quietly = TRUE)
rss_gb <- function() if (has_ps) ps::ps_memory_info()[["rss"]] / 1e9 else NA_real_

bench <- data.frame(section = character(), what = character(),
                    elapsed_s = numeric(), peak_r_mb = numeric(),
                    rss_gb = numeric(), n_windows = integer())
note <- function(section, what, t, peak_mb = NA_real_, n = NA_integer_) {
  bench <<- rbind(bench, data.frame(
    section = section, what = what, elapsed_s = round(unname(t["elapsed"]), 1),
    peak_r_mb = round(peak_mb), rss_gb = round(rss_gb(), 2), n_windows = n))
  message(sprintf("[bench] %-18s %-28s %7.1f s", section, what,
                  unname(t["elapsed"])))
}
peak_mb <- function(g) sum(g[, "max used"] * c(56, 8)) / 1e6  # rough MB

run_chr <- function(chr, label, warm = FALSE) {
  chr_len <- sl[[chr]]
  base::gc(reset = TRUE)
  t1 <- system.time(
    bins <- make_genomiccoord(bsgenome = pkg, chromosomes = chr,
                              window = 200L, slide = 200L,
                              start = 1, end = chr_len, strand = "+",
                              verbose = FALSE))
  note(label, "make_genomiccoord", t1)
  t2 <- system.time(
    gr <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
                                 method = "preload_chr"))
  note(label, "to_genomic_ranges_fast", t2)
  base::gc(reset = TRUE)
  t3 <- system.time(
    tm <- suppressWarnings(
      tm_calculate(gr, method = "tm_nn",
                   nn_table = "DNA_NN_SantaLucia_2004", Na = 50)))
  note(label, "tm_calculate (serial)", t3,
       peak_mb = peak_mb(base::gc()), n = length(tm$gr))
  invisible(list(gr = gr, tm = tm))
}

## -- 1. chr21 warm-up: cold then warm --------------------------------------
res21  <- run_chr("chr21", "chr21 cold")
res21b <- run_chr("chr21", "chr21 warm")
rm(res21, res21b); base::gc()

## -- 2. chr1: full pipeline (genome file now cached = 'warm disk') ---------
res1 <- run_chr("chr1", "chr1")

## -- 3. chr1: within-call BPPARAM, same session, back to back --------------
gr1 <- res1$gr
t_serial <- system.time(
  tm_s <- suppressWarnings(tm_calculate(gr1, method = "tm_nn", Na = 50)))
note("chr1 BPPARAM", "SerialParam (default)", t_serial)
t_snow <- system.time(
  tm_p <- suppressWarnings(tm_calculate(gr1, method = "tm_nn", Na = 50,
                                        BPPARAM = SnowParam(n_workers))))
note("chr1 BPPARAM", sprintf("SnowParam(%d)", n_workers), t_snow)
stopifnot(all.equal(tm_s$gr$Tm, tm_p$gr$Tm))   # results identical
rm(res1, gr1, tm_s, tm_p); base::gc()

## -- 4. chromosome-level parallelism: chr10-chr19 --------------------------
chr_worker <- function(chr, pkg) {
  suppressPackageStartupMessages(library(TmCalculator))
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  genome  <- get(pkg, envir = asNamespace(pkg))
  chr_len <- GenomeInfoDb::seqlengths(genome)[[chr]]
  bins <- make_genomiccoord(bsgenome = pkg, chromosomes = chr,
                            window = 200L, slide = 200L,
                            start = 1, end = chr_len, strand = "+",
                            verbose = FALSE)
  gr  <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
                                method = "preload_chr")
  out <- suppressWarnings(
    tm_calculate(gr, method = "tm_nn",
                 nn_table = "DNA_NN_SantaLucia_2004", Na = 50))$gr
  out$sequence <- NULL; out$complement <- NULL
  attr(out, "rss_gb") <- if (requireNamespace("ps", quietly = TRUE))
    ps::ps_memory_info()[["rss"]] / 1e9 else NA_real_
  out
}
chrs10 <- names(sort(sl[paste0("chr", 10:19)], decreasing = TRUE))

t_ser10 <- system.time(
  res_ser <- lapply(chrs10, chr_worker, pkg = pkg))
note("chr10-19", "serial (lapply)", t_ser10,
     n = sum(vapply(res_ser, length, 1L)))
rm(res_ser); base::gc()

t_par10 <- system.time(
  res_par <- bplapply(chrs10, chr_worker, pkg = pkg,
                      BPPARAM = SnowParam(n_workers)))
note("chr10-19", sprintf("SnowParam(%d)", n_workers), t_par10,
     n = sum(vapply(res_par, length, 1L)))
message("worker RSS (GB): ",
        paste(round(vapply(res_par, attr, 1.0, "rss_gb"), 2), collapse = " "))
rm(res_par); base::gc()

## -- 5. whole genome (optional) --------------------------------------------
if (run_genome) {
  chrs_all <- names(sort(sl[paste0("chr", c(1:22, "X", "Y"))],
                         decreasing = TRUE))
  t_gen <- system.time(
    res_gen <- bplapply(chrs_all, chr_worker, pkg = pkg,
                        BPPARAM = SnowParam(n_workers)))
  tm_genome <- unlist(GenomicRanges::GRangesList(res_gen))
  note("genome", sprintf("SnowParam(%d)", n_workers), t_gen,
       n = length(tm_genome))
  print(summary(tm_genome$Tm))
}

## -- 6. Report -------------------------------------------------------------
utils::write.csv(bench, "benchmark_hg38.csv", row.names = FALSE)
cat("\n\n")
cat(sprintf("| %s | %s | %s | %s | %s |\n",
            "section", "step", "elapsed (s)", "peak R mem (MB)", "RSS (GB)"))
cat("|---|---|---|---|---|\n")
for (i in seq_len(nrow(bench))) {
  cat(sprintf("| %s | %s | %.1f | %s | %s |\n",
              bench$section[i], bench$what[i], bench$elapsed_s[i],
              ifelse(is.na(bench$peak_r_mb[i]), "", bench$peak_r_mb[i]),
              ifelse(is.na(bench$rss_gb[i]), "", bench$rss_gb[i])))
}
cat("\nSaved to benchmark_hg38.csv -- paste the table above back for the vignette.\n")
