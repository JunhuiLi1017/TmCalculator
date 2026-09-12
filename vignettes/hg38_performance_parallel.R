## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo    = TRUE,
  eval    = FALSE,   # requires BSgenome.Hsapiens.UCSC.hg38 and long runtimes
  message = FALSE,
  warning = FALSE,
  fig.retina = 1,
  dpi        = 72
)

## ----libraries----------------------------------------------------------------
# library(TmCalculator)
# library(BiocParallel)
# 
# pkg <- "BSgenome.Hsapiens.UCSC.hg38"
# suppressPackageStartupMessages(library(pkg, character.only = TRUE))
# genome <- get(pkg, envir = asNamespace(pkg))

## ----measuring----------------------------------------------------------------
# gc(reset = TRUE)                      # reset the "max used" high-water mark
# ## ... run the code to be measured ...
# gc()                                  # read peak from the "max used" column
# 
# ps::ps_memory_info()[["rss"]] / 1e9   # current process resident set, GB

## ----chr21-warmup-------------------------------------------------------------
# chr_len21 <- GenomeInfoDb::seqlengths(genome)[["chr21"]]
# 
# t21 <- system.time({
#   bins21 <- make_genomiccoord(bsgenome = pkg, chromosomes = "chr21",
#                               window = 200L, slide = 200L,
#                               start = 1, end = chr_len21, strand = "+")
#   gr21   <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins21),
#                                    method = "preload_chr")
#   tm21   <- tm_calculate(gr21, method = "tm_nn",
#                          nn_table = "DNA_NN_SantaLucia_2004", Na = 50)
# })
# t21["elapsed"]        # ~10-15 s cold, ~6 s warm, on the test machine
# head(tm21$gr)

## ----chr1-serial--------------------------------------------------------------
# chr_len <- GenomeInfoDb::seqlengths(genome)[["chr1"]]
# 
# t_coord <- system.time({
#   bins <- make_genomiccoord(bsgenome = pkg, chromosomes = "chr1",
#                             window = 200L, slide = 200L,
#                             start = 1, end = chr_len, strand = "+")
# })
# 
# t_extract <- system.time({
#   gr_batch <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
#                                      method = "preload_chr")
# })
# 
# base::gc(reset = TRUE)
# t_tm <- system.time({
#   tm_chr1 <- tm_calculate(gr_batch, method = "tm_nn",
#                           nn_table = "DNA_NN_SantaLucia_2004", Na = 50)
# })
# base::gc()   # "max used" = peak R memory during the Tm step
# 
# rbind(t_coord, t_extract, t_tm)[, "elapsed"]
# ## t_coord t_extract      t_tm
# ##     5.5      15.4      30.8     (freshly booted, otherwise idle machine)

## ----chr1-parallel------------------------------------------------------------
# system.time({
#   tm_serial <- tm_calculate(gr_batch, method = "tm_nn", Na = 50)
# })
# ## elapsed ~ 30 s
# 
# system.time({
#   tm_snow <- tm_calculate(gr_batch, method = "tm_nn", Na = 50,
#                           BPPARAM = SnowParam(workers = 5))
# })
# ## elapsed ~ 55 s on the idle machine (33-80 s across repeated sessions)
# ## -- NEVER faster than serial

## ----genome-parallel----------------------------------------------------------
# chrs <- paste0("chr", c(1:22, "X", "Y"))
# ## Better: sort largest-first for load balance -- see "Choosing the
# ## worker count" below for the one-liner.
# n_workers <- 5   # see "Choosing the worker count" below
# 
# runtime <- system.time({
#   res_list <- bplapply(chrs, function(chr, pkg) {
#     ## SnowParam workers are fresh R processes: load packages HERE,
#     ## not in the manager session.
#     suppressPackageStartupMessages(library(TmCalculator))
#     suppressPackageStartupMessages(library(pkg, character.only = TRUE))
# 
#     genome  <- get(pkg, envir = asNamespace(pkg))
#     chr_len <- GenomeInfoDb::seqlengths(genome)[[chr]]
# 
#     bins <- make_genomiccoord(bsgenome = pkg, chromosomes = chr,
#                               window = 200L, slide = 200L,
#                               start = 1, end = chr_len, strand = "+",
#                               verbose = FALSE)
#     gr  <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
#                                   method = "preload_chr")
#     out <- tm_calculate(gr, method = "tm_nn",
#                         nn_table = "DNA_NN_SantaLucia_2004",
#                         Na = 50)$gr        # serial inside the worker
# 
#     ## Drop sequence columns before returning: Tm/GC are what we keep,
#     ## and this cuts per-chromosome serialization from ~500 MB to a few MB.
#     out$sequence   <- NULL
#     out$complement <- NULL
# 
#     ## Record this worker's peak resident memory (GB) for budgeting.
#     attr(out, "worker_rss_gb") <- ps::ps_memory_info()[["rss"]] / 1e9
#     out
#   }, pkg = pkg, BPPARAM = SnowParam(workers = n_workers))
# 
#   tm_genome <- unlist(GenomicRanges::GRangesList(res_list))
# })
# 
# runtime
# sapply(res_list, attr, "worker_rss_gb")   # per-worker memory check
# length(tm_genome)                         # 14,687,330 windows
# summary(tm_genome$Tm)
# ## Min. 47.2  1st Qu. 68.9  Median 71.5  Mean 72.0  3rd Qu. 74.8  Max. 101.1

## ----strategies---------------------------------------------------------------
# sl   <- GenomeInfoDb::seqlengths(genome)[paste0("chr", c(1:22, "X", "Y"))]
# chrs <- names(sort(sl, decreasing = TRUE))          # largest first
# 
# tasks_chrom <- lapply(chrs, function(ch)
#   list(chr = ch, start = 1L, end = as.integer(sl[[ch]])))
# 
# ## Segment boundaries must be multiples of `slide`, or the window grid
# ## shifts between segments and the result stops matching a whole-
# ## chromosome run. 50 Mb = 250,000 x 200 bp.
# seg_size  <- 50e6
# tasks_seg <- unlist(lapply(names(sl), function(ch) {
#   st <- seq(1, sl[[ch]], by = seg_size)
#   lapply(st, function(s)
#     list(chr = ch, start = as.integer(s),
#          end = as.integer(min(s + seg_size - 1, sl[[ch]]))))
# }), recursive = FALSE)
# 
# length(tasks_chrom)   # 24
# length(tasks_seg)     # 73
# 
# n <- 5
# BPPARAM_static  <- SnowParam(workers = n)
# BPPARAM_dynamic <- SnowParam(workers = n, tasks = length(tasks_chrom))
# BPPARAM_segment <- SnowParam(workers = n, tasks = length(tasks_seg))

## ----sweep-data, eval=TRUE----------------------------------------------------
# The numbers below are the measured sweep. They are held in a data frame
# rather than typed into a markdown table so that the table and the figure
# that follow cannot drift apart: both are rendered from this object.
sweep <- data.frame(
  strategy = rep(c("static", "dynamic", "segment"), each = 6),
  workers  = rep(1:6, 3),
  tasks    = rep(c(24L, 24L, 73L), each = 6),
  wall_s   = c(590.5, 431.9, 328.0, 270.0, 251.0, 244.6,
               589.1, 318.2, 237.4, 209.8, 236.4, 262.0,
               578.5, 320.2, 227.9, 193.9, 180.0, 194.0),
  startup_s = c(0.0, 7.1, 7.9, 8.1, 9.0, 9.7,
                0.0, 7.4, 7.2, 8.2, 8.7, 9.9,
                0.0, 6.9, 7.8, 8.3, 9.0, 9.8),
  work_s   = c(590.5, 613.9, 638.2, 656.5, 730.0, 906.3,
               589.1, 620.7, 683.9, 795.8, 1093.9, 1481.9,
               578.5, 629.1, 665.3, 736.9, 842.6, 1085.9),
  longest_task_s = c(55.2, 57.1, 60.8, 61.9, 70.7, 92.4,
                     53.4, 59.4, 63.6, 69.7, 88.2, 127.5,
                     11.5, 13.8, 15.6, 15.7, 17.5, 22.6),
  rss_gb   = c(4.75, 4.06, 4.06, 3.85, 3.70, 3.35,
               4.74, 3.86, 3.38, 2.90, 2.67, 2.15,
               4.87, 2.58, 2.58, 2.34, 1.86, 1.77),
  stringsAsFactors = FALSE)

# Derived rather than transcribed. Speedup is the one-worker wall clock of
# the SAME strategy over the configuration's wall clock: segmenting changes
# how much total work there is, so a ratio formed against a configuration's
# own summed task times would credit a strategy for doing more work and the
# three could not be compared. Efficiency follows from that speedup, not
# from the task times.
serial <- with(sweep, tapply(wall_s[workers == 1], strategy[workers == 1], I))
sweep$speedup    <- serial[sweep$strategy] / sweep$wall_s
sweep$efficiency <- sweep$speedup / sweep$workers
serial_work <- with(sweep, tapply(work_s[workers == 1], strategy[workers == 1], I))
sweep$work_ratio <- sweep$work_s / serial_work[sweep$strategy]
sweep$strategy <- factor(sweep$strategy, levels = c("static", "dynamic", "segment"))
sweep <- sweep[order(sweep$strategy, sweep$workers), ]

## ----sweep-table, eval=TRUE---------------------------------------------------
tab <- data.frame(
  Strategy   = ifelse(duplicated(sweep$strategy), "", as.character(sweep$strategy)),
  Workers    = sweep$workers,
  Tasks      = sweep$tasks,
  `Wall (s)` = round(sweep$wall_s, 1),
  `Start-up (s)` = round(sweep$startup_s, 1),
  `Total task time (s)` = round(sweep$work_s, 1),
  `Work ratio` = round(sweep$work_ratio, 2),
  `Longest task (s)` = round(sweep$longest_task_s, 1),
  Speedup    = round(sweep$speedup, 2),
  Efficiency = round(sweep$efficiency, 2),
  `Peak RSS per worker (GB)` = round(sweep$rss_gb, 2),
  check.names = FALSE, stringsAsFactors = FALSE)
knitr::kable(tab, row.names = FALSE,
             caption = "Task-partitioning strategy and worker count.")

## ----sweep-figure, eval=TRUE, fig.width=10.5, fig.height=3.9, fig.cap="Parallel performance under three task-partitioning strategies. (A) Speedup against worker count, relative to the one-worker run of the same strategy; the dashed line marks linear speedup. (B) Total task time, summed over all tasks and measured inside the workers; values above the dashed serial reference indicate that a configuration performed more work than the serial run rather than dividing the same work among more processes. (C) Peak resident set size of the heaviest worker."----
# Three panels, one argument. Panel A shows that every strategy turns over
# before the cores run out, which invites the usual explanation of a ragged
# schedule. Panel B rules that out: the total task time grows with the worker
# count, so the later configurations are doing more work, not dividing a
# fixed amount badly. Panel C gives the reason and the remedy together, since
# the strategy that holds its memory down is the one whose work grows least.
pal  <- c(static = "#1B5E9C", dynamic = "#C0392B", segment = "#5C6B73")
pchs <- c(static = 16,        dynamic = 17,        segment = 15)
lv   <- levels(sweep$strategy)

op <- par(mfrow = c(1, 3), mar = c(4.4, 4.5, 2.2, 0.8), las = 1,
          cex = 0.8, mgp = c(2.8, 0.7, 0))

series <- function(col, ylab, ref = NULL, ref_lab = NULL, ymax = NULL) {
  plot(NA, xlim = range(sweep$workers),
       ylim = c(0, if (is.null(ymax)) max(sweep[[col]]) * 1.05 else ymax),
       bty = "n", xaxt = "n", xlab = "Workers", ylab = ylab)
  axis(1, at = sort(unique(sweep$workers)))
  if (!is.null(ref)) {
    abline(h = ref, lty = 2, col = "grey55")
    text(max(sweep$workers), ref, ref_lab, adj = c(1.1, -0.5),
         cex = 0.75, col = "grey40")
  }
  for (s in lv) {
    d <- sweep[sweep$strategy == s, ]
    lines(d$workers, d[[col]], col = pal[s], lwd = 2)
    points(d$workers, d[[col]], col = pal[s], pch = pchs[s])
  }
}

series("speedup", "Speedup", ymax = max(sweep$workers) * 1.02)
abline(a = 0, b = 1, lty = 2, col = "grey55")
text(max(sweep$workers), max(sweep$workers), "linear",
     adj = c(1.1, -0.4), cex = 0.75, col = "grey40")
for (s in lv) {                      # redraw over the reference line
  d <- sweep[sweep$strategy == s, ]
  lines(d$workers, d$speedup, col = pal[s], lwd = 2)
  points(d$workers, d$speedup, col = pal[s], pch = pchs[s])
}
mtext("A", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)
legend("topleft", bty = "n", legend = lv, col = pal[lv], pch = pchs[lv],
       lwd = 2, seg.len = 1.3, cex = 0.85)

series("work_s", "Total task time (s)",
       ref = median(sweep$work_s[sweep$workers == 1]), ref_lab = "serial")
mtext("B", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)

series("rss_gb", "Peak resident set size per worker (GB)")
mtext("C", side = 3, adj = 0, font = 2, line = 0.7, cex = 1.05)

par(op)

## ----per-chrom----------------------------------------------------------------
# tk <- utils::read.csv("bench_parallel_strategy_tasks.csv")
# tk <- tk[tk$strategy == "dynamic", ]        # one task per chromosome
# 
# ## Median over repetitions, one column per worker count.
# w <- stats::reshape(
#   stats::aggregate(secs ~ chr + n_workers, data = tk, FUN = stats::median),
#   idvar = "chr", timevar = "n_workers", direction = "wide")
# names(w) <- sub("^secs\\.", "w", names(w))
# 
# w$inflation <- w$w6 / w$w1                  # 6 workers vs serial
# w <- w[order(-w$w1), ]
# w

## ----segment-parallel---------------------------------------------------------
# chr_len <- GenomeInfoDb::seqlengths(genome)[["chr1"]]
# 
# ## Segment length must be a multiple of `slide` so the window grid stays
# ## aligned across segment boundaries (50 Mb = 250,000 x 200).
# seg_starts <- seq(1, chr_len, by = 50e6)
# seg <- data.frame(start = seg_starts,
#                   end   = pmin(seg_starts + 50e6 - 1, chr_len))
# 
# t_seg <- system.time({
#   res_seg <- bplapply(seq_len(nrow(seg)), function(i, seg, pkg) {
#     suppressPackageStartupMessages(library(TmCalculator))
#     suppressPackageStartupMessages(library(pkg, character.only = TRUE))
#     bins <- make_genomiccoord(bsgenome = pkg, chromosomes = "chr1",
#                               window = 200L, slide = 200L,
#                               start = seg$start[i], end = seg$end[i],
#                               strand = "+", trim_N = "none",
#                               verbose = FALSE)
#     gr <- to_genomic_ranges_fast(list(pkg_name = pkg, seq = bins),
#                                  method = "preload_chr")
#     out <- tm_calculate(gr, method = "tm_nn",
#                         nn_table = "DNA_NN_SantaLucia_2004", Na = 50)$gr
#     out$sequence <- NULL; out$complement <- NULL
#     out
#   }, seg = seg, pkg = pkg, BPPARAM = SnowParam(workers = 5))
# 
#   tm_chr1_seg <- sort(unlist(GenomicRanges::GRangesList(res_seg)))
# })
# 
# t_seg["elapsed"]      # 38 s measured, vs ~52 s for the serial FULL pipeline
# length(tm_chr1_seg)   # 1,152,300

## ----workers------------------------------------------------------------------
# mem_gb    <- 16                                     # your machine
# cores     <- parallel::detectCores(logical = FALSE) # physical cores
# per_worker_gb <- 2                                  # 50 Mb segment task
# n_workers <- min(cores - 1L, floor((mem_gb - 4) / per_worker_gb))
# n_workers
# ## 16 GB / 6 cores  -> 5 workers, which is what the sweep found fastest
# ## 32 GB / 8 cores  -> 7 workers
# ## 64 GB / 10 cores -> 9 workers

## ----sort-chrs----------------------------------------------------------------
# sl   <- GenomeInfoDb::seqlengths(genome)
# chrs <- names(sort(sl[paste0("chr", c(1:22, "X", "Y"))], decreasing = TRUE))

## ----sessioninfo, eval=TRUE---------------------------------------------------
sessionInfo()

