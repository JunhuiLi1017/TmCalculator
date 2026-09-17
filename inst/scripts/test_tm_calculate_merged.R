#!/usr/bin/env Rscript
# Exercise the merged tm_calculate() before it goes into the package.
#
#   cd .../TmCalculator/R && Rscript test_tm_calculate_merged.R
#
# Needs TmCalculator INSTALLED, not just load_all'ed: the parallel section
# ships task closures to fresh R processes and those load the installed
# package. The two files under test, tm_calculate.R and tm_source.R, are
# sourced here so they can be checked before installing.
#
# Everything runs on half a megabyte of chr21 and finishes in about a minute.

suppressPackageStartupMessages({
  library(TmCalculator)
  library(GenomicRanges)
})
## -- find the two files under test ------------------------------------------
# Tried in order: an explicit TM_SRC, the directory this script sits in
# resolved back to the package's R/, the working directory, and R/ below it.
# The script lives in inst/scripts/, so ../../R is the checkout's R/ whether
# it is run from the package root, from inst/scripts/, or by absolute path.
NEEDED <- c("tm_source.R", "tm_calculate.R")
here <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- if (length(here)) dirname(normalizePath(here[1])) else NA_character_
cand <- c(Sys.getenv("TM_SRC"),
          if (!is.na(here)) file.path(here, "..", "..", "R"),
          if (!is.na(here)) file.path(here, "..", "..", "..", "R"),
          ".", "R", "../R")
cand <- unique(cand[nzchar(cand)])
SRC  <- NULL
for (d in cand)
  if (all(file.exists(file.path(d, NEEDED)))) { SRC <- d; break }
if (is.null(SRC))
  stop("Cannot find ", paste(NEEDED, collapse = " and "), ".\n",
       "  Run this from the package root, or set TM_SRC to its R/ directory:\n",
       "    TM_SRC=/path/to/TmCalculator/R Rscript ",
       if (!is.na(here)) file.path(here, "test_tm_calculate_merged.R")
       else "test_tm_calculate_merged.R", "\n",
       "  Looked in: ", paste(normalizePath(cand, mustWork = FALSE),
                              collapse = "\n             "))
cat("sources under test : ", normalizePath(SRC), "\n", sep = "")
for (f in NEEDED) source(file.path(SRC, f))

PKG <- "BSgenome.Hsapiens.UCSC.hg38"
if (!requireNamespace(PKG, quietly = TRUE))
  stop("Install ", PKG, " first: BiocManager::install(\"", PKG, "\")")

# Past the telomeric N run, so every task has real sequence in it.
REG  <- "chr21:10000001-10500000"
ARGS <- list(method = "tm_nn", nn_table = "DNA_NN_SantaLucia_2004", Na = 50)

ok <- 0L; bad <- 0L
check <- function(label, expr) {
  res <- tryCatch(isTRUE(expr),
                  error = function(e) structure(FALSE, msg = conditionMessage(e)))
  if (isTRUE(res)) { ok <<- ok + 1L; cat(sprintf("  PASS  %s\n", label)) }
  else { bad <<- bad + 1L
         cat(sprintf("  FAIL  %s%s\n", label,
                     if (!is.null(attr(res, "msg")))
                       paste0("  [", attr(res, "msg"), "]") else "")) }
}
err   <- function(e) inherits(try(e, silent = TRUE), "try-error")
warns <- function(e) tryCatch({ force(e); FALSE }, warning = function(w) TRUE)
same  <- function(a, b)
  length(a) == length(b) &&
  identical(as.character(seqnames(a)), as.character(seqnames(b))) &&
  identical(start(a), start(b)) &&
  isTRUE(all.equal(a$Tm, b$Tm)) && isTRUE(all.equal(a$GC, b$GC))
prof  <- function(...) do.call(tm_calculate, c(list(...), ARGS))$gr
# Assignments outside check() used to halt the run at the first failure, so a
# single bug hid every check after it. safe() turns an error into a value the
# later checks can fail on individually.
safe  <- function(e) tryCatch(e, error = function(err) {
  cat(sprintf("  ERROR %s\n", conditionMessage(err))); err })

## -- 1. The direct route is untouched ---------------------------------------
cat("\n1. Sequences in, one Tm each: the behaviour that must not change\n")

seqs <- c("ATGCGATGCGAAGGCGATGGCGTGTAGAATAGATCACATACTGCATAGCTGATC",
          "ATGCGATGCGCCCGGAGATAGAAGGCGTAGATACAGATCAGTAGCACCTTGAGAC")
old <- safe(tm_calculate(seqs))
check("returns a TmCalculator object", inherits(old, "TmCalculator"))
check("regression values unchanged",
      isTRUE(all.equal(round(old$gr$Tm, 5), c(67.06562, 69.64434), tolerance = 0.1)))
check("sequence columns kept for supplied sequences",
      all(c("sequence", "complement") %in% names(mcols(old$gr))))

## -- 2. A genome, and the invariant that matters ----------------------------
cat("\n2. A BSgenome source: regions, tiling, segmenting\n")

one <- safe(prof(PKG, regions = REG, window = 200, slide = 200, unit = "region",
            verbose = FALSE))
bins <- safe(make_genomiccoord(bsgenome = PKG, chromosomes = "chr21",
                               window = 200, slide = 200, start = 10000001,
                               end = 10500000, strand = "+", trim_N = "none",
                               verbose = FALSE))
manual <- safe(do.call(tm_calculate,
                       c(list(input_seq = to_genomic_ranges_fast(
                         list(pkg_name = PKG, seq = bins), method = "preload_chr")),
                         ARGS))$gr)
check("agrees with the hand-written three-step route", same(one, manual))
check("sequence columns dropped for a genome source",
      !any(c("sequence", "complement") %in% names(mcols(one))))

seg <- safe(prof(PKG, regions = REG, window = 200, slide = 200, unit = "segment",
            segment_size = 200e3, verbose = FALSE))
check("segmenting does not change the windows", same(seg, one))
# A segment size that is not a multiple of slide is rounded down, which is
# the case most likely to go wrong silently.
seg2 <- safe(prof(PKG, regions = REG, window = 200, slide = 200, unit = "segment",
                  segment_size = 123456, verbose = FALSE))
check("ragged segment_size still aligns to the grid", same(seg2, one))

## -- 2b. The sequence really is the sequence at those coordinates ----------
# Nothing else here would catch a profile that is correct in shape and
# wrong in position. Every other check compares one run of this code with
# another, and a coordinate error inside window construction or sequence
# extraction moves both alike: the segmented and unsegmented runs would
# agree with each other and both be shifted along the chromosome. This is
# the only check with an outside witness, so it stays whatever the
# extraction path does internally.
cat("\n2b. Window sequences against getSeq at the same coordinates\n")

gobj <- safe(get(sub("^BSgenome\\.([^.]+)\\..*$", "\\1", PKG),
                 envir = asNamespace(PKG)))
truth <- function(gr)
  as.character(Biostrings::getSeq(gobj, as.character(seqnames(gr)),
                                  start = start(gr), end = end(gr)))

# Two regions with different starts: an off-by-origin error survives one of
# these only when the region happens to begin at 1.
for (rg in c("chr21:10200001-10203000", "chr21:31415927-31418927")) {
  w <- safe(prof(PKG, regions = rg, window = 200, slide = 200,
                 keep_sequence = TRUE, verbose = FALSE))
  check(sprintf("%s: sequences match the genome", rg),
        identical(as.character(w$sequence), truth(w)))
  check(sprintf("%s: first window starts where asked", rg),
        start(w)[1] == as.numeric(sub(".*:(\\d+)-.*", "\\1", rg)))
}

# The same region reached as one task and as three must agree, and both must
# agree with the genome: agreeing only with each other is what a shifted
# profile also does.
s1 <- safe(prof(PKG, regions = "chr21:10200001-10203000", window = 200,
                slide = 200, unit = "region", keep_sequence = TRUE,
                verbose = FALSE))
s3 <- safe(prof(PKG, regions = "chr21:10200001-10203000", window = 200,
                slide = 200, unit = "segment", segment_size = 1000,
                keep_sequence = TRUE, verbose = FALSE))
check("one task and three give the same sequences", same(s1, s3) &&
        identical(as.character(s1$sequence), as.character(s3$sequence)))
check("three tasks still match the genome",
      identical(as.character(s3$sequence), truth(s3)))

## -- 3. regions, in every form ----------------------------------------------
cat("\n3. Region forms\n")

small <- c("chr21:10000001-10100000", "chr22:15000001-15100000")
a <- safe(prof(PKG, regions = small,      unit = "region", window = 200, verbose = FALSE))
b <- safe(prof(PKG, regions = rev(small), unit = "region", window = 200, verbose = FALSE))
check("two regions on two chromosomes", length(unique(seqnames(a))) == 2L)
check("result order is genomic, not the order typed", same(a, b))

c1 <- safe(prof(PKG, regions = "chr21:10000001-10100000", window = 200, verbose = FALSE))
c2 <- safe(prof(PKG, regions = "chr21:10,000,001-10,100,000", window = 200, verbose = FALSE))
c3 <- safe(prof(PKG, regions = GRanges("chr21", IRanges(10000001, 10100000)),
                window = 200, verbose = FALSE))
check("commas in coordinates parse", same(c1, c2))
check("a GRanges selects the same region", same(c1, c3))

src <- .tm_source(PKG)
check("regions = 21 and \"chr21\" agree",
      identical(.tm_regions(21, src), .tm_regions("chr21", src)))
check("a bare chromosome is whole, a sub-region is not",
      all(.tm_regions("chr21", src)$whole) &&
        isFALSE(.tm_regions("chr21:1-1000", src)$whole))
check("the default is standardChromosomes, chrM included",
      "chrM" %in% .tm_regions(NULL, src)$name)

cat("\n   Errors and warnings\n")
check("unknown chromosome errors",    err(.tm_regions("chr99", src)))
check("end beyond chromosome errors", err(.tm_regions("chr21:1-999999999", src)))
check("start after end errors",       err(.tm_regions("chr21:500-100", src)))
check("malformed region errors",      err(.tm_regions("chr21:100", src)))
check("overlapping regions warn",
      warns(.tm_regions(c("chr21:1-2000", "chr21:1500-3000"), src)))
check("adjacent regions do not warn",
      !warns(.tm_regions(c("chr21:1-2000", "chr21:2001-3000"), src)))
check("window = NULL on a whole chromosome is refused",
      err(tm_calculate(PKG, regions = "chr21", verbose = FALSE)))

## -- 4. Sequences as a source, with regions ---------------------------------
cat("\n4. Sequences as a source\n")

set.seed(1)
oligos <- vapply(seq_len(40), function(i)
  paste(sample(c("A", "C", "G", "T"), 60, TRUE), collapse = ""), character(1))
named <- oligos; names(named) <- sprintf("oligo_%02d", seq_along(named))

ref <- safe(tm_calculate(oligos)$gr)
v1 <- safe(prof(named, verbose = FALSE))
check("one row per sequence", length(v1) == length(named))
check("Tm unchanged by the staging round trip",
      isTRUE(all.equal(sort(v1$Tm), sort(ref$Tm))))
check("names become seqnames", all(as.character(seqnames(v1)) %in% names(named)))
check("no temporary file left behind",
      !length(list.files(tempdir(), pattern = "^tm_.*\\.fa$")))

v2 <- safe(prof(named, regions = "oligo_03", verbose = FALSE))
check("regions selects by the caller's name", length(v2) == 1L)
v3 <- safe(prof(unname(oligos), regions = "3:1-20", verbose = FALSE))
check("an unnamed vector is addressed by position",
      length(v3) == 1L && width(v3) == 20L)
check("\"3:1-20\" is the third sequence's first 20 bases",
      isTRUE(all.equal(v3$Tm, tm_calculate(substr(oligos[3], 1, 20))$gr$Tm)))

## -- 4b. A GRanges as the source: regions means overlap ---------------------
cat("\n4b. A GRanges source, selected by overlap\n")

gsrc <- GRanges(c("chrA", "chrA", "chrB"),
                IRanges(c(1, 1001, 1), width = c(60, 60, 60)),
                sequence = unname(oligos[1:3]))
g0 <- safe(prof(gsrc, verbose = FALSE))
check("all ranges when regions is NULL", length(g0) == 3L)

check("a seqname takes every range on it",
      length(prof(gsrc, regions = "chrA", verbose = FALSE)) == 2L)
# Staging turns each range into a FASTA record that starts at 1, so without
# an offset the second range would come back at 1-60 instead of 1001-1060.
check("a staged GRanges keeps its own coordinates",
      identical(start(prof(gsrc, regions = "chrA", verbose = FALSE)),
                c(1L, 1001L)))
check("a sequence-only GRanges needs no complement column",
      inherits(safe(tm_calculate(gsrc, method = "tm_nn")), "TmCalculator"))
check("an interval takes only what it overlaps",
      length(prof(gsrc, regions = "chrA:1-100", verbose = FALSE)) == 1L)
check("a one-base overlap still counts",
      length(prof(gsrc, regions = "chrA:1060-2000", verbose = FALSE)) == 1L)
check("a gap between ranges selects nothing",
      err(prof(gsrc, regions = "chrA:200-900", verbose = FALSE)))
check("a GRanges query overlaps the same way",
      length(prof(gsrc, regions = GRanges("chrA", IRanges(1, 100)),
                  verbose = FALSE)) == 1L)
check("several regions union rather than duplicate",
      length(prof(gsrc, regions = c("chrA", "chrA:1-100"), verbose = FALSE)) == 2L)
check("bare integers address ranges by position",
      length(prof(gsrc, regions = c(1, 3), verbose = FALSE)) == 2L)
check("an unknown seqname errors",
      err(.tm_regions("chrZ", .tm_source(gsrc))))
check("selected ranges keep their whole sequence, not a clipped piece",
      isTRUE(all.equal(prof(gsrc, regions = "chrA:1-10", verbose = FALSE)$Tm,
                       tm_calculate(unname(oligos[1]))$gr$Tm)))
check("a GRanges without a sequence column is refused",
      err(.tm_source(GRanges("chrA", IRanges(1, 60)))))

## -- 5. Parallel must not change the answer ---------------------------------
cat("\n5. Parallel execution\n")
if (!requireNamespace("BiocParallel", quietly = TRUE)) {
  cat("  SKIP  BiocParallel not installed\n")
} else {
  suppressPackageStartupMessages(library(BiocParallel))
  par <- safe(prof(PKG, regions = REG, window = 200, slide = 200,
                   unit = "segment", segment_size = 200e3,
                   BPPARAM = SnowParam(workers = 2), verbose = FALSE))
  check("two workers reproduce the serial result", same(par, one))
  pseq <- safe(prof(named, BPPARAM = SnowParam(workers = 2), verbose = FALSE))
  check("staged sequences survive two workers",
        isTRUE(all.equal(sort(pseq$Tm), sort(ref$Tm))))
}

## -- 6. The deprecated alias ------------------------------------------------
cat("\n6. tm_profile() still works, once, with a warning\n")
check("warns", warns(tm_profile(PKG, regions = REG, window = 200,
                                slide = 200, verbose = FALSE)))
check("returns a bare GRanges",
      is(suppressWarnings(tm_profile(PKG, regions = REG, window = 200,
                                     slide = 200, verbose = FALSE)), "GRanges"))

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0L) quit(status = 1L)
