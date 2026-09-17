#' Melting temperature profile over a genome, a FASTA file, or a set of sequences
#'
#' Builds a Tm profile from a sequence source: a \pkg{BSgenome} package, a
#' FASTA file, or a vector of sequences. The requested regions are tiled
#' into windows and a melting temperature is computed for each. This is the
#' batch and parallel entry point of the package;
#' \code{\link{tm_calculate}} is the single-call one.
#'
#' Parallelism pays here because a task is self-contained: it opens the
#' source itself and builds its own windows, fetches its own sequence and
#' assembles its own result, so only a name and a coordinate pair cross
#' between processes. A vector of sequences is staged as a temporary FASTA
#' file for the same reason, so that the workers read it rather than receive
#' it. Dividing the sequences of a single \code{tm_calculate()} call among
#' workers is the arrangement that does not pay, because it parallelises the
#' inner loop alone and leaves window construction, sequence retrieval and
#' result assembly in the calling process. See
#' \code{vignette("hg38_performance_parallel", package = "TmCalculator")}.
#'
#' @param seq_source Where the sequence comes from. One of:
#'   \itemize{
#'     \item the name of an installed \pkg{BSgenome} package, e.g.
#'       \code{"BSgenome.Hsapiens.UCSC.hg38"}. A package name rather than a
#'       loaded object, because each worker opens the source itself;
#'     \item the path to a FASTA file, optionally gzipped;
#'     \item a character vector of two or more sequences, which is staged as
#'       a temporary FASTA file under \code{tmpdir} and removed on exit. For
#'       a single sequence, and for a few thousand short ones, the staging
#'       and the worker start-up cost more than the calculation: use
#'       \code{\link{tm_calculate}} there.
#'   }
#' @param regions What to tile. Any of:
#'   \itemize{
#'     \item chromosome or record names, or numbers: \code{1:2},
#'       \code{c("chr1", "chr2")}, \code{c("contig_7", "contig_9")};
#'     \item coordinate strings \code{"name:start-end"}, e.g.
#'       \code{c("chr1:1-10000000", "chrX:5,000,000-6,000,000")}. Commas and
#'       scientific notation in the numbers are accepted;
#'     \item a mixture of the two, \code{c("chr21", "chr1:1-10e6")};
#'     \item a \code{GRanges} object, whose ranges are used directly.
#'   }
#'   For a BSgenome source the \code{chr} prefix is added or removed as the
#'   genome requires, so \code{1:2} works on a UCSC genome and
#'   \code{c("chr1", "chr2")} on an Ensembl one; FASTA record names are
#'   matched exactly, since they are arbitrary. Default: every standard
#'   chromosome of the genome, or every record of the FASTA file.
#' @param window Window width in base pairs. \code{NULL} means one window
#'   per region, which is what short records such as array probes or primers
#'   call for, and is the default when \code{seq_source} is a vector of
#'   sequences. Records shorter than \code{window} are returned whole.
#' @param slide Step between window starts, defaulting to \code{window},
#'   which gives a non-overlapping tiling. Ignored when \code{window} is
#'   \code{NULL}.
#' @param unit How the requested regions are turned into tasks.
#'   \code{"segment"} cuts them into pieces of about \code{segment_size} bp;
#'   \code{"region"} makes one task per requested region. Segments are
#'   faster and need less memory per worker, because no worker ever holds a
#'   whole large chromosome.
#' @param segment_size Task size in base pairs when \code{unit = "segment"}.
#'   Rounded down to a multiple of \code{slide}, and offsets are measured
#'   from each region's own start, so the window grid is identical to the one
#'   a single unsegmented task would produce. Default 50 Mb.
#' @param BPPARAM A \code{BiocParallelParam} from \pkg{BiocParallel}, e.g.
#'   \code{SnowParam(workers = 5)}. \code{NULL} (the default) runs the tasks
#'   sequentially in this process with no dependency on \pkg{BiocParallel}.
#'   Tasks are dispatched longest first, one at a time, which is what keeps
#'   the workers busy to the end of the run.
#' @param keep_sequence Keep the \code{sequence} and \code{complement}
#'   columns in the result. \code{FALSE} by default: they are roughly 500 MB
#'   per large chromosome, and returning them from a worker costs more than
#'   the calculation.
#' @param tmpdir Directory for the temporary FASTA file written when
#'   \code{seq_source} is a vector of sequences. Worth setting on a cluster,
#'   where the default \code{tempdir()} is often a small partition.
#' @param verbose Report the task count and the number of windows.
#' @param ... Passed to \code{\link{tm_calculate}}: \code{method},
#'   \code{nn_table}, \code{Na} and the other model arguments.
#'
#' @return A \code{GRanges} object in genomic order, carrying \code{Tm} and
#'   \code{GC} metadata columns, with the chromosome or record name in its
#'   \code{seqnames}.
#'
#'   When \code{seq_source} is a vector of sequences the result is keyed by
#'   input position, \code{seq1}, \code{seq2} and so on, rather than by the
#'   caller's names, so that a sequence dropped for containing \code{N} can
#'   still be identified; the caller's names, if the vector had any, are
#'   returned in a \code{name} column.
#'
#' @section Assembly gaps:
#'   A task covering a whole chromosome of a BSgenome trims the leading and
#'   trailing assembly gaps, since a telomeric run of N carries no windows.
#'   A task covering a user-specified region, or any FASTA record, does not
#'   trim, because the requested start is the requested start. Windows
#'   containing N are dropped in every case.
#'
#' @section Choosing the worker count:
#'   Memory, not cores, is usually the binding constraint: each worker holds
#'   the sequences of its own task. Measured on GRCh38 at 200 bp, peak
#'   resident memory per worker is 4.74 to 4.87 GB when the unit is a whole
#'   chromosome and 1.86 GB when it is a 50 Mb segment. A workable rule is
#'   the smaller of the number of physical cores less one, and the available
#'   memory less 4 GB divided by 2 GB per worker.
#'
#' @examples
#' \dontrun{
#' library(BiocParallel)
#' hg38 <- "BSgenome.Hsapiens.UCSC.hg38"
#'
#' # Whole genome, 50 Mb tasks, five workers
#' tm <- tm_profile(hg38, window = 200, slide = 200,
#'                  method = "tm_nn", nn_table = "DNA_NN_SantaLucia_2004",
#'                  Na = 50, BPPARAM = SnowParam(workers = 5))
#'
#' # Chromosomes 1 and 2 only. These are equivalent.
#' tm_profile(hg38, regions = 1:2, method = "tm_nn")
#' tm_profile(hg38, regions = c("chr1", "chr2"), method = "tm_nn")
#'
#' # Explicit regions, and a mixture of regions and whole chromosomes
#' tm_profile(hg38, regions = c("chr1:1-10000000", "chrX:5,000,000-6,000,000"))
#' tm_profile(hg38, regions = c("chr21", "chr1:1-10e6"))
#'
#' # Regions from an existing object, for example a set of promoters
#' tm_profile(hg38, regions = promoters_gr, window = 50, slide = 25)
#'
#' # A FASTA assembly, tiled the same way
#' tm_profile("contigs.fa.gz", window = 200, slide = 200,
#'            BPPARAM = SnowParam(workers = 4))
#'
#' # Short records, one Tm each: array probes, primers, synthetic oligos
#' probe_tm <- tm_profile("probes.fa", window = NULL, method = "tm_nn", Na = 50)
#'
#' # Sequences already in R. Staged to a temporary FASTA so that the workers
#' # read them; worth it for large sets, not for a handful.
#' tm_profile(oligos, method = "tm_nn", Na = 50,
#'            BPPARAM = SnowParam(workers = 5))
#' }
#'
#' @seealso \code{\link{tm_calculate}}, which this function calls once per
#'   task, for sequences already in hand; \code{\link{make_genomiccoord}} for
#'   the tiling and \code{\link{integrate_granges}} for combining the
#'   resulting profile with other genomic data.
#' @importFrom Biostrings readBStringSet writeXStringSet fasta.seqlengths
#' @importFrom Biostrings BStringSet subseq
#' @export
tm_profile <- function(seq_source,
                       regions       = NULL,
                       window        = 200L,
                       slide         = window,
                       unit          = c("segment", "region"),
                       segment_size  = 50e6,
                       BPPARAM       = NULL,
                       keep_sequence = FALSE,
                       tmpdir        = tempdir(),
                       verbose       = TRUE,
                       ...) {
  unit <- match.arg(unit)

  # -- Sequences given directly ----------------------------------------------
  # A vector of sequences is written to a temporary FASTA file and then read
  # back by the workers. That looks wasteful and is not: it moves the whole
  # pipeline into the worker, window construction, complement generation and
  # result assembly included, where sending the sequences over a socket would
  # parallelise only the inner loop and leave the rest in this process. One
  # disk round trip buys that, and it is the same round trip R's own
  # serialisation would pay to reach the workers anyway. Below a few thousand
  # sequences the round trip and the worker start-up cost more than the
  # calculation; use tm_calculate() there.
  spilled <- NULL
  if (is.character(seq_source) && length(seq_source) > 1L) {
    if (anyNA(seq_source)) stop("'seq_source' contains NA.")
    if (missing(window)) window <- NULL          # one Tm per input sequence
    spilled    <- .spill_fasta(seq_source, tmpdir)
    seq_source <- spilled
    on.exit(unlink(spilled), add = TRUE)
    if (verbose)
      message(sprintf("tm_profile: %s sequences staged in %s",
                      format(attr(spilled, "n"), big.mark = ","), spilled))
  }

  if (!is.null(window)) {
    window <- as.integer(window)
    slide  <- as.integer(if (is.null(slide)) window else slide)
    if (window < 1L || slide < 1L) stop("'window' and 'slide' must be positive.")
  } else {
    slide <- NULL
  }
  if (!is.character(seq_source) || length(seq_source) != 1L)
    stop("'seq_source' must be a single string: an installed BSgenome package ",
         "name, or the path to a FASTA file.")

  # A path that exists is a file; a BSgenome package name never is.
  is_fasta <- file.exists(seq_source) && !dir.exists(seq_source)
  if (is_fasta) {
    lens <- .fasta_lengths(seq_source)
  } else {
    if (!requireNamespace(seq_source, quietly = TRUE))
      stop("'seq_source' is neither a file nor an installed BSgenome package: ",
           seq_source, "\n  A single sequence belongs in tm_calculate(); ",
           "two or more are accepted here as a vector.")
    genome <- get(seq_source, envir = asNamespace(seq_source))
    lens   <- GenomeInfoDb::seqlengths(genome)
    if (is.null(regions)) {
      regions <- GenomeInfoDb::standardChromosomes(genome)
      if (length(regions) == 0L) regions <- names(lens)
    }
  }
  if (is.null(regions)) regions <- names(lens)
  req <- .parse_regions(regions, lens, seq_source, fuzzy = !is_fasta)

  # -- Tasks ------------------------------------------------------------------
  # A task is one region, or one segment of a long region. For FASTA the
  # short regions are then merged: reading record n means skipping the first
  # n - 1, so a task per record would rescan the file once per record, which
  # is quadratic on a file of many short records such as a probe set.
  step  <- if (is.null(slide)) 1L else slide
  req$idx <- match(req$chr, names(lens))
  pieces <- list()
  for (i in seq_len(nrow(req))) {
    base <- list(name = req$chr[i], idx = req$idx[i],
                 whole = req$whole[i] && !is_fasta)
    span <- req$end[i] - req$start[i] + 1
    if (unit == "region" || span <= segment_size) {
      pieces[[length(pieces) + 1L]] <-
        c(base, list(start = req$start[i], end = req$end[i], entire = req$whole[i]))
    } else {
      seg <- max(step, floor(segment_size / step) * step)
      for (s in seq(req$start[i], req$end[i], by = seg))
        pieces[[length(pieces) + 1L]] <-
          c(base, list(start = s, end = min(s + seg - 1L, req$end[i]),
                       whole = FALSE, entire = FALSE))
    }
  }

  if (!is_fasta) {
    tasks <- pieces
  } else {
    tasks <- list(); run <- NULL
    flush <- function() if (!is.null(run)) tasks[[length(tasks) + 1L]] <<- run
    for (p in pieces) {
      # Both sides must be entire records. readBStringSet() reads whole
      # records, so appending one to a task that only wants a slice of an
      # earlier record would pull that whole record into the worker as well.
      mergeable <- isTRUE(p$entire) && unit == "segment" &&
        !is.null(run) && isTRUE(run$mergeable) && p$idx == run$idx2 + 1L &&
        (run$bases + (p$end - p$start + 1)) <= segment_size
      if (mergeable) {
        run$idx2   <- p$idx
        run$names  <- c(run$names, p$name)
        run$starts <- c(run$starts, p$start)
        run$ends   <- c(run$ends, p$end)
        run$bases  <- run$bases + (p$end - p$start + 1)
      } else {
        flush()
        run <- list(idx1 = p$idx, idx2 = p$idx, names = p$name,
                    starts = p$start, ends = p$end,
                    bases = p$end - p$start + 1,
                    mergeable = isTRUE(p$entire))
      }
    }
    flush()
  }
  if (!length(tasks)) stop("No regions to tile.")
  bases <- vapply(tasks, function(z)
    if (is.null(z$bases)) z$end - z$start + 1 else z$bases, numeric(1))
  ord <- order(bases, decreasing = TRUE)

  if (verbose)
    message(sprintf("tm_profile: %d %s task%s over %d %s",
                    length(tasks), unit, if (length(tasks) == 1L) "" else "s",
                    length(unique(req$chr)),
                    if (is_fasta) "FASTA record(s)" else "chromosome(s)"))

  # -- Dispatch ---------------------------------------------------------------
  worker <- if (is_fasta) .tm_task_fasta else .tm_task_bsgenome
  args <- list(src = seq_source, window = window, slide = slide,
               keep_sequence = keep_sequence, dots = list(...))
  if (is.null(BPPARAM)) {
    res <- lapply(tasks[ord], function(tk) do.call(worker, c(list(tk), args)))
  } else {
    if (!requireNamespace("BiocParallel", quietly = TRUE))
      stop("BPPARAM was supplied but BiocParallel is not installed. ",
           "Install it, or leave BPPARAM = NULL to run the tasks serially.")
    # One task at a time rather than one pre-split block per worker: with
    # tasks of unequal length, pre-splitting leaves workers idle at the end.
    BiocParallel::bptasks(BPPARAM) <- length(tasks)
    res <- BiocParallel::bplapply(tasks[ord], worker, BPPARAM = BPPARAM,
                                  src = seq_source, window = window, slide = slide,
                                  keep_sequence = keep_sequence, dots = list(...))
  }

  res  <- res[order(ord)]
  keep <- vapply(res, length, integer(1)) > 0L
  if (!any(keep))
    stop("No windows were produced. Check 'window', 'slide' and the region widths.")
  out <- unlist(GenomicRanges::GRangesList(res[keep]), use.names = FALSE)

  if (!is.null(spilled)) {
    # Sequences are keyed by input position, not by their own name, so that a
    # sequence dropped for containing N can still be identified. The caller's
    # names, if any, come back as a column.
    idx <- as.integer(sub("^seq", "", as.character(GenomeInfoDb::seqnames(out))))
    out <- out[order(idx)]
    nm  <- attr(spilled, "names_in")
    if (!is.null(nm))
      GenomicRanges::mcols(out)$name <- nm[sort(idx)]
  }
  if (verbose)
    message(sprintf("tm_profile: %s windows", format(length(out), big.mark = ",")))
  out
}


# -- Stage a vector of sequences as a FASTA file -----------------------------
# Written in blocks so that the file, rather than a second copy of the
# sequences, is what holds the data: a BStringSet over the whole vector would
# double peak memory at exactly the input size where this path is worth using.
#' @keywords internal
.spill_fasta <- function(seqs, tmpdir = tempdir()) {
  if (!dir.exists(tmpdir))
    stop("'tmpdir' does not exist: ", tmpdir)
  path <- tempfile(pattern = "tm_profile_", tmpdir = tmpdir, fileext = ".fa")
  n    <- length(seqs)
  nm   <- names(seqs)
  # The record name carries the input position. Anything after the first word
  # is a comment as far as this package is concerned, so the caller's own name
  # rides along without being used for matching.
  block <- 1e5L
  for (from in seq(1L, n, by = block)) {
    to  <- min(from + block - 1L, n)
    ii  <- from:to
    x   <- Biostrings::BStringSet(unname(seqs[ii]))
    names(x) <- if (is.null(nm)) sprintf("seq%d", ii)
                else sprintf("seq%d %s", ii, nm[ii])
    Biostrings::writeXStringSet(x, path, append = from > 1L)
  }
  attr(path, "n")        <- n
  attr(path, "names_in") <- nm
  path
}


# -- One task against a BSgenome ---------------------------------------------
#' @keywords internal
.tm_task_bsgenome <- function(task, src, window, slide, keep_sequence, dots) {
  suppressPackageStartupMessages({
    requireNamespace("TmCalculator", quietly = TRUE)
    requireNamespace(src, quietly = TRUE)
  })
  w <- if (is.null(window)) task$end - task$start + 1L else window
  s <- if (is.null(slide))  w                          else slide
  bins <- TmCalculator::make_genomiccoord(
    bsgenome = src, chromosomes = task$name,
    window = w, slide = s,
    start = task$start, end = task$end, strand = "+",
    trim_N = if (isTRUE(task$whole)) "ends" else "none",
    verbose = FALSE)
  gr <- TmCalculator::to_genomic_ranges_fast(
    list(pkg_name = src, seq = bins), method = "preload_chr")
  .tm_task_finish(gr, keep_sequence, dots)
}

# -- One task against a FASTA file -------------------------------------------
# The worker reads only its own records, using skip/nrec, so the sequence
# never travels between processes; what is sent is a record index.
#' @keywords internal
.tm_task_fasta <- function(task, src, window, slide, keep_sequence, dots) {
  suppressPackageStartupMessages(requireNamespace("TmCalculator", quietly = TRUE))
  n_rec <- task$idx2 - task$idx1 + 1L
  recs  <- Biostrings::readBStringSet(src, skip = task$idx1 - 1L, nrec = n_rec)
  if (length(recs) != n_rec)
    stop("FASTA records ", task$idx1, "-", task$idx2, " could not be read from ", src)

  seqs <- character(0)
  for (k in seq_len(n_rec)) {
    x <- recs[[k]]
    a <- task$starts[k]; b <- task$ends[k]
    if (a > 1L || b < length(x)) x <- Biostrings::subseq(x, start = a, end = b)
    n <- length(x)
    w <- if (is.null(window)) n else min(window, n)  # a short record stays whole
    s <- if (is.null(slide))  w else slide
    starts <- if (n <= w) 1L else seq(1L, n - w + 1L, by = s)
    piece  <- as.character(Biostrings::subseq(rep(x, length(starts)),
                                              start = starts, width = w))
    # Absolute coordinates, so a segmented record reassembles correctly.
    off <- a - 1L
    names(piece) <- sprintf("%s:%d-%d", task$names[k],
                            off + starts, off + starts + w - 1L)
    seqs <- c(seqs, piece)
  }
  gr <- TmCalculator::to_genomic_ranges(seqs)
  .tm_task_finish(gr, keep_sequence, dots)
}

#' @keywords internal
.tm_task_finish <- function(gr, keep_sequence, dots) {
  out <- do.call(TmCalculator::tm_calculate, c(list(input_seq = gr), dots))$gr
  if (!keep_sequence) {
    out$sequence   <- NULL
    out$complement <- NULL
  }
  out
}

# -- Record lengths of a FASTA file, without reading the sequences -----------
#' @keywords internal
.fasta_lengths <- function(path) {
  lens <- Biostrings::fasta.seqlengths(path)
  if (!length(lens)) stop("No sequences found in the FASTA file: ", path)
  names(lens) <- sub("\\s.*$", "", names(lens))   # first word, as elsewhere
  if (anyDuplicated(names(lens)))
    stop("Duplicated record names in ", path, ": ",
         paste(unique(names(lens)[duplicated(names(lens))]), collapse = ", "))
  lens
}


#' Resolve the regions argument of tm_profile()
#'
#' Accepts chromosome or record names or numbers, \code{"name:start-end"}
#' strings, a mixture of the two, or a \code{GRanges}, and returns a data
#' frame of resolved regions in source order with a flag marking those that
#' cover a whole chromosome or record. An unresolvable name or an
#' out-of-bounds coordinate is an error rather than a silent drop: quietly
#' skipping a region would return a profile that is short by that region
#' with nothing to say so.
#' @param x The user's \code{regions} argument.
#' @param lens Named vector of chromosome or record lengths.
#' @param label Source name, used in error messages.
#' @param fuzzy Allow the \code{chr} prefix to be added or removed when
#'   matching. True for a BSgenome, whose naming style is a convention;
#'   false for a FASTA file, whose record names are arbitrary.
#' @return A data frame with columns \code{chr}, \code{start}, \code{end},
#'   \code{whole}.
#' @keywords internal
.parse_regions <- function(x, lens, label = "the source", fuzzy = TRUE) {
  if (methods::is(x, "GRanges")) {
    chr <- as.character(GenomeInfoDb::seqnames(x))
    st  <- as.numeric(BiocGenerics::start(x))
    en  <- as.numeric(BiocGenerics::end(x))
  } else {
    x <- as.character(x)
    if (!length(x)) stop("'regions' is empty.")
    if (anyNA(x))   stop("'regions' contains NA.")

    # A bare integer means the nth record when it is not itself a name.
    idx <- suppressWarnings(as.integer(x))
    bare <- !is.na(idx) & !(x %in% names(lens)) & idx >= 1 & idx <= length(lens)
    if (!fuzzy && any(bare)) x[bare] <- names(lens)[idx[bare]]

    has_coord <- grepl(":", x, fixed = TRUE)
    chr <- ifelse(has_coord, sub(":.*$", "", x), x)
    st  <- en <- rep(NA_real_, length(x))

    if (any(has_coord)) {
      # Commas are stripped so that "chr1:5,000,000-6,000,000" parses, and
      # the numbers go through as.numeric so that 10e6 does too.
      coord <- gsub(",", "", sub("^[^:]*:", "", x[has_coord]), fixed = TRUE)
      parts <- strsplit(coord, "-", fixed = TRUE)
      bad   <- vapply(parts, length, integer(1)) != 2L
      if (any(bad))
        stop("Region(s) not in \"name:start-end\" form: ",
             paste(x[has_coord][bad], collapse = ", "))
      s <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 1L)))
      e <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2L)))
      if (anyNA(s) || anyNA(e))
        stop("Region(s) with non-numeric coordinates: ",
             paste(x[has_coord][is.na(s) | is.na(e)], collapse = ", "))
      st[has_coord] <- s
      en[has_coord] <- e
    }
  }

  chr <- .match_seqlevels(chr, names(lens), label, fuzzy = fuzzy)

  full <- is.na(st) & is.na(en)          # a bare name means the whole record
  st[full] <- 1
  en[full] <- as.numeric(lens[chr[full]])

  if (any(st < 1)) stop("Region start below 1: ", paste(chr[st < 1], collapse = ", "))
  if (any(st > en))
    stop("Region start after its end: ",
         paste(sprintf("%s:%.0f-%.0f", chr[st > en], st[st > en], en[st > en]),
               collapse = ", "))
  over <- en > as.numeric(lens[chr])
  if (any(over))
    stop("Region end beyond the length of its chromosome or record: ",
         paste(sprintf("%s:%.0f-%.0f (length %.0f)", chr[over], st[over], en[over],
                       as.numeric(lens[chr[over]])), collapse = ", "))

  out <- data.frame(chr = chr, start = st, end = en,
                    whole = st == 1 & en == as.numeric(lens[chr]),
                    stringsAsFactors = FALSE)
  out <- out[order(match(out$chr, names(lens)), out$start), , drop = FALSE]

  # Overlapping requests would count the same windows twice, and a profile
  # with silently duplicated windows is worse than one that complains.
  if (nrow(out) > 1L) {
    same <- out$chr[-1] == out$chr[-nrow(out)]
    if (any(same & out$start[-1] <= out$end[-nrow(out)]))
      warning("Requested regions overlap; the overlapping windows will appear ",
              "more than once in the result.", call. = FALSE)
  }
  rownames(out) <- NULL
  out
}


#' Resolve names against the records a source offers
#'
#' Accepts what a user is likely to type. With \code{fuzzy = TRUE},
#' \code{1:2} and \code{c("1", "2")} resolve on a UCSC genome and
#' \code{c("chr1", "chr2")} on an Ensembl one, element by element. With
#' \code{fuzzy = FALSE} names must match exactly, which is right for FASTA
#' records: their names are arbitrary, and guessing at a prefix could match
#' the wrong record. Element order is preserved; the caller decides how to
#' sort.
#' @param x Names to resolve.
#' @param available Names present in the source.
#' @param label Source name, used in error messages.
#' @param fuzzy Allow the \code{chr} prefix to be added or removed.
#' @return A character vector of resolved names, the same length as \code{x}.
#' @keywords internal
.match_seqlevels <- function(x, available, label = "the source", fuzzy = TRUE) {
  x <- as.character(x)
  hit <- if (!fuzzy) ifelse(x %in% available, x, NA_character_) else ifelse(
    x %in% available, x,
    ifelse(paste0("chr", x) %in% available, paste0("chr", x),
           ifelse(sub("^chr", "", x) %in% available, sub("^chr", "", x),
                  NA_character_)))
  if (anyNA(hit))
    stop("Not found in ", label, ": ",
         paste(unique(x[is.na(hit)]), collapse = ", "),
         "\nAvailable: ", paste(utils::head(available, 8), collapse = ", "),
         if (length(available) > 8) " ..." else "")
  hit
}
