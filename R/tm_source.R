# ===========================================================================
# Where sequence comes from, and which part of it to take.
#
# tm_calculate() accepts four kinds of input and one `regions` argument that
# means the same thing in all four. The machinery for that lives here so the
# user-facing function stays readable.
#
#   source                 regions selects                 parallel by
#   ---------------------  ------------------------------  ---------------
#   BSgenome package name  chromosomes, by name or number  region or segment
#   FASTA path             records, by name or number      region or segment
#   character vector       elements, by name or position   staged to FASTA
#   GRanges with sequence  ranges, by overlap or seqname   staged to FASTA
#
# The rule that makes `regions` uniform: the identifier before the colon is
# resolved against whatever names the source itself offers, and falls back to
# position when the source is unnamed or the identifier is a bare integer
# that is not one of those names. So "chr1" is a chromosome in a BSgenome, a
# record in a FASTA and a seqname in a GRanges, and "1:1-200" is the first
# sequence's first 200 bases in an unnamed character vector.
# ===========================================================================

#' Classify the input and describe what it offers
#'
#' @param x The user's \code{input_seq}.
#' @param complement_seq Passed through to \code{\link{to_genomic_ranges}}
#'   when \code{x} is a character vector or a file.
#' @return A list with \code{kind}, the reference needed to reopen the source
#'   (\code{pkg}, \code{path} or \code{gr}), and \code{lens}, the length of
#'   every chromosome, record or range it contains.
#' @keywords internal
.tm_source <- function(x, complement_seq = NULL) {
  if (methods::is(x, "GRanges")) {
    if (is.null(GenomicRanges::mcols(x)$sequence))
      stop("A GRanges given as 'input_seq' must carry a 'sequence' column. ",
           "To profile coordinates against a genome, pass the genome as ",
           "'input_seq' and the coordinates as 'regions'.")
    lens <- stats::setNames(nchar(as.character(GenomicRanges::mcols(x)$sequence)),
                            as.character(GenomeInfoDb::seqnames(x)))
    return(list(kind = "granges", gr = x, lens = lens))
  }
  if (!is.character(x) || !length(x))
    stop("'input_seq' must be a character vector, a FASTA path, the name of ",
         "an installed BSgenome package, or a GRanges carrying sequences.")

  if (length(x) == 1L && file.exists(x) && !dir.exists(x))
    return(list(kind = "fasta", path = x, lens = .fasta_lengths(x)))

  if (length(x) == 1L && requireNamespace(x, quietly = TRUE) &&
      exists(sub("^BSgenome\\.([^.]+)\\..*$", "\\1", x), envir = asNamespace(x))) {
    g <- get(sub("^BSgenome\\.([^.]+)\\..*$", "\\1", x), envir = asNamespace(x))
    return(list(kind = "bsgenome", pkg = x, obj = g,
                lens = GenomeInfoDb::seqlengths(g)))
  }
  if (length(x) == 1L && requireNamespace(x, quietly = TRUE)) {
    # A BSgenome package whose object is not named after the organism.
    nm <- ls(asNamespace(x))
    hit <- nm[vapply(nm, function(n)
      methods::is(get(n, envir = asNamespace(x)), "BSgenome"), logical(1))]
    if (length(hit)) {
      g <- get(hit[1], envir = asNamespace(x))
      return(list(kind = "bsgenome", pkg = x, obj = g,
                  lens = GenomeInfoDb::seqlengths(g)))
    }
  }
  # Anything else that is character is a vector of sequences. A single string
  # that is neither a file nor a package lands here too, which is right: one
  # sequence is a legitimate input.
  lens <- stats::setNames(nchar(x), if (is.null(names(x))) rep("", length(x)) else names(x))
  list(kind = "sequences", seqs = x, complement = complement_seq, lens = lens)
}


#' Where a source's records sit in the coordinate system the caller cares about
#'
#' Staging a \code{GRanges} to a FASTA turns each range into a record that
#' starts at 1, so without this the profile of a range at chr1:1001-1060
#' would come back as chr1:1-60. Every other source already numbers from the
#' coordinate the caller asked about, so its offset is zero.
#'
#' @param idx Row indices into the source.
#' @param src A source from \code{\link{.tm_source}}.
#' @return A numeric vector of offsets to add to record-relative positions.
#' @keywords internal
.tm_offsets <- function(idx, src) {
  if (src$kind != "granges") return(rep(0, length(idx)))
  as.numeric(BiocGenerics::start(src$gr))[idx] - 1
}


#' Resolve `regions` against a source
#'
#' Accepts names, numbers, \code{"name:start-end"} strings, a mixture, or a
#' \code{GRanges}. Returns one row per requested region, in source order.
#' An unresolvable name or an out-of-bounds coordinate is an error rather
#' than a silent drop: quietly skipping a region would return a profile that
#' is short by that region with nothing to say so.
#'
#' @param regions The user's \code{regions} argument, or \code{NULL} for all.
#' @param src A source from \code{\link{.tm_source}}.
#' @return A data frame with \code{idx}, \code{name}, \code{start},
#'   \code{end}, \code{whole}.
#' @keywords internal
.tm_regions <- function(regions, src) {
  lens  <- src$lens
  named <- nzchar(names(lens))
  # A BSgenome with no regions asked for means the standard chromosomes;
  # every other source means everything it contains.
  if (is.null(regions)) {
    idx <- if (src$kind == "bsgenome") {
      std <- GenomeInfoDb::standardChromosomes(src$obj)
      match(intersect(std, names(lens)), names(lens))
    } else seq_along(lens)
    return(data.frame(idx = idx, name = names(lens)[idx],
                      start = 1, end = as.numeric(lens[idx]), whole = TRUE,
                      offset = .tm_offsets(idx, src),
                      stringsAsFactors = FALSE))
  }

  # A GRanges source already carries coordinates, so `regions` is a question
  # about overlap rather than about extraction: which of these ranges does
  # the caller want. Whole ranges come back, not clipped pieces of them, on
  # the reasoning that the sequences are already fixed and a clipped range
  # would need its sequence recut to stay honest.
  if (src$kind == "granges") {
    query <- if (methods::is(regions, "GRanges")) regions
             else .tm_as_granges(regions, src)
    if (is.null(query)) {                      # positional: the nth ranges
      idx <- .tm_match(as.character(regions), names(lens), src)
    } else {
      idx <- sort(unique(S4Vectors::queryHits(
        suppressWarnings(GenomicRanges::findOverlaps(src$gr, query)))))
      if (!length(idx)) stop("No range of 'input_seq' overlaps 'regions'.")
    }
    return(data.frame(idx = idx, name = names(lens)[idx], start = 1,
                      end = as.numeric(lens[idx]), whole = TRUE,
                      offset = .tm_offsets(idx, src),
                      stringsAsFactors = FALSE))
  }

  if (methods::is(regions, "GRanges")) {
    chr <- as.character(GenomeInfoDb::seqnames(regions))
    st  <- as.numeric(BiocGenerics::start(regions))
    en  <- as.numeric(BiocGenerics::end(regions))
  } else {
    r <- as.character(regions)
    if (!length(r)) stop("'regions' is empty.")
    if (anyNA(r))   stop("'regions' contains NA.")
    has <- grepl(":", r, fixed = TRUE)
    chr <- ifelse(has, sub(":.*$", "", r), r)
    st  <- en <- rep(NA_real_, length(r))
    if (any(has)) {
      co <- gsub(",", "", sub("^[^:]*:", "", r[has]), fixed = TRUE)
      pr <- strsplit(co, "-", fixed = TRUE)
      bad <- vapply(pr, length, integer(1)) != 2L
      if (any(bad))
        stop("Region(s) not in \"name:start-end\" form: ",
             paste(r[has][bad], collapse = ", "))
      a <- suppressWarnings(as.numeric(vapply(pr, `[`, character(1), 1L)))
      b <- suppressWarnings(as.numeric(vapply(pr, `[`, character(1), 2L)))
      if (anyNA(a) || anyNA(b))
        stop("Region(s) with non-numeric coordinates: ",
             paste(r[has][is.na(a) | is.na(b)], collapse = ", "))
      st[has] <- a; en[has] <- b
    }
  }

  idx <- .tm_match(chr, names(lens), src)
  full <- is.na(st) & is.na(en)
  st[full] <- 1; en[full] <- as.numeric(lens[idx[full]])

  if (any(st < 1))
    stop("Region start below 1: ", paste(chr[st < 1], collapse = ", "))
  if (any(st > en))
    stop("Region start after its end: ",
         paste(sprintf("%s:%.0f-%.0f", chr[st > en], st[st > en], en[st > en]),
               collapse = ", "))
  over <- en > as.numeric(lens[idx])
  if (any(over))
    stop("Region end beyond the length of its sequence: ",
         paste(sprintf("%s:%.0f-%.0f (length %.0f)", chr[over], st[over],
                       en[over], as.numeric(lens[idx[over]])), collapse = ", "))

  out <- data.frame(idx = idx, name = names(lens)[idx], start = st, end = en,
                    whole = st == 1 & en == as.numeric(lens[idx]),
                    offset = .tm_offsets(idx, src),
                    stringsAsFactors = FALSE)
  out$name[!nzchar(out$name)] <- sprintf("seq%d", out$idx[!nzchar(out$name)])
  out <- out[order(out$idx, out$start), , drop = FALSE]
  if (nrow(out) > 1L) {
    same <- out$idx[-1] == out$idx[-nrow(out)]
    if (any(same & out$start[-1] <= out$end[-nrow(out)]))
      warning("Requested regions overlap; the overlapping windows will ",
              "appear more than once in the result.", call. = FALSE)
  }
  rownames(out) <- NULL
  out
}


#' Read `regions` as an interval query against a GRanges source
#'
#' \code{"chr1"} means the whole of that \code{seqname} and
#' \code{"chr1:1000-2000"} means that interval of it. A vector of bare
#' integers is positional instead, the nth ranges, which is what the same
#' input means for a FASTA file or for sequences; \code{NULL} is returned in
#' that case so the caller resolves it by position.
#' @param regions Character or numeric regions.
#' @param src A GRanges source from \code{\link{.tm_source}}.
#' @return A \code{GRanges} to overlap against, or \code{NULL} for positional.
#' @keywords internal
.tm_as_granges <- function(regions, src) {
  r <- as.character(regions)
  lv <- unique(as.character(GenomeInfoDb::seqnames(src$gr)))
  # Bare integers that are not seqnames address ranges by position.
  n <- suppressWarnings(as.integer(r))
  if (all(!is.na(n)) && !any(r %in% lv)) return(NULL)

  has <- grepl(":", r, fixed = TRUE)
  chr <- ifelse(has, sub(":.*$", "", r), r)
  st  <- rep(1, length(r))
  en  <- rep(.Machine$integer.max, length(r))
  if (any(has)) {
    co <- gsub(",", "", sub("^[^:]*:", "", r[has]), fixed = TRUE)
    pr <- strsplit(co, "-", fixed = TRUE)
    bad <- vapply(pr, length, integer(1)) != 2L
    if (any(bad))
      stop("Region(s) not in \"name:start-end\" form: ",
           paste(r[has][bad], collapse = ", "))
    a <- suppressWarnings(as.numeric(vapply(pr, `[`, character(1), 1L)))
    b <- suppressWarnings(as.numeric(vapply(pr, `[`, character(1), 2L)))
    if (anyNA(a) || anyNA(b))
      stop("Region(s) with non-numeric coordinates: ",
           paste(r[has][is.na(a) | is.na(b)], collapse = ", "))
    st[has] <- a; en[has] <- b
  }
  miss <- setdiff(chr, lv)
  if (length(miss))
    stop("Not found among the seqnames of 'input_seq': ",
         paste(unique(miss), collapse = ", "),
         "\nAvailable: ", paste(utils::head(lv, 8), collapse = ", "),
         if (length(lv) > 8) " ..." else "")
  GenomicRanges::GRanges(chr, IRanges::IRanges(as.integer(st), as.integer(en)))
}


#' Resolve one identifier against a source's own names
#'
#' Names first, then position. On a BSgenome the \code{chr} prefix is added
#' or removed as the genome requires, since that is a naming convention
#' rather than information; FASTA records and \code{seqnames} are matched
#' exactly, because theirs are arbitrary and guessing could match the wrong
#' one.
#' @param x Identifiers to resolve.
#' @param available The source's names, possibly empty strings.
#' @param src The source, used for its kind and for error messages.
#' @return Integer positions into \code{available}.
#' @keywords internal
.tm_match <- function(x, available, src) {
  x <- as.character(x)
  hit <- match(x, available)
  if (src$kind == "bsgenome") {
    miss <- is.na(hit)
    hit[miss] <- match(paste0("chr", x[miss]), available)
    miss <- is.na(hit)
    hit[miss] <- match(sub("^chr", "", x[miss]), available)
  }
  # A bare integer that is not a name means the nth sequence.
  miss <- is.na(hit)
  if (any(miss)) {
    n <- suppressWarnings(as.integer(x[miss]))
    ok <- !is.na(n) & n >= 1 & n <= length(available)
    hit[miss][ok] <- n[ok]
  }
  if (anyNA(hit)) {
    lab <- switch(src$kind, bsgenome = src$pkg, fasta = src$path,
                  granges = "the GRanges", "the sequences")
    shown <- available[nzchar(available)]
    stop("Not found in ", lab, ": ",
         paste(unique(x[is.na(hit)]), collapse = ", "),
         if (length(shown))
           paste0("\nAvailable: ", paste(utils::head(shown, 8), collapse = ", "),
                  if (length(shown) > 8) " ..." else "")
         else sprintf("\nThe source has %d unnamed sequences; use a number.",
                      length(available)))
  }
  hit
}


#' Cut the requested regions into tasks
#'
#' One task is what a single worker owns from start to finish. Long regions
#' are split into \code{segment_size} pieces; short whole records are merged,
#' because reading record n of a file means skipping the first n - 1, so a
#' task per record would rescan the file once per record.
#' @param req Regions from \code{\link{.tm_regions}}.
#' @param src A source from \code{\link{.tm_source}}.
#' @param unit,segment_size,step Task granularity.
#' @return A list of tasks.
#' @keywords internal
.tm_tasks <- function(req, src, unit, segment_size, step) {
  pieces <- list()
  for (i in seq_len(nrow(req))) {
    base <- list(idx = req$idx[i], name = req$name[i],
                 offset = if (is.null(req$offset)) 0 else req$offset[i],
                 whole = req$whole[i] && src$kind == "bsgenome")
    if (unit == "region" || (req$end[i] - req$start[i] + 1) <= segment_size) {
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
  if (src$kind == "bsgenome") return(pieces)

  # Record-based sources: merge consecutive whole records into one read.
  tasks <- list(); run <- NULL
  flush <- function() if (!is.null(run)) tasks[[length(tasks) + 1L]] <<- run
  for (p in pieces) {
    joinable <- isTRUE(p$entire) && unit == "segment" && !is.null(run) &&
      isTRUE(run$mergeable) && p$idx == run$idx2 + 1L &&
      (run$bases + (p$end - p$start + 1)) <= segment_size
    if (joinable) {
      run$idx2   <- p$idx
      run$names  <- c(run$names, p$name)
      run$offs   <- c(run$offs, p$offset)
      run$starts <- c(run$starts, p$start)
      run$ends   <- c(run$ends, p$end)
      run$bases  <- run$bases + (p$end - p$start + 1)
    } else {
      flush()
      run <- list(idx1 = p$idx, idx2 = p$idx, names = p$name,
                  offs = p$offset, starts = p$start, ends = p$end,
                  bases = p$end - p$start + 1, mergeable = isTRUE(p$entire))
    }
  }
  flush()
  tasks
}


#' Record lengths of a FASTA file, without reading the sequences
#' @param path Path to a FASTA file, optionally gzipped.
#' @return A named numeric vector of record lengths.
#' @keywords internal
.fasta_lengths <- function(path) {
  lens <- Biostrings::fasta.seqlengths(path)
  if (!length(lens)) stop("No sequences found in the FASTA file: ", path)
  names(lens) <- sub("\\s.*$", "", names(lens))
  if (anyDuplicated(names(lens)))
    stop("Duplicated record names in ", path, ": ",
         paste(unique(names(lens)[duplicated(names(lens))]), collapse = ", "))
  lens
}


#' Stage sequences as a FASTA file so that workers read rather than receive them
#'
#' Written in blocks: an XStringSet over the whole vector would double peak
#' memory at exactly the input size where this path is worth using. The
#' record name carries the input position, so a sequence dropped later for
#' containing \code{N} can still be identified.
#' @param seqs Character vector of sequences.
#' @param tmpdir Directory for the file.
#' @return The path, carrying the input count and names as attributes.
#' @keywords internal
.spill_fasta <- function(seqs, tmpdir = tempdir()) {
  if (!dir.exists(tmpdir)) stop("'tmpdir' does not exist: ", tmpdir)
  path <- tempfile(pattern = "tm_", tmpdir = tmpdir, fileext = ".fa")
  n <- length(seqs); nm <- names(seqs); block <- 1e5L
  # Record names are the caller's own when they are usable, so that
  # regions = "myseq" resolves the way it does for a FASTA or a genome.
  # Otherwise they are positions, which is also what regions = 3 means.
  usable <- !is.null(nm) && all(nzchar(nm)) && !anyDuplicated(nm)
  for (from in seq(1L, n, by = block)) {
    ii <- from:min(from + block - 1L, n)
    x  <- Biostrings::BStringSet(unname(seqs[ii]))
    names(x) <- if (usable) nm[ii] else sprintf("seq%d", ii)
    Biostrings::writeXStringSet(x, path, append = from > 1L)
  }
  attr(path, "n") <- n; attr(path, "names_in") <- if (usable) NULL else nm
  path
}


# ---------------------------------------------------------------- the runner

#' Run the tasks and reassemble one profile
#'
#' Each task opens the source itself, so what crosses between processes is a
#' name and a coordinate pair rather than any sequence. That is what makes
#' parallelism pay here, and why an in-memory source is staged to a file
#' first rather than sent over a socket.
#' @param tasks From \code{\link{.tm_tasks}}.
#' @param src A source from \code{\link{.tm_source}}.
#' @param window,slide Tiling; \code{NULL} window means one window per region.
#' @param model Model arguments to hand to \code{\link{tm_calculate}}.
#' @param BPPARAM A \code{BiocParallelParam}, or \code{NULL} for this process.
#' @param keep_sequence Keep the sequence columns in the result.
#' @param verbose Report task and window counts.
#' @return A \code{GRanges} in source order.
#' @importFrom BiocParallel bplapply bptasks<-
#' @keywords internal
.tm_run <- function(tasks, src, window, slide, model, BPPARAM,
                    keep_sequence, verbose) {
  if (!length(tasks)) stop("No regions to tile.")
  bases <- vapply(tasks, function(z)
    if (is.null(z$bases)) z$end - z$start + 1 else z$bases, numeric(1))
  ord <- order(bases, decreasing = TRUE)      # longest first: fills the tail
  if (verbose)
    message(sprintf("tm_calculate: %d task%s over %s",
                    length(tasks), if (length(tasks) == 1L) "" else "s",
                    switch(src$kind, bsgenome = "the genome",
                           fasta = "the FASTA file", "the sequences")))

  worker <- if (src$kind == "bsgenome") .tm_task_bsgenome else .tm_task_fasta
  ref    <- if (src$kind == "bsgenome") src$pkg else src$path
  args   <- list(src = ref, window = window, slide = slide,
                 keep_sequence = keep_sequence, model = model)

  if (is.null(BPPARAM)) {
    res <- lapply(tasks[ord], function(tk) do.call(worker, c(list(tk), args)))
  } else {
    # One task at a time rather than one pre-split block per worker: with
    # tasks of unequal length, pre-splitting leaves workers idle at the end.
    BiocParallel::bptasks(BPPARAM) <- length(tasks)
    res <- BiocParallel::bplapply(tasks[ord], worker, BPPARAM = BPPARAM,
                                  src = ref, window = window, slide = slide,
                                  keep_sequence = keep_sequence, model = model)
  }
  res  <- res[order(ord)]
  keep <- vapply(res, length, integer(1)) > 0L
  if (!any(keep))
    stop("No windows were produced. Check 'window', 'slide' and the region widths.")
  out <- unlist(GenomicRanges::GRangesList(res[keep]), use.names = FALSE)
  if (verbose)
    message(sprintf("tm_calculate: %s windows",
                    format(length(out), big.mark = ",")))
  out
}

#' One task against a BSgenome
#' @param task A task from \code{\link{.tm_tasks}}.
#' @param src BSgenome package name.
#' @param window,slide,keep_sequence,model As in \code{\link{.tm_run}}.
#' @return A \code{GRanges}.
#' @keywords internal
.tm_task_bsgenome <- function(task, src, window, slide, keep_sequence, model) {
  suppressPackageStartupMessages({
    requireNamespace("TmCalculator", quietly = TRUE)
    requireNamespace(src, quietly = TRUE)
  })
  w <- if (is.null(window)) task$end - task$start + 1L else window
  s <- if (is.null(slide))  w                          else slide
  bins <- TmCalculator::make_genomiccoord(
    bsgenome = src, chromosomes = task$name, window = w, slide = s,
    start = task$start, end = task$end, strand = "+",
    # A whole chromosome has its telomeric N runs trimmed; a region the
    # caller named starts where they said it does.
    trim_N = if (isTRUE(task$whole)) "ends" else "none", verbose = FALSE)
  gr <- TmCalculator::to_genomic_ranges_fast(
    list(pkg_name = src, seq = bins), method = "preload_chr")
  .tm_finish(gr, keep_sequence, model)
}

#' One task against a FASTA file
#'
#' Reads only its own records, with \code{skip} and \code{nrec}, so the
#' sequence never travels between processes.
#' @inheritParams .tm_task_bsgenome
#' @return A \code{GRanges}.
#' @keywords internal
.tm_task_fasta <- function(task, src, window, slide, keep_sequence, model) {
  suppressPackageStartupMessages(requireNamespace("TmCalculator", quietly = TRUE))
  n_rec <- task$idx2 - task$idx1 + 1L
  recs  <- Biostrings::readBStringSet(src, skip = task$idx1 - 1L, nrec = n_rec)
  if (length(recs) != n_rec)
    stop("FASTA records ", task$idx1, "-", task$idx2, " unreadable from ", src)

  seqs <- character(0)
  for (k in seq_len(n_rec)) {
    x <- recs[[k]]
    a <- task$starts[k]; b <- task$ends[k]
    if (a > 1L || b < length(x)) x <- Biostrings::subseq(x, start = a, end = b)
    n <- length(x)
    w <- if (is.null(window)) n else min(window, n)  # short records stay whole
    s <- if (is.null(slide))  w else slide
    starts <- if (n <= w) 1L else seq(1L, n - w + 1L, by = s)
    piece  <- as.character(Biostrings::subseq(rep(x, length(starts)),
                                              start = starts, width = w))
    # Record-relative position, plus where the record itself begins in the
    # caller's coordinates. The first term makes segments of one record
    # reassemble; the second keeps a staged GRanges on its original ranges.
    off <- a - 1L + (if (is.null(task$offs)) 0 else task$offs[k])
    names(piece) <- sprintf("%s:%.0f-%.0f", task$names[k],
                            off + starts, off + starts + w - 1L)
    seqs <- c(seqs, piece)
  }
  .tm_finish(TmCalculator::to_genomic_ranges(seqs), keep_sequence, model)
}

#' Compute Tm for one task's windows
#' @param gr Windows with sequences.
#' @param keep_sequence Keep the sequence columns.
#' @param model Model arguments for \code{\link{tm_calculate}}.
#' @return A \code{GRanges}.
#' @keywords internal
.tm_finish <- function(gr, keep_sequence, model) {
  out <- do.call(TmCalculator::tm_calculate, c(list(input_seq = gr), model))$gr
  if (!keep_sequence) { out$sequence <- NULL; out$complement <- NULL }
  out
}
