#' Generate sliding-window genomic coordinate strings for Tm calculation
#'
#' Tiles one or more chromosomes from a \pkg{BSgenome} object into
#' overlapping windows and returns a named character vector of coordinate
#' strings in the format
#' \code{"chr:start-end:strand:genome_pkg:region_id"}.
#' This vector is the primary input for downstream Tm calculation
#' functions (\code{tm_nn}, \code{tm_gc}, \code{tm_wallace}) that
#' accept genomic coordinate strings.
#'
#' @section Coordinate string format:
#' Each element of the returned vector follows the pattern:
#' \preformatted{
#'   chr1:10001-10200:+:BSgenome.Hsapiens.UCSC.hg38:region1
#'   chr1:10051-10250:+:BSgenome.Hsapiens.UCSC.hg38:region2
#'   ...
#' }
#' Fields (colon-separated):
#' \enumerate{
#'   \item \strong{chromosome}  --  e.g. \code{chr1}
#'   \item \strong{start-end}   --  1-based, inclusive coordinates
#'   \item \strong{strand}      --  \code{+} or \code{-}
#'   \item \strong{genome}      --  BSgenome package name (character)
#'   \item \strong{region_id}   --  unique label \code{regionN}
#' }
#'
#' @section N-base trimming:
#' Human (and most mammalian) chromosomes begin and end with long
#' stretches of \code{N} bases that represent assembly gaps.  Windows
#' that consist entirely or predominantly of \code{N}s produce
#' meaningless Tm values.  The function offers two trimming strategies
#' controlled by \code{trim_N}:
#' \describe{
#'   \item{\code{"none"}}{No trimming.  Windows start at \code{start}
#'     and end at \code{end} (or chromosome length).}
#'   \item{\code{"ends"}}{Detect the first and last positions of
#'     non-\code{N} bases on each chromosome using
#'     \code{Biostrings::letterFrequency()} on a coarse block scan,
#'     then trim \code{start} and \code{end} to those positions.
#'     Efficient: reads the chromosome sequence only once in blocks.}
#'   \item{\code{"filter"}}{Generate windows across the full
#'     \code{start}-\code{end} range but then remove any window whose
#'     N-base fraction exceeds \code{max_N_frac} (default 0.1, i.e.
#'     10\%).  More granular than \code{"ends"} but slower because it
#'     reads each window sequence.}
#' }
#'
#' @param bsgenome    A \code{BSgenome} object (e.g.
#'   \code{BSgenome.Hsapiens.UCSC.hg38}) \strong{or} a character
#'   string with the BSgenome package name (loaded automatically).
#' @param chromosomes Character vector of chromosome names to tile.
#'   Must be present in \code{seqlevels(bsgenome)}.  Default:
#'   the 24 standard human chromosomes
#'   \code{paste0("chr", c(1:22, "X", "Y"))}.
#' @param window      Integer.  Width of each sliding window in base
#'   pairs.  Default \code{200L}.
#' @param slide       Integer.  Step size between consecutive window
#'   starts in base pairs.  \code{slide == window} gives non-overlapping
#'   tiling; \code{slide < window} gives overlapping windows.
#'   Default \code{50L}.
#' @param start       Integer \strong{or} \code{NULL}.  Override the
#'   start position for every chromosome.  When \code{NULL} (default)
#'   the start is chromosome position 1 (adjusted for N-trimming if
#'   \code{trim_N != "none"}).  Can be a named integer vector to set
#'   different starts per chromosome.
#' @param end         Integer \strong{or} \code{NULL}.  Override the
#'   end position for every chromosome.  When \code{NULL} (default)
#'   the end is the chromosome length (adjusted for N-trimming).
#'   Can be a named integer vector.
#' @param strand      Character.  Strand label embedded in the
#'   coordinate string.  One of \code{"+"} (default) or \code{"-"}.
#' @param trim_N      Character.  N-base handling strategy.  One of
#'   \code{"ends"} (default), \code{"filter"}, or \code{"none"}.
#'   See section \emph{N-base trimming} above.
#' @param max_N_frac  Numeric in \eqn{[0, 1]}.  Maximum fraction of
#'   \code{N} bases tolerated per window.  Windows exceeding this
#'   threshold are dropped.  Only used when \code{trim_N = "filter"}.
#'   Default \code{0.1}.
#' @param N_scan_block Integer.  Block size (bp) for the coarse N-end
#'   detection scan used by \code{trim_N = "ends"}.  Larger values are
#'   faster but less precise.  Default \code{10000L}.
#' @param region_prefix Character.  Prefix for region IDs embedded in
#'   the coordinate string.  Default \code{"region"} (producing
#'   \code{region1}, \code{region2}, ...).
#' @param genome_pkg_name Character \strong{or} \code{NULL}.  The
#'   genome package name to embed in the coordinate string (4th field).
#'   When \code{NULL} (default), the name is extracted automatically
#'   from the \code{bsgenome} object via \code{S4Vectors::metadata(bsgenome)}.
#'   Supply explicitly when using a custom BSgenome object whose
#'   metadata name differs from the canonical package name.
#' @param as_vector   Logical.  When \code{TRUE} (default), return a
#'   plain named character vector.  When \code{FALSE}, return a
#'   \code{data.frame} with columns \code{coord}, \code{chr},
#'   \code{start}, \code{end}, \code{strand}, \code{genome},
#'   \code{region_id}, \code{chr_start_used}, \code{chr_end_used}
#'    --  useful for downstream GRanges construction.
#' @param verbose     Logical.  Print per-chromosome progress messages.
#'   Default \code{TRUE}.
#'
#' @return When \code{as_vector = TRUE} (default): a named character
#'   vector of coordinate strings, one element per window.  Names are
#'   the region IDs (\code{region1}, \code{region2}, ...).
#'
#'   When \code{as_vector = FALSE}: a \code{data.frame} with columns
#'   \code{coord}, \code{chr}, \code{win_start}, \code{win_end},
#'   \code{strand}, \code{genome}, \code{region_id},
#'   \code{chr_start_used}, \code{chr_end_used}.
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' ## -- Basic usage: tile chr1 with 200 bp windows, 50 bp slide ----------
#' coords <- make_genomiccoord(
#'   bsgenome    = BSgenome.Hsapiens.UCSC.hg38,
#'   chromosomes = "chr1",
#'   window      = 200L,
#'   slide       = 50L
#' )
#' length(coords)      # number of windows on chr1
#' head(coords, 3)
#' # region1 "chr1:10001-10200:+:BSgenome.Hsapiens.UCSC.hg38:region1"
#' # region2 "chr1:10051-10250:+:BSgenome.Hsapiens.UCSC.hg38:region2"
#' # region3 "chr1:10101-10300:+:BSgenome.Hsapiens.UCSC.hg38:region3"
#'
#' ## -- Non-overlapping tiling (slide == window) -------------------------
#' coords_nonoverlap <- make_genomiccoord(
#'   bsgenome    = BSgenome.Hsapiens.UCSC.hg38,
#'   chromosomes = paste0("chr", 1:22),
#'   window      = 200L,
#'   slide       = 200L    # no overlap
#' )
#' length(coords_nonoverlap)   # ~15 million windows across autosomes
#'
#' ## -- Custom start/end (e.g. a specific sub-region) --------------------
#' coords_sub <- make_genomiccoord(
#'   bsgenome    = BSgenome.Hsapiens.UCSC.hg38,
#'   chromosomes = "chr1",
#'   window      = 200L,
#'   slide       = 50L,
#'   start       = 1000000L,
#'   end         = 2000000L
#' )
#' length(coords_sub)   # 19,981 windows in 1 Mb region
#'
#' ## -- No N-trimming (use full chromosome length) ------------------------
#' coords_noN <- make_genomiccoord(
#'   bsgenome    = BSgenome.Hsapiens.UCSC.hg38,
#'   chromosomes = "chr1",
#'   window      = 200L,
#'   slide       = 50L,
#'   trim_N      = "none"
#' )
#'
#' ## -- Per-window N-filtering (removes windows with >10% N) -------------
#' coords_filt <- make_genomiccoord(
#'   bsgenome    = BSgenome.Hsapiens.UCSC.hg38,
#'   chromosomes = "chr1",
#'   window      = 200L,
#'   slide       = 50L,
#'   trim_N      = "filter",
#'   max_N_frac  = 0.10
#' )
#'
#' ## -- Get data.frame output for GRanges construction -------------------
#' df <- make_genomiccoord(
#'   bsgenome    = BSgenome.Hsapiens.UCSC.hg38,
#'   chromosomes = "chr1",
#'   window      = 200L,
#'   slide       = 50L,
#'   as_vector   = FALSE
#' )
#' gr <- GenomicRanges::GRanges(
#'   seqnames = df$chr,
#'   ranges   = IRanges::IRanges(start = df$win_start, end = df$win_end),
#'   strand   = df$strand
#' )
#'
#' ## -- Pass directly to tm_nn -------------------------------------------
#' coords <- make_genomiccoord(
#'   bsgenome    = BSgenome.Hsapiens.UCSC.hg38,
#'   chromosomes = "chr1",
#'   window      = 200L,
#'   slide       = 200L,
#'   start       = 1000000L,
#'   end         = 1010000L
#' )
#' tm_results <- tm_nn(coords, Na = 50)
#' }
#'
#' @author Junhui Li
#' @seealso \code{\link{tm_nn}}, \code{\link{tm_gc}},
#'   \code{\link{tm_wallace}}, \code{\link[BSgenome]{BSgenome}},
#'   \code{\link[Biostrings]{letterFrequency}}
#' @importFrom methods is
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom GenomeInfoDb seqlengths seqlevels genome
#' @importFrom BSgenome provider organism
#' @importFrom S4Vectors metadata
#' @export
make_genomiccoord <- function(
    bsgenome,
    chromosomes    = NULL,
    window         = 200L,
    slide          = 50L,
    start          = NULL,
    end            = NULL,
    strand         = "+",
    trim_N         = c("ends", "filter", "none"),
    max_N_frac     = 0.10,
    N_scan_block   = window,
    region_prefix  = "region",
    genome_pkg_name = NULL,
    as_vector      = TRUE,
    verbose        = TRUE
) {

  ## -- 0. Argument validation -----------------------------------------------
  trim_N <- match.arg(trim_N)

  window <- as.integer(window)
  slide  <- as.integer(slide)
  if (window < 1L) stop("'window' must be a positive integer.")
  if (slide  < 1L) stop("'slide' must be a positive integer.")
  if (slide  > window)
    warning("'slide' > 'window': windows will not be contiguous.",
            call. = FALSE)
  if (!strand %in% c("+", "-"))
    stop("'strand' must be \"+\" or \"-\".")
  if (max_N_frac < 0 || max_N_frac > 1)
    stop("'max_N_frac' must be in [0, 1].")

  ## -- 1. Resolve BSgenome object -------------------------------------------
  bsgenome <- .resolve_bsgenome(bsgenome)

  ## -- 2. Resolve genome package name for coordinate string -----------------
  pkg_name <- .resolve_pkg_name(bsgenome, genome_pkg_name)

  ## -- 3. Validate chromosomes -----------------------------------------------
  all_seqs   <- GenomeInfoDb::seqlevels(bsgenome)
  bad_chrs   <- setdiff(chromosomes, all_seqs)
  if (length(bad_chrs) > 0)
    stop("Chromosomes not found in BSgenome: ", paste(bad_chrs, collapse = ", "), "\nAvailable: ", paste(utils::head(all_seqs, 10), collapse = ", "), " ...")
  chr_lens <- GenomeInfoDb::seqlengths(bsgenome)[chromosomes]

  ## -- 4. Normalise start / end overrides -----------------------------------
  ## Allow scalar (applied to all chromosomes) or named vector
  start_override <- .normalise_override(start, chromosomes, "start")
  end_override   <- .normalise_override(end,   chromosomes, "end")

  ## -- 5. Main loop over chromosomes ----------------------------------------
  result_list <- vector("list", length(chromosomes))
  global_idx  <- 1L   # running region counter across all chromosomes

  for (ci in seq_along(chromosomes)) {

    chr     <- chromosomes[ci]
    chr_len <- chr_lens[chr]

    if (verbose)
      message(sprintf("[make_genomiccoord] Processing %s  (length = %s bp)", chr, format(chr_len, big.mark = ",")))

    ## -- 5a. Determine effective start / end for this chromosome ---------
    chr_start_raw <- if (!is.null(start_override)) start_override[chr] else 1L
    chr_end_raw   <- if (!is.null(end_override))   end_override[chr]   else chr_len

    chr_start_raw <- max(1L, as.integer(chr_start_raw))
    chr_end_raw   <- min(chr_len, as.integer(chr_end_raw))

    if (chr_start_raw >= chr_end_raw)
      stop("start >= end for chromosome ", chr, " (start=", chr_start_raw, ", end=", chr_end_raw, ")")

    ## -- 5b. N-end trimming ("ends" strategy) ----------------------------
    if (trim_N == "ends") {
      n_bounds <- .find_N_bounds(bsgenome, chr, chr_start_raw, chr_end_raw, N_scan_block)
      chr_start_eff <- n_bounds$first_nonN
      chr_end_eff   <- n_bounds$last_nonN

      if (verbose)
        message(sprintf("  N-trimming: effective range %s - %s", format(chr_start_eff, big.mark = ","), format(chr_end_eff,   big.mark = ",")))
    } else {
      chr_start_eff <- chr_start_raw
      chr_end_eff   <- chr_end_raw
    }

    if (chr_start_eff >= chr_end_eff - window + 1L) {
      warning("No valid windows on ", chr, " after N-trimming. Skipping.",
              call. = FALSE)
      result_list[[ci]] <- data.frame(
        coord = character(0), chr = character(0),
        win_start = integer(0), win_end = integer(0),
        strand = character(0), genome = character(0),
        region_id = character(0),
        chr_start_used = integer(0), chr_end_used = integer(0),
        stringsAsFactors = FALSE
      )
      next
    }

    ## -- 5c. Generate window start positions -----------------------------
    win_starts <- seq(from = chr_start_eff,
                      to   = chr_end_eff - window + 1L,
                      by   = slide)
    win_ends   <- win_starts + window - 1L

    ## Clip last window if it overhangs the effective end
    overhangs       <- win_ends > chr_end_eff
    win_ends[overhangs] <- chr_end_eff

    ## Drop any windows that became shorter than window after clipping
    ## (only happens at the very last window  --  keep it if >= 1 bp)
    valid_windows <- (win_ends - win_starts + 1L) >= 1L
    win_starts    <- win_starts[valid_windows]
    win_ends      <- win_ends[valid_windows]

    n_wins <- length(win_starts)
    if (n_wins == 0L) {
      warning("No valid windows on ", chr, ". Skipping.", call. = FALSE)
      result_list[[ci]] <- data.frame(
        coord = character(0), chr = character(0),
        win_start = integer(0), win_end = integer(0),
        strand = character(0), genome = character(0),
        region_id = character(0),
        chr_start_used = integer(0), chr_end_used = integer(0),
        stringsAsFactors = FALSE
      )
      next
    }

    ## -- 5d. Build coordinate strings ------------------------------------
    region_ids <- paste0(region_prefix, seq(global_idx, global_idx + n_wins - 1L))
    coords_chr <- sprintf("%s:%d-%d:%s:%s:%s",
                          chr, win_starts, win_ends, strand,
                          pkg_name, region_ids)

    ## -- 5e. Per-window N filtering ("filter" strategy) ------------------
    if (trim_N == "filter") {
      if (verbose)
        message(sprintf("  Filtering windows with N fraction > %.0f%%",
                        max_N_frac * 100))
      keep <- .filter_N_windows(bsgenome, chr, win_starts, win_ends,
                                 max_N_frac, window)
      if (verbose)
        message(sprintf("  %d / %d windows retained after N-filter",
                        sum(keep), n_wins))
      win_starts  <- win_starts[keep]
      win_ends    <- win_ends[keep]
      region_ids  <- region_ids[keep]
      coords_chr  <- coords_chr[keep]
      n_wins      <- sum(keep)

      ## Re-number region IDs to be globally contiguous
      region_ids <- paste0(region_prefix, seq(global_idx, global_idx + n_wins - 1L))
      coords_chr <- sprintf("%s:%d-%d:%s:%s:%s",
                            chr, win_starts, win_ends, strand,
                            pkg_name, region_ids)
    }

    if (verbose)
      message(sprintf("  Generated %s windows",
                      format(n_wins, big.mark = ",")))

    result_list[[ci]] <- data.frame(
      coord          = coords_chr,
      chr            = chr,
      win_start      = win_starts,
      win_end        = win_ends,
      strand         = strand,
      genome         = pkg_name,
      region_id      = region_ids,
      chr_start_used = chr_start_eff,
      chr_end_used   = chr_end_eff,
      stringsAsFactors = FALSE
    )

    global_idx <- global_idx + n_wins
  }

  ## -- 6. Combine and return ---------------------------------------------
  df_all <- do.call(rbind, result_list)
  rownames(df_all) <- NULL

  if (verbose)
    message(sprintf("[make_genomiccoord] Done. Total windows: %s",
                    format(nrow(df_all), big.mark = ",")))

  if (as_vector) {
    out <- df_all$coord
    names(out) <- df_all$region_id
    return(out)
  }

  df_all
}


# ============================================================================
# Internal helpers
# ============================================================================

#' @keywords internal
.resolve_bsgenome <- function(x) {
  if (is.character(x)) {
    if (!requireNamespace(x, quietly = TRUE))
      stop("BSgenome package '", x, "' is not installed.\n",
           "Install with: BiocManager::install(\"", x, "\")")
    return(get(x, envir = asNamespace(x)))
  }
  if (!methods::is(x, "BSgenome"))
    stop("'bsgenome' must be a BSgenome object or a package name string.")
  x
}


#' @keywords internal
.resolve_pkg_name <- function(bsgenome, user_name) {
  if (!is.null(user_name)) return(as.character(user_name))
  # Try metadata slot first (most reliable)
  meta <- tryCatch(
    S4Vectors::metadata(bsgenome),
    error = function(e) NULL
  )
  if (!is.null(meta) && "Package" %in% names(meta))
    return(meta[["Package"]])
  # Fallback: class name or provider + genome
  cls <- class(bsgenome)
  if (length(cls) > 0 && nchar(cls[1]) > 5) return(cls[1])
  paste0("BSgenome.", BSgenome::organism(bsgenome), ".",
         BSgenome::provider(bsgenome), ".", GenomeInfoDb::genome(bsgenome))
}


#' @keywords internal
.normalise_override <- function(val, chromosomes, arg_name) {
  if (is.null(val)) return(NULL)
  val <- as.integer(val)
  if (length(val) == 1L) {
    # Scalar -> replicate for all chromosomes
    out <- rep(val, length(chromosomes))
    names(out) <- chromosomes
    return(out)
  }
  if (is.null(names(val)))
    stop("When '", arg_name, "' is a vector with length > 1, it must be ",
         "a named integer vector with chromosome names as names.")
  missing_chrs <- setdiff(chromosomes, names(val))
  if (length(missing_chrs) > 0)
    stop("'", arg_name, "' does not cover all requested chromosomes. ",
         "Missing: ", paste(missing_chrs, collapse = ", "))
  val[chromosomes]
}


#' Detect the first and last non-N positions on a chromosome
#'
#' Uses a coarse block scan to find the positions of the first and last
#' non-N base on a chromosome. Efficient because it only reads
#' N_scan_block bp at a time from each end.
#'
#' @keywords internal
.find_N_bounds <- function(bsgenome, chr, chr_start, chr_end, block_size) {

  chr_seq_len <- chr_end - chr_start + 1L
  block_size  <- min(block_size, chr_seq_len)

  ## -- Scan from the LEFT to find first non-N position ---------------------
  first_nonN <- chr_start   # default: start right away
  scan_start <- chr_start

  while (scan_start <= chr_end) {
    scan_end <- min(scan_start + block_size - 1L, chr_end)
    block    <- Biostrings::getSeq(
      bsgenome,
      names  = chr,
      start  = scan_start,
      end    = scan_end,
      as.character = FALSE
    )
    ## Count non-N bases in the block
    acgt <- Biostrings::letterFrequency(block, letters = "ACGT", as.prob = FALSE)
    if (sum(acgt) > 0L) {
      ## Found non-N in this block  --  find the exact position
      block_str <- as.character(block)
      block_chars <- strsplit(block_str, "", fixed = TRUE)[[1]]
      rel_pos    <- which(block_chars != "N")[1L]
      first_nonN <- scan_start + rel_pos - 1L
      break
    }
    scan_start <- scan_end + 1L
  }

  ## -- Scan from the RIGHT to find last non-N position ----------------------
  last_nonN <- chr_end   # default: end right at chr_end
  scan_end2 <- chr_end

  while (scan_end2 >= chr_start) {
    scan_start2 <- max(scan_end2 - block_size + 1L, chr_start)
    block2 <- Biostrings::getSeq(
      bsgenome,
      names  = chr,
      start  = scan_start2,
      end    = scan_end2,
      as.character = FALSE
    )
    acgt2 <- Biostrings::letterFrequency(block2, letters = "ACGT", as.prob = FALSE)
    if (sum(acgt2) > 0L) {
      block_str2  <- as.character(block2)
      block_chars2 <- strsplit(block_str2, "", fixed = TRUE)[[1]]
      ## Last non-N: search from the end of the block
      rel_pos2   <- max(which(block_chars2 != "N"))
      last_nonN  <- scan_start2 + rel_pos2 - 1L
      break
    }
    scan_end2 <- scan_start2 - 1L
  }

  list(first_nonN = first_nonN, last_nonN = last_nonN)
}


#' Filter windows with too many N bases
#'
#' Reads each window sequence and drops those where the N fraction
#' exceeds max_N_frac. Uses batched getSeq for efficiency.
#'
#' @keywords internal
.filter_N_windows <- function(bsgenome, chr, win_starts, win_ends,
                               max_N_frac, window) {

  n_wins <- length(win_starts)
  if (n_wins == 0L) return(logical(0))

  ## Batch-fetch all window sequences in one getSeq call
  gr_query <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = win_starts, end = win_ends)
  )
  seqs <- Biostrings::getSeq(bsgenome, gr_query)  # DNAStringSet

  ## Count N bases per window
  n_counts  <- Biostrings::letterFrequency(seqs, letters = "N",
                                            as.prob = FALSE)[, 1]
  win_widths <- win_ends - win_starts + 1L
  n_frac     <- n_counts / win_widths

  n_frac <= max_N_frac
}
