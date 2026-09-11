#' Integrate a Tm GRanges with multi-omic feature ranges
#'
#' @description
#' Combines the output of \code{\link{tm_calculate}} (a GRanges object with
#' \code{Tm} and \code{GC} columns) with a second GRanges carrying
#' arbitrary multi-omic metadata (ChIP-seq peaks, ATAC-seq signal, methylation
#' sites, gene annotations, etc.) using one of four positional strategies:
#'
#' \describe{
#'   \item{\code{"overlap"}}{Each tm range is annotated with the aggregated
#'     metadata of all feature ranges it directly overlaps.}
#'   \item{\code{"nearest"}}{Each tm range is annotated with the metadata of
#'     its single closest feature range, plus an added distance column.}
#'   \item{\code{"window"}}{Each tm range is expanded symmetrically by
#'     \code{window_size} bp and annotated with aggregated metadata from all
#'     features that fall within the expanded window.}
#'   \item{\code{"bin"}}{The genomic space covered by the data is tiled into
#'     equal-width bins. Each bin is annotated with the mean tm / GC of
#'     overlapping tm ranges \emph{and} the aggregated feature values - suitable
#'     for joint heatmaps and genome-wide correlation analyses.}
#' }
#'
#' For strategies \code{"overlap"} and \code{"window"}, when a single Tm range
#' matches multiple features the default behaviour is to \emph{summarise}:
#' numeric columns are aggregated via \code{agg_fun} (default \code{mean}),
#' and categorical columns are collapsed to a comma-separated string of unique
#' values.
#'
#' @param gr_tm A \code{GRanges} object produced by \code{tm_calculate()} (or
#'   \code{tm_calculate()$gr}). Must contain at least a \code{Tm} metadata
#'   column. A \code{gc} column is used automatically when present.
#' @param gr_features A \code{GRanges} object with multi-omic feature ranges.
#'   All (or a subset of) its metadata columns are transferred / aggregated.
#' @param strategy Character. Integration strategy. One of
#'   \code{"overlap"} (default), \code{"nearest"}, \code{"window"}, or
#'   \code{"bin"}.
#' @param feature_cols Character vector. Names of metadata columns in
#'   \code{gr_features} to transfer. \code{NULL} (default) transfers all
#'   metadata columns.
#' @param prefix Character. Prefix prepended to transferred column names to
#'   avoid clashes with existing columns in \code{gr_tm}. Default: \code{""}.
#'   Use e.g. \code{"feat_"} if there are naming conflicts.
#' @param window_size Integer. Half-width (bp) of the symmetric window added
#'   around each Tm range in \code{"window"} mode. Default: \code{1000}.
#' @param bin_size Integer. Width (bp) of genomic bins in \code{"bin"} mode.
#'   Default: \code{1e6} (1 Mb). Smaller values give finer resolution but
#'   sparser coverage.
#' @param agg_fun Function. Applied to numeric feature values when multiple
#'   features map to the same Tm range / bin. It is called as
#'   \code{agg_fun(values, na.rm = TRUE)}, so it must accept an
#'   \code{na.rm} argument; \code{mean}, \code{median}, \code{sum} and
#'   \code{max} all qualify, whereas \code{function(x) x[1]} does not.
#'   Ignored for character columns, which are always joined as a
#'   comma-separated list of their unique values. Default: \code{mean}.
#' @param weight Character. How features are combined within a range.
#'   \code{"none"} (default) gives every feature that passes the overlap test
#'   the same weight, whatever the length of its overlap. \code{"overlap"}
#'   computes a mean weighted by the number of overlapping base pairs,
#'   \eqn{\sum_j w_{ij} v_j / \sum_j w_{ij}} with \eqn{w_{ij}} the width of
#'   the intersection of range \eqn{i} and feature \eqn{j}.
#'
#'   The distinction matters wherever a signal changes sharply. A 200 bp
#'   window whose first 190 bp are covered at depth 5 and whose last 10 bp are
#'   covered at depth 200 has a true mean depth of 14.75; unweighted
#'   aggregation returns 102.5, because the 10 bp feature counts as much as
#'   the 190 bp one. Raising \code{min_overlap} does not fix this, it only
#'   reverses the sign of the bias by discarding the short feature entirely.
#'
#'   \code{"overlap"} requires \code{agg_fun = mean}: an overlap-weighted
#'   maximum has no accepted definition, and quietly ignoring \code{agg_fun}
#'   would return an unweighted value that looks weighted. It does not apply
#'   to \code{strategy = "nearest"}, which performs no aggregation.
#' @param report_coverage Logical. Add a \code{covered_frac} column giving the
#'   fraction of each range covered by at least one feature, computed after
#'   reducing the features so that overlapping ones are not counted twice.
#'   A weighted mean normalises by the covered bases, so a value derived from
#'   a quarter of a window is indistinguishable from one derived from all of
#'   it; this column is what makes the difference visible. Produced by the
#'   \code{"overlap"} and \code{"window"} strategies. Default: \code{FALSE}.
#' @param min_overlap Integer. Minimum overlap in base pairs required between
#'   a Tm range and a feature range. It is a filter and not a weight: a
#'   feature either qualifies or does not, and a qualifying feature counts in
#'   full. Applies to the \code{"overlap"} strategy only; \code{"window"} and
#'   \code{"bin"} require a single overlapping base pair. Default:
#'   \code{1}.
#' @param ignore_strand Logical. If \code{TRUE} (default), strand is ignored
#'   when finding overlaps / nearest neighbours.
#' @param keep_unmatched Logical. In \code{"overlap"} mode only: if
#'   \code{TRUE} (default) Tm ranges with no overlapping feature are retained
#'   with \code{NA} in the transferred columns. If \code{FALSE}, unmatched Tm
#'   ranges are dropped.
#' @param distance_col Character. Name of the distance column added in
#'   \code{"nearest"} mode. Default: \code{"distance_to_feature"}.
#'
#' @section Aggregating continuous signals:
#' When several features map to the same range, numeric columns are summarised
#' by \code{agg_fun} and character columns are joined as their unique values.
#' Which features take part is decided by \code{min_overlap}, and how much each
#' one counts is decided by \code{weight}. The two are easy to confuse, and the
#' distinction is what determines whether a coverage-like signal is summarised
#' correctly at a boundary.
#'
#' Let \eqn{R_i} be the \eqn{i}-th range of \code{gr_tm}, \eqn{F_j} the
#' \eqn{j}-th feature, \eqn{x_j} its value, and
#' \eqn{w_{ij} = |R_i \cap F_j|} the number of base pairs they share. The set
#' of features entering the summary of \eqn{R_i} is
#' \eqn{S_i = \{ j : w_{ij} \ge \code{min_overlap} \}}, and
#'
#' \deqn{v_i = \mathrm{agg\_fun}(\{x_j : j \in S_i\})}{
#'       v_i = agg_fun({x_j : j in S_i})}
#'
#' with \code{weight = "none"}, or
#'
#' \deqn{v_i = \frac{\sum_{j \in S_i} w_{ij} x_j}{\sum_{j \in S_i} w_{ij}}}{
#'       v_i = sum_j w_ij x_j / sum_j w_ij}
#'
#' with \code{weight = "overlap"}.
#'
#' The unweighted form gives a feature that overlaps by one base pair the same
#' influence as one that spans the whole range. This is harmless where a signal
#' is flat and wrong where it steps, which is to say at exon boundaries, peak
#' edges and promoters. Raising \code{min_overlap} does not repair it: the
#' threshold is a filter, so a short feature is either counted in full or
#' discarded in full, and the bias changes sign rather than disappearing. The
#' example below shows both failures against a case with a known answer.
#'
#' The weighted mean normalises by the covered base pairs, not by the width of
#' the range, so a value derived from a quarter of a window is
#' indistinguishable from one derived from all of it.
#' \code{report_coverage = TRUE} adds a \code{covered_frac} column giving the
#' fraction of each range covered by at least one feature, computed after
#' \code{\link[GenomicRanges]{reduce}}-ing the features so that overlapping
#' ones are not double counted. Multiply by it to convert a mean over covered
#' bases into a mean over the range.
#'
#' Weighting is not the default. Enabling it changes numeric output, and
#' existing analyses should stay reproducible unless their author decides
#' otherwise.
#'
#' @return
#' \itemize{
#'   \item \code{"overlap"}, \code{"nearest"}, \code{"window"}: A
#'     \code{GRanges} object with the same ranges as \code{gr_tm} (minus
#'     unmatched ranges if \code{keep_unmatched = FALSE}), with additional
#'     metadata columns from \code{gr_features}.
#'   \item \code{"bin"}: A new \code{GRanges} of genomic bins. Each bin
#'     carries \code{Tm_mean}, \code{GC_mean} (if available),
#'     \code{n_tm_ranges}, \code{n_features}, and one aggregated column per
#'     requested feature column.
#' }
#'
#' @examples
#' ## Aggregation: a coverage track with a known answer -----------------------
#' ## Two 200 bp windows over a signal that steps sharply inside the first.
#' ##
#' ##   window 1  [  1 .. 200]        window 2  [401 .. 600]
#' ##   depth     [  1 .. 190] = 5
#' ##             [191 .. 400] = 200
#' ##                                           [401 .. 450] = 60
#' library(GenomicRanges)
#'
#' win <- GRanges("chr1", IRanges(start = c(1, 401), width = 200),
#'                Tm = c(70, 72))
#' cov <- GRanges("chr1", IRanges(start = c(1, 191, 401),
#'                                end   = c(190, 400, 450)),
#'                cov = c(5, 200, 60))
#'
#' ## Window 1 truly averages (5 * 190 + 200 * 10) / 200 = 14.75.
#' integrate_granges(win, cov, strategy = "overlap")$cov
#' ## 102.5  60   the 10 bp feature counts as much as the 190 bp one
#'
#' integrate_granges(win, cov, strategy = "overlap",
#'                   weight = "overlap")$cov
#' ## 14.75  60   overlap-weighted mean recovers the true depth
#'
#' integrate_granges(win, cov, strategy = "overlap", min_overlap = 20L)$cov
#' ## 5      60   the threshold discards the short feature; bias reverses
#'
#' ## Window 2 is only a quarter covered. Weighting cannot show that, because
#' ## it normalises by the covered bases; the coverage column can.
#' res <- integrate_granges(win, cov, strategy = "overlap",
#'                          weight = "overlap", report_coverage = TRUE)
#' res$cov                      # 14.75  60
#' res$covered_frac             # 1.00   0.25
#' res$cov * res$covered_frac   # 14.75  15    mean over the whole window
#'
#' \dontrun{
#' library(GenomicRanges)
#'
#' # -- Sample data ----------------------------------------------------------
#' set.seed(42)
#' gr_tm <- GRanges(
#'   seqnames = c(rep("chr1", 60), rep("chr2", 30)),
#'   ranges   = IRanges(
#'     start = c(sort(sample(1:249e6, 60)),
#'               sort(sample(1:243e6, 30))),
#'     width = sample(50:200, 90, replace = TRUE)
#'   ),
#'   Tm = runif(90, 55, 85),
#'   GC = runif(90, 30, 70)
#' )
#'
#' gr_features <- GRanges(
#'   seqnames = c(rep("chr1", 40), rep("chr2", 20)),
#'   ranges   = IRanges(
#'     start = c(sort(sample(1:249e6, 40)),
#'               sort(sample(1:243e6, 20))),
#'     width = sample(500:5000, 60, replace = TRUE)
#'   ),
#'   score      = runif(60, 0, 100),
#'   peak_type  = sample(c("narrow", "broad"), 60, replace = TRUE),
#'   signal     = rnorm(60, 5, 2)
#' )
#'
#' # Strategy 1: overlap - annotate Tm ranges with overlapping peak features
#' res_overlap <- integrate_granges(gr_tm, gr_features,
#'                                   strategy = "overlap")
#'
#' # Strategy 2: nearest - every Tm range gets its closest peak + distance
#' res_nearest <- integrate_granges(gr_tm, gr_features,
#'                                   strategy = "nearest")
#' head(res_nearest$distance_to_feature)
#'
#' # Strategy 3: window - 5 kb window around each probe
#' res_window <- integrate_granges(gr_tm, gr_features,
#'                                  strategy = "window", window_size = 5000)
#'
#' # Strategy 4: bin - 500 kb genome bins with mean Tm and aggregated signal
#' res_bin <- integrate_granges(gr_tm, gr_features,
#'                               strategy = "bin", bin_size = 5e5)
#' as.data.frame(res_bin) |> head()
#'
#' # Use a subset of feature columns and add a prefix
#' integrate_granges(gr_tm, gr_features,
#'                   strategy    = "overlap",
#'                   feature_cols = c("score", "peak_type"),
#'                   prefix       = "chip_")
#' }
#'
#' @importFrom GenomicRanges findOverlaps distanceToNearest resize tile GRanges
#'   mcols seqnames start end width
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @encoding UTF-8
#' @author Junhui Li
#' @export

integrate_granges <- function(
    gr_tm,
    gr_features,
    strategy       = c("overlap", "nearest", "window", "bin"),
    feature_cols   = NULL,
    prefix         = "",
    window_size    = 1000L,
    bin_size       = 1e6,
    agg_fun        = mean,
    weight         = c("none", "overlap"),
    report_coverage = FALSE,
    min_overlap    = 1L,
    ignore_strand  = TRUE,
    keep_unmatched = TRUE,
    distance_col   = "distance_to_feature"
) {

  strategy <- match.arg(strategy)
  weight   <- match.arg(weight)
  weighted <- identical(weight, "overlap")

  # Weighting is defined for a mean and for nothing else: there is no
  # canonical weighted maximum, and a weighted median needs a definition
  # chosen rather than assumed. Silently ignoring agg_fun would be the worst
  # outcome, since the result would look weighted and would not be.
  if (weighted && !identical(agg_fun, mean))
    stop("weight = \"overlap\" computes an overlap-length-weighted mean and ",
         "cannot honour a different `agg_fun`. Use agg_fun = mean, or ",
         "weight = \"none\".", call. = FALSE)
  if (weighted && strategy == "nearest")
    stop("weight = \"overlap\" does not apply to strategy = \"nearest\", ",
         "which transfers one feature per range and performs no aggregation.",
         call. = FALSE)

  # -- Input validation -------------------------------------------------------
  if (!inherits(gr_tm, "GRanges"))
    stop("'gr_tm' must be a GRanges object.")
  if (!inherits(gr_features, "GRanges"))
    stop("'gr_features' must be a GRanges object.")
  if (!"Tm" %in% names(GenomicRanges::mcols(gr_tm)))
    stop("'gr_tm' must contain a 'Tm' metadata column (output of tm_calculate()).")

  all_feat_cols <- names(GenomicRanges::mcols(gr_features))
  if (is.null(feature_cols)) {
    feature_cols <- all_feat_cols
  } else {
    bad <- setdiff(feature_cols, all_feat_cols)
    if (length(bad) > 0)
      stop(sprintf("feature_cols not found in gr_features: %s",
                   paste(bad, collapse = ", ")))
  }
  if (length(feature_cols) == 0)
    stop("gr_features has no metadata columns to transfer.")

  # Prefixed output column names
  out_names <- paste0(prefix, feature_cols)

  # -- Aggregation helpers ----------------------------------------------------
  # Overlap width for each (query, subject) pair of a Hits object. Computed
  # arithmetically rather than with pintersect(): findOverlaps() already
  # returns parallel indices, so the widths are one vectorised subtraction and
  # no intermediate GRanges of the same length has to be allocated, which
  # matters at genome scale.
  .ov_width <- function(gq, gs, qi, si)
    pmin(GenomicRanges::end(gq)[qi],   GenomicRanges::end(gs)[si]) -
    pmax(GenomicRanges::start(gq)[qi], GenomicRanges::start(gs)[si]) + 1L

  # Summarises the values mapped to one range. Character columns are joined
  # regardless of `weight`, since a weighted mean of labels is meaningless;
  # this is documented rather than signalled, because a table of features
  # normally carries both kinds of column and refusing the whole call for one
  # of them would be unhelpful.
  .agg <- function(vals, wts = NULL) {
    keep <- !is.na(vals)
    if (!any(keep)) return(NA)
    vals <- vals[keep]
    if (!is.numeric(vals))
      return(paste(sort(unique(as.character(vals))), collapse = ","))
    if (weighted) {
      wts <- as.numeric(wts[keep])
      tot <- sum(wts)
      if (!is.finite(tot) || tot <= 0) return(NA_real_)
      sum(vals * wts) / tot
    } else agg_fun(vals, na.rm = TRUE)
  }

  # Fraction of a range covered by at least one feature. The features are
  # reduced first, so overlapping features are not counted twice and the
  # result cannot exceed one; without that step the quantity would not be a
  # fraction at all.
  .covered_fraction <- function(gr_query, gr_feat) {
    red <- GenomicRanges::reduce(gr_feat, ignore.strand = ignore_strand)
    h   <- GenomicRanges::findOverlaps(gr_query, red,
                                       ignore.strand = ignore_strand)
    out <- rep(0, length(gr_query))
    if (length(h) == 0L) return(out)
    qi <- S4Vectors::queryHits(h); si <- S4Vectors::subjectHits(h)
    w  <- .ov_width(gr_query, red, qi, si)
    tot <- tapply(w, qi, sum)
    out[as.integer(names(tot))] <- as.numeric(tot)
    pmin(out / GenomicRanges::width(gr_query), 1)
  }

  # Build a data.frame of aggregated feature columns given Hits object
  # query = gr_query, subject = gr_features (already subset to feature_cols)
  .aggregate_hits <- function(hits, n_query, gr_feat_sub, gr_query = NULL) {
    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)
    ov_w  <- if (weighted && length(q_idx))
      .ov_width(gr_query, gr_feat_sub, q_idx, s_idx) else NULL
    meta  <- as.data.frame(GenomicRanges::mcols(gr_feat_sub)[, feature_cols,
                                                               drop = FALSE],
                            stringsAsFactors = FALSE)

    result <- lapply(seq_along(feature_cols), function(j) {
      col_vals <- meta[[j]]
      # Initialise with NA of the correct type
      if (is.numeric(col_vals)) out <- rep(NA_real_, n_query)
      else                      out <- rep(NA_character_, n_query)

      if (length(q_idx) == 0) return(out)

      agg_vals <- tapply(seq_along(q_idx), q_idx, function(i)
        .agg(col_vals[s_idx[i]], ov_w[i]))
      idx      <- as.integer(names(agg_vals))
      out[idx] <- as.vector(agg_vals)
      out
    })
    names(result) <- out_names
    as.data.frame(result, stringsAsFactors = FALSE)
  }


  # ════════════════════════════════════════════════════════════════════════════
  # STRATEGY 1: OVERLAP
  # ════════════════════════════════════════════════════════════════════════════
  if (strategy == "overlap") {

    hits <- GenomicRanges::findOverlaps(
      gr_tm, gr_features,
      minoverlap    = min_overlap,
      ignore.strand = ignore_strand
    )

    feat_df <- .aggregate_hits(hits, length(gr_tm), gr_features, gr_tm)

    # Attach aggregated columns to gr_tm
    for (j in seq_along(out_names))
      GenomicRanges::mcols(gr_tm)[[out_names[j]]] <- feat_df[[out_names[j]]]

    # Reported alongside the summary rather than instead of it: a weighted
    # mean normalises by the covered bases, so a value computed from a
    # quarter of a window looks exactly like one computed from all of it.
    if (report_coverage)
      GenomicRanges::mcols(gr_tm)[[paste0(prefix, "covered_frac")]] <-
        .covered_fraction(gr_tm, gr_features)

    if (!keep_unmatched) {
      matched <- unique(S4Vectors::queryHits(hits))
      gr_tm   <- gr_tm[matched]
    }

    return(gr_tm)
  }


  # ════════════════════════════════════════════════════════════════════════════
  # STRATEGY 2: NEAREST
  # ════════════════════════════════════════════════════════════════════════════
  if (strategy == "nearest") {

    nn <- GenomicRanges::distanceToNearest(
      gr_tm, gr_features,
      ignore.strand = ignore_strand
    )

    q_idx    <- S4Vectors::queryHits(nn)
    s_idx    <- S4Vectors::subjectHits(nn)
    dist_vec <- S4Vectors::mcols(nn)$distance

    # Add distance column (NA for any Tm range that had no feature on same chr)
    d_out <- rep(NA_integer_, length(gr_tm))
    d_out[q_idx] <- dist_vec
    GenomicRanges::mcols(gr_tm)[[distance_col]] <- d_out

    # Transfer one metadata value per Tm range (no aggregation needed)
    meta <- as.data.frame(GenomicRanges::mcols(gr_features)[, feature_cols,
                                                              drop = FALSE],
                           stringsAsFactors = FALSE)
    for (j in seq_along(feature_cols)) {
      col_vals <- meta[[j]]
      if (is.numeric(col_vals)) out <- rep(NA_real_,      length(gr_tm))
      else                      out <- rep(NA_character_,  length(gr_tm))
      out[q_idx] <- col_vals[s_idx]
      GenomicRanges::mcols(gr_tm)[[out_names[j]]] <- out
    }

    return(gr_tm)
  }


  # ════════════════════════════════════════════════════════════════════════════
  # STRATEGY 3: WINDOW
  # ════════════════════════════════════════════════════════════════════════════
  if (strategy == "window") {

    # Expand each Tm range symmetrically by window_size bp
    gr_expanded <- GenomicRanges::resize(
      gr_tm,
      width = GenomicRanges::width(gr_tm) + 2L * as.integer(window_size),
      fix   = "center"
    )
    # Clamp to chromosome start (avoid negative coordinates)
    GenomicRanges::start(gr_expanded) <- pmax(1L,
                                               GenomicRanges::start(gr_expanded))

    hits    <- GenomicRanges::findOverlaps(
      gr_expanded, gr_features,
      ignore.strand = ignore_strand
    )
    feat_df <- .aggregate_hits(hits, length(gr_tm), gr_features, gr_expanded)

    for (j in seq_along(out_names))
      GenomicRanges::mcols(gr_tm)[[out_names[j]]] <- feat_df[[out_names[j]]]

    # Relative to the EXPANDED range, which is what the overlaps were found
    # against; a fraction of the original range would not describe the
    # summary sitting beside it.
    if (report_coverage)
      GenomicRanges::mcols(gr_tm)[[paste0(prefix, "covered_frac")]] <-
        .covered_fraction(gr_expanded, gr_features)

    return(gr_tm)
  }


  # ════════════════════════════════════════════════════════════════════════════
  # STRATEGY 4: BIN
  # ════════════════════════════════════════════════════════════════════════════
  if (strategy == "bin") {

    has_gc  <- "GC" %in% names(GenomicRanges::mcols(gr_tm))
    chrs    <- intersect(
      unique(as.character(GenomicRanges::seqnames(gr_tm))),
      unique(as.character(GenomicRanges::seqnames(gr_features)))
    )
    if (length(chrs) == 0) {
      chrs <- unique(as.character(GenomicRanges::seqnames(gr_tm)))
      message("No shared chromosomes between gr_tm and gr_features; ",
              "bins will have NA for feature columns.")
    }

    bin_list <- lapply(chrs, function(chr_i) {

      gr_tm_c  <- gr_tm[GenomicRanges::seqnames(gr_tm)      == chr_i]
      gr_ft_c  <- gr_features[GenomicRanges::seqnames(gr_features) == chr_i]

      # Bin extent: union of data ranges on this chromosome
      chr_from <- min(GenomicRanges::start(gr_tm_c))
      chr_to   <- max(GenomicRanges::end(gr_tm_c))
      if (length(gr_ft_c) > 0) {
        chr_from <- min(chr_from, min(GenomicRanges::start(gr_ft_c)))
        chr_to   <- max(chr_to,   max(GenomicRanges::end(gr_ft_c)))
      }

      # Create bins via tile()
      chr_range <- GenomicRanges::GRanges(
        seqnames = chr_i,
        ranges   = IRanges::IRanges(start = chr_from, end = chr_to)
      )
      bins <- unlist(GenomicRanges::tile(chr_range, width = as.integer(bin_size)))

      n_bins <- length(bins)

      # -- Aggregate Tm (and GC) per bin -----------------------------------
      tm_hits  <- GenomicRanges::findOverlaps(bins, gr_tm_c,
                                               ignore.strand = ignore_strand)
      q_tm     <- S4Vectors::queryHits(tm_hits)
      s_tm     <- S4Vectors::subjectHits(tm_hits)
      tm_vals  <- GenomicRanges::mcols(gr_tm_c)$Tm

      Tm_mean    <- rep(NA_real_, n_bins)
      n_tm       <- integer(n_bins)

      ov_tm <- if (weighted && length(q_tm))
        .ov_width(bins, gr_tm_c, q_tm, s_tm) else NULL

      if (length(q_tm) > 0) {
        agg_tm <- tapply(seq_along(q_tm), q_tm, function(i)
          .agg(tm_vals[s_tm[i]], ov_tm[i]))
        Tm_mean[as.integer(names(agg_tm))] <- as.numeric(agg_tm)
        n_tm <- tabulate(q_tm, nbins = n_bins)
      }
      GenomicRanges::mcols(bins)$Tm_mean     <- Tm_mean
      GenomicRanges::mcols(bins)$n_tm_ranges <- n_tm

      if (has_gc) {
        gc_vals  <- GenomicRanges::mcols(gr_tm_c)$GC
        GC_mean  <- rep(NA_real_, n_bins)
        if (length(q_tm) > 0) {
          agg_gc <- tapply(seq_along(q_tm), q_tm, function(i)
            .agg(gc_vals[s_tm[i]], ov_tm[i]))
          GC_mean[as.integer(names(agg_gc))] <- as.numeric(agg_gc)
        }
        GenomicRanges::mcols(bins)$GC_mean <- GC_mean
      }

      # -- Aggregate feature columns per bin -------------------------------
      n_feat_col    <- rep(0L, n_bins)

      if (length(gr_ft_c) > 0) {
        feat_hits <- GenomicRanges::findOverlaps(bins, gr_ft_c,
                                                  ignore.strand = ignore_strand)
        q_ft <- S4Vectors::queryHits(feat_hits)
        s_ft <- S4Vectors::subjectHits(feat_hits)
        meta <- as.data.frame(
          GenomicRanges::mcols(gr_ft_c)[, feature_cols, drop = FALSE],
          stringsAsFactors = FALSE
        )

        ov_ft <- if (weighted && length(q_ft))
          .ov_width(bins, gr_ft_c, q_ft, s_ft) else NULL

        for (j in seq_along(feature_cols)) {
          col_vals <- meta[[j]]
          if (is.numeric(col_vals)) out <- rep(NA_real_,     n_bins)
          else                      out <- rep(NA_character_, n_bins)
          if (length(q_ft) > 0) {
            agg_vals <- tapply(seq_along(q_ft), q_ft, function(i)
              .agg(col_vals[s_ft[i]], ov_ft[i]))
            out[as.integer(names(agg_vals))] <- as.vector(agg_vals)
          }
          GenomicRanges::mcols(bins)[[out_names[j]]] <- out
        }

        n_feat_col <- tabulate(q_ft, nbins = n_bins)
      } else {
        # No features on this chromosome: fill with NA
        for (j in seq_along(feature_cols)) {
          col_vals <- GenomicRanges::mcols(gr_features)[[feature_cols[j]]]
          GenomicRanges::mcols(bins)[[out_names[j]]] <-
            if (is.numeric(col_vals)) rep(NA_real_, n_bins)
            else                      rep(NA_character_, n_bins)
        }
      }
      GenomicRanges::mcols(bins)$n_features <- n_feat_col

      bins
    })

    result <- do.call(c, bin_list)
    return(result)
  }
}
