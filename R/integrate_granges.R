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
#'   features map to the same Tm range / bin. Must accept a numeric vector and
#'   an \code{na.rm} argument (e.g. \code{mean}, \code{median}, \code{sum},
#'   \code{max}). Default: \code{mean}.
#' @param min_overlap Integer. Minimum overlap in base pairs required between
#'   a Tm range and a feature range in \code{"overlap"} mode. Default:
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
    min_overlap    = 1L,
    ignore_strand  = TRUE,
    keep_unmatched = TRUE,
    distance_col   = "distance_to_feature"
) {

  strategy <- match.arg(strategy)

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

  # -- Aggregation helper -----------------------------------------------------
  # Summarises a vector: mean (numeric) or comma-joined unique values (character)
  .agg <- function(vals) {
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA)
    if (is.numeric(vals)) agg_fun(vals, na.rm = TRUE)
    else paste(sort(unique(as.character(vals))), collapse = ",")
  }

  # Build a data.frame of aggregated feature columns given Hits object
  # query = gr_query, subject = gr_features (already subset to feature_cols)
  .aggregate_hits <- function(hits, n_query, gr_feat_sub) {
    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)
    meta  <- as.data.frame(GenomicRanges::mcols(gr_feat_sub)[, feature_cols,
                                                               drop = FALSE],
                            stringsAsFactors = FALSE)

    result <- lapply(seq_along(feature_cols), function(j) {
      col_vals <- meta[[j]]
      # Initialise with NA of the correct type
      if (is.numeric(col_vals)) out <- rep(NA_real_, n_query)
      else                      out <- rep(NA_character_, n_query)

      if (length(q_idx) == 0) return(out)

      agg_vals <- tapply(col_vals[s_idx], q_idx, .agg)
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

    feat_df <- .aggregate_hits(hits, length(gr_tm), gr_features)

    # Attach aggregated columns to gr_tm
    for (j in seq_along(out_names))
      GenomicRanges::mcols(gr_tm)[[out_names[j]]] <- feat_df[[out_names[j]]]

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
    feat_df <- .aggregate_hits(hits, length(gr_tm), gr_features)

    for (j in seq_along(out_names))
      GenomicRanges::mcols(gr_tm)[[out_names[j]]] <- feat_df[[out_names[j]]]

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

      if (length(q_tm) > 0) {
        agg_tm <- tapply(tm_vals[s_tm], q_tm, .agg)
        Tm_mean[as.integer(names(agg_tm))] <- as.numeric(agg_tm)
        n_tm <- tabulate(q_tm, nbins = n_bins)
      }
      GenomicRanges::mcols(bins)$Tm_mean     <- Tm_mean
      GenomicRanges::mcols(bins)$n_tm_ranges <- n_tm

      if (has_gc) {
        gc_vals  <- GenomicRanges::mcols(gr_tm_c)$GC
        GC_mean  <- rep(NA_real_, n_bins)
        if (length(q_tm) > 0) {
          agg_gc <- tapply(gc_vals[s_tm], q_tm, .agg)
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

        for (j in seq_along(feature_cols)) {
          col_vals <- meta[[j]]
          if (is.numeric(col_vals)) out <- rep(NA_real_,     n_bins)
          else                      out <- rep(NA_character_, n_bins)
          if (length(q_ft) > 0) {
            agg_vals <- tapply(col_vals[s_ft], q_ft, .agg)
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
