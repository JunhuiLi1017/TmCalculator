utils::globalVariables(c("group", "value", "rank"))

#' Compare Tm distributions across groups
#'
#' Visualise how melting temperature (and optionally other numeric metadata)
#' differs between region classes such as peak vs. non-peak, mutant vs.
#' wild-type, or any categorical annotation stored in a GRanges object.
#' Designed as the visual companion to \code{\link{compare_groups}()}.
#'
#' @param gr A \code{GRanges} object containing at least a numeric \code{Tm}
#'   metadata column and a categorical grouping column.
#' @param group Character. Name of the metadata column that defines the groups
#'   to compare (e.g. \code{"in_mutH"}).
#' @param value Character. Numeric metadata column to plot on the y-axis.
#'   Default: \code{"Tm"}. Can be any numeric column such as \code{"GC"}.
#' @param plot_type Character. Visualisation style:
#'   \describe{
#'     \item{\code{"box"}}{Box plot (default). Shows median, quartiles, and
#'       outliers per group.}
#'     \item{\code{"violin"}}{Violin plot with embedded box plot. Shows the
#'       full density shape of each group.}
#'     \item{\code{"rank"}}{Rank plot. All values sorted ascending on the
#'       x-axis, coloured by group. Visually demonstrates whether one group
#'       clusters at higher or lower values — a graphical Wilcoxon test.}
#'     \item{\code{"density"}}{Overlaid density curves per group.}
#'     \item{\code{"ridgeline"}}{Stacked density ridges per group (requires
#'       \pkg{ggridges}).}
#'     \item{\code{"ecdf"}}{Empirical cumulative distribution function per
#'       group. Useful for visualising distributional shifts.}
#'     \item{\code{"sina"}}{Sina (strip) plot — jittered points shaped by
#'       the local density, combining the resolution of a strip chart with the
#'       density information of a violin plot.}
#'   }
#' @param color_palette Character. Viridis palette for group colours:
#'   \code{"viridis"} (default), \code{"magma"}, \code{"plasma"},
#'   \code{"inferno"}, or \code{"cividis"}.
#' @param show_points Logical. Overlay jittered points on box / violin plots.
#'   Default: \code{FALSE}. Ignored for rank, density, ridgeline, ecdf.
#' @param point_size Numeric. Point size. Default: \code{1}.
#' @param point_alpha Numeric in \code{[0, 1]}. Point transparency.
#'   Default: \code{0.3}.
#' @param notch Logical. Draw notched box plots (approximate 95\% CI for
#'   median). Default: \code{FALSE}.
#' @param add_mean Logical. Add a diamond marker at the group mean on box /
#'   violin plots. Default: \code{FALSE}.
#' @param facet_by Character or \code{NULL}. Optional metadata column to facet
#'   the plot by (e.g. \code{"chromosome"}). Default: \code{NULL}.
#' @param title Character or \code{NULL}. Plot title. A sensible default is
#'   generated when \code{NULL}.
#' @param ylab Character. Y-axis label. Default: \code{"Tm (\u00B0C)"}.
#' @param xlab Character. X-axis label. Default: group column name for
#'   box/violin/sina, \code{"Rank"} for rank, \code{value} for density/ecdf.
#' @param show_pvalue Logical. Annotate the plot with Wilcoxon rank-sum test
#'   p-values. For two groups a single p-value is shown; for three or more
#'   groups all pairwise comparisons are displayed. P-values are placed as
#'   bracket annotations on box / violin / sina plots, or as a subtitle on
#'   density / ecdf / rank / ridgeline plots. Default: \code{FALSE}.
#' @param p_adjust_method Character. Method for p-value adjustment when there
#'   are more than two groups (passed to \code{\link[stats]{p.adjust}}).
#'   Default: \code{"BH"} (Benjamini-Hochberg).
#' @param show_legend Logical. Show colour legend. Default: \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The function pairs naturally with \code{\link{integrate_granges}()} and
#' \code{\link{compare_groups}()}:
#'
#' \enumerate{
#'   \item Use \code{integrate_granges()} to annotate Tm windows with feature
#'     overlaps (e.g. ChIP-seq peaks, repeat classes).
#'   \item Use \code{compare_groups()} for the statistical test.
#'   \item Use \code{plot_tm()} to visualise the distributional difference.
#' }
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' data(ecoli_rep_hotspots)
#' 
#' library("BSgenome.Ecoli.NCBI.ASM584v2")
#' genome_name <- "BSgenome.Ecoli.NCBI.ASM584v2"
#' chr_name    <- "U00096.3"
#' genome <- get(genome_name, envir = asNamespace(genome_name))
#' chr_length  <- length(genome[[chr_name]])

#' genome_name="BSgenome.Ecoli.NCBI.ASM584v2"
#' bins_gc <- make_genomiccoord(
#'   bsgenome    = genome_name,
#'   chromosomes = chr_name,
#'   window      = 200L,
#'   slide       = 200L,
#'   start       = 1,
#'   end         = chr_length,
#'   strand      = "+"
#' )
#' input_new <- list(pkg_name = genome_name, seq = bins_gc)
#' gr_batch <- to_genomic_ranges_fast(input_new)

#' tm_ASM584v2 <- tm_calculate(
#'   gr_batch,
#'   method   = "tm_nn"
#' )
#' gr_tm <- tm_ASM584v2$gr
#'
#' # Annotate with MutL-AR peak membership
#' mutH_peaks <- GRanges(
#'   seqnames = ecoli_rep_hotspots$all_peaks_IP_mutH$chr,
#'   ranges   = IRanges(start = ecoli_rep_hotspots$all_peaks_IP_mutH$start,
#'                      end   = ecoli_rep_hotspots$all_peaks_IP_mutH$end)
#' )
#' mutH_peaks$peak_id <- paste0("mutH_", seq_along(mutH_peaks))
#' gr_annot <- integrate_granges(
#'   gr_tm = gr_tm, gr_features = mutH_peaks,
#'   strategy = "overlap", feature_cols = "peak_id",
#'   keep_unmatched = TRUE
#' )
#' gr_annot$in_mutH <- ifelse(is.na(gr_annot$peak_id), "non_peak", "peak")
#'
#' # Box plot — peak vs non-peak Tm
#' plot_tm(gr_annot, group = "in_mutH")
#'
#' # Box plot with Wilcoxon p-value bracket
#' plot_tm(gr_annot, group = "in_mutH", show_pvalue = TRUE)
#'
#' # Violin plot with p-value and jittered points
#' plot_tm(gr_annot, group = "in_mutH", plot_type = "violin",
#'         show_points = TRUE, show_pvalue = TRUE)
#'
#' # Rank plot — visual Wilcoxon test with p-value subtitle
#' plot_tm(gr_annot, group = "in_mutH", plot_type = "rank",
#'         show_pvalue = TRUE)
#'
#' # Density overlay with p-value
#' plot_tm(gr_annot, group = "in_mutH", plot_type = "density",
#'         show_pvalue = TRUE)
#'
#' # ECDF — cumulative distribution comparison
#' plot_tm(gr_annot, group = "in_mutH", plot_type = "ecdf",
#'         show_pvalue = TRUE)
#'
#' # Compare GC content instead of Tm
#' plot_tm(gr_annot, group = "in_mutH", value = "GC", ylab = "GC content",
#'         show_pvalue = TRUE)
#'
#' # Notched box plot with mean markers and p-value
#' plot_tm(gr_annot, group = "in_mutH", notch = TRUE, add_mean = TRUE,
#'         show_pvalue = TRUE)
#' }
#'
#' @importFrom stats wilcox.test pairwise.wilcox.test
#' @importFrom GenomicRanges mcols seqnames start end
#'
#' @export

plot_tm <- function(
    gr,
    group,
    value          = "Tm",
    plot_type      = c("box", "violin", "rank", "density",
                        "ridgeline", "ecdf", "sina"),
    color_palette  = c("viridis", "magma", "plasma", "inferno", "cividis"),
    show_points    = FALSE,
    point_size     = 1,
    point_alpha    = 0.3,
    notch          = FALSE,
    add_mean       = FALSE,
    facet_by       = NULL,
    show_pvalue    = FALSE,
    p_adjust_method = "BH",
    title          = NULL,
    ylab           = "Tm (\u00B0C)",
    xlab           = NULL,
    show_legend    = TRUE
) {

  # -- Check optional dependencies ---------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_tm(). ",
         "Install it with: install.packages(\"ggplot2\")", call. = FALSE)

  # -- Input validation -------------------------------------------------------
  plot_type     <- match.arg(plot_type)
  color_palette <- match.arg(color_palette)

  if (!inherits(gr, "GRanges"))
    stop("'gr' must be a GRanges object.")

  meta_cols <- names(GenomicRanges::mcols(gr))

  if (!value %in% meta_cols)
    stop(sprintf("'%s' not found in metadata columns of gr.", value))
  if (!is.numeric(GenomicRanges::mcols(gr)[[value]]))
    stop(sprintf("'%s' must be a numeric metadata column.", value))
  if (!group %in% meta_cols)
    stop(sprintf("group column '%s' not found in metadata columns of gr.", group))
  if (!is.null(facet_by) && !facet_by %in% c("chromosome", meta_cols))
    stop(sprintf("facet_by column '%s' not found.", facet_by))

  # -- Build data.frame -------------------------------------------------------
  df <- data.frame(
    value      = GenomicRanges::mcols(gr)[[value]],
    group      = as.character(GenomicRanges::mcols(gr)[[group]]),
    chromosome = as.character(GenomicRanges::seqnames(gr)),
    stringsAsFactors = FALSE
  )

  if (!is.null(facet_by) && facet_by != "chromosome" && facet_by != group) {
    df[[facet_by]] <- as.character(GenomicRanges::mcols(gr)[[facet_by]])
  }

  df <- df[!is.na(df$value) & !is.na(df$group), ]
  if (nrow(df) == 0) stop("No non-NA data to plot.")

  # -- Group palette -----------------------------------------------------------
  grp_levels <- sort(unique(df$group))
  n_grp      <- length(grp_levels)
  grp_colors <- stats::setNames(
    grDevices::hcl.colors(n_grp, palette = paste0(toupper(substring(color_palette, 1, 1)),
                                                        substring(color_palette, 2))),
    grp_levels
  )
  df$group <- factor(df$group, levels = grp_levels)

  # -- Defaults ----------------------------------------------------------------
  if (is.null(title)) {
    title <- sprintf("%s by %s", value, group)
  }
  if (is.null(xlab)) {
    xlab <- switch(plot_type,
      rank    = "Rank (ascending)",
      density = value,
      ecdf    = value,
      ridgeline = value,
      group  # box, violin, sina
    )
  }

  # -- Shared summary layer ----------------------------------------------------
  mean_layer <- if (add_mean && plot_type %in% c("box", "violin", "sina")) {
    ggplot2::stat_summary(
      fun = mean, geom = "point", shape = 18, size = 3, color = "black",
      show.legend = FALSE
    )
  }

  # -- Build plot by type ------------------------------------------------------
  p <- switch(plot_type,

    # ---- BOX ----------------------------------------------------------------
    box = {
      pp <- ggplot2::ggplot(df, ggplot2::aes(
        x = group, y = value, fill = group
      )) +
        ggplot2::geom_boxplot(
          notch = notch, alpha = 0.7, outlier.size = 0.8,
          outlier.alpha = 0.5
        )
      if (show_points)
        pp <- pp + ggplot2::geom_jitter(
          width = 0.2, size = point_size, alpha = point_alpha
        )
      pp
    },

    # ---- VIOLIN -------------------------------------------------------------
    violin = {
      pp <- ggplot2::ggplot(df, ggplot2::aes(
        x = group, y = value, fill = group
      )) +
        ggplot2::geom_violin(alpha = 0.7, trim = FALSE) +
        ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.8,
                              outlier.shape = NA)
      if (show_points)
        pp <- pp + ggplot2::geom_jitter(
          width = 0.2, size = point_size, alpha = point_alpha
        )
      pp
    },

    # ---- RANK ---------------------------------------------------------------
    rank = {
      df <- df[order(df$value), ]
      df$rank <- seq_len(nrow(df))
      ggplot2::ggplot(df, ggplot2::aes(
        x = rank, y = value, color = group
      )) +
        ggplot2::geom_point(size = point_size, alpha = 0.6)
    },

    # ---- DENSITY ------------------------------------------------------------
    density = {
      ggplot2::ggplot(df, ggplot2::aes(
        x = value, fill = group, color = group
      )) +
        ggplot2::geom_density(alpha = 0.35, linewidth = 0.6)
    },

    # ---- RIDGELINE ----------------------------------------------------------
    ridgeline = {
      if (!requireNamespace("ggridges", quietly = TRUE))
        stop("Package 'ggridges' is required for ridgeline plots. ",
             "Install it with: install.packages(\"ggridges\")", call. = FALSE)
      ggplot2::ggplot(df, ggplot2::aes(
        x = value, y = group, fill = group
      )) +
        ggridges::geom_density_ridges(alpha = 0.6, scale = 1.2)
    },

    # ---- ECDF ---------------------------------------------------------------
    ecdf = {
      ggplot2::ggplot(df, ggplot2::aes(
        x = value, color = group
      )) +
        ggplot2::stat_ecdf(linewidth = 0.8)
    },

    # ---- SINA ---------------------------------------------------------------
    sina = {
      if (!requireNamespace("ggforce", quietly = TRUE))
        stop("Package 'ggforce' is required for sina plots. ",
             "Install it with: install.packages(\"ggforce\")", call. = FALSE)
      pp <- ggplot2::ggplot(df, ggplot2::aes(
        x = group, y = value, color = group
      )) +
        ggforce::geom_sina(size = point_size, alpha = 0.5)
      pp
    }
  )

  # -- Mean markers (box / violin / sina) --------------------------------------
  if (!is.null(mean_layer)) p <- p + mean_layer

  # -- Color / fill scales -----------------------------------------------------
  if (plot_type %in% c("box", "violin", "density", "ridgeline")) {
    p <- p + ggplot2::scale_fill_manual(values = grp_colors, name = group)
  }
  if (plot_type %in% c("rank", "density", "ecdf", "sina")) {
    p <- p + ggplot2::scale_color_manual(values = grp_colors, name = group)
  }

  # -- Faceting ----------------------------------------------------------------
  if (!is.null(facet_by)) {
    facet_var <- if (facet_by == "chromosome") "chromosome" else facet_by
    p <- p + ggplot2::facet_wrap(facet_var, scales = "free_y")
  }

  # -- Wilcoxon p-value annotation ----------------------------------------------
  pval_subtitle <- NULL

  if (show_pvalue && n_grp >= 2) {
    # Compute pairwise Wilcoxon tests
    pw <- stats::pairwise.wilcox.test(
      df$value, df$group,
      p.adjust.method = p_adjust_method,
      exact = FALSE
    )

    # Format p-value matrix into label strings
    pmat <- pw$p.value
    pairs <- which(!is.na(pmat), arr.ind = TRUE)
    pval_labels <- vapply(seq_len(nrow(pairs)), function(i) {
      r <- pairs[i, 1]
      cc <- pairs[i, 2]
      pv <- pmat[r, cc]
      sprintf("%s vs %s: p = %s",
              rownames(pmat)[r], colnames(pmat)[cc], format.pval(pv, digits = 3))
    }, character(1))

    pval_text <- paste(pval_labels, collapse = "   |   ")

    if (plot_type %in% c("box", "violin", "sina") && n_grp == 2) {
      # Bracket annotation for two-group box/violin/sina
      pv <- stats::wilcox.test(
        value ~ group, data = df, exact = FALSE
      )$p.value

      y_max <- max(df$value, na.rm = TRUE)
      y_range <- diff(range(df$value, na.rm = TRUE))
      bracket_y <- y_max + y_range * 0.05
      label_y   <- bracket_y + y_range * 0.03

      p <- p +
        ggplot2::annotate(
          "segment",
          x = 1, xend = 2, y = bracket_y, yend = bracket_y,
          linewidth = 0.4
        ) +
        ggplot2::annotate(
          "segment",
          x = 1, xend = 1, y = bracket_y, yend = bracket_y - y_range * 0.01,
          linewidth = 0.4
        ) +
        ggplot2::annotate(
          "segment",
          x = 2, xend = 2, y = bracket_y, yend = bracket_y - y_range * 0.01,
          linewidth = 0.4
        ) +
        ggplot2::annotate(
          "text",
          x = 1.5, y = label_y,
          label = paste0("p = ", format.pval(pv, digits = 3)),
          size = 3.5, hjust = 0.5
        )
    } else if (plot_type %in% c("box", "violin", "sina") && n_grp > 2) {
      # Multiple groups: stacked brackets
      y_max   <- max(df$value, na.rm = TRUE)
      y_range <- diff(range(df$value, na.rm = TRUE))

      for (i in seq_len(nrow(pairs))) {
        r  <- pairs[i, 1]
        cc <- pairs[i, 2]
        pv <- pmat[r, cc]

        # Map group names to x positions (factor levels)
        x1 <- which(grp_levels == colnames(pmat)[cc])
        x2 <- which(grp_levels == rownames(pmat)[r])

        bracket_y <- y_max + y_range * (0.05 + 0.08 * (i - 1))
        label_y   <- bracket_y + y_range * 0.025

        p <- p +
          ggplot2::annotate(
            "segment",
            x = x1, xend = x2, y = bracket_y, yend = bracket_y,
            linewidth = 0.4
          ) +
          ggplot2::annotate(
            "segment",
            x = x1, xend = x1,
            y = bracket_y, yend = bracket_y - y_range * 0.01,
            linewidth = 0.4
          ) +
          ggplot2::annotate(
            "segment",
            x = x2, xend = x2,
            y = bracket_y, yend = bracket_y - y_range * 0.01,
            linewidth = 0.4
          ) +
          ggplot2::annotate(
            "text",
            x = (x1 + x2) / 2, y = label_y,
            label = paste0("p = ", format.pval(pv, digits = 3)),
            size = 3, hjust = 0.5
          )
      }
    } else {
      # density / ecdf / rank / ridgeline: add as subtitle
      pval_subtitle <- pval_text
    }
  }

  # -- Labels and theme --------------------------------------------------------
  legend_pos <- if (show_legend) "right" else "none"

  # Swap axis labels for horizontal plot types
  if (plot_type %in% c("density", "ecdf", "ridgeline")) {
    p <- p + ggplot2::labs(title = title, subtitle = pval_subtitle,
                           x = ylab, y = xlab)
  } else {
    p <- p + ggplot2::labs(title = title, subtitle = pval_subtitle,
                           x = xlab, y = ylab)
  }

  p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle   = ggplot2::element_text(hjust = 0.5, size = 9),
      legend.position = legend_pos,
      strip.text      = ggplot2::element_text(face = "bold")
    )

  p
}
