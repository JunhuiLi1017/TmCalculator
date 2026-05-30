#' Compare numeric GRanges metadata across groups
#'
#' Test whether one or more numeric metadata columns (e.g. \code{Tm}, \code{GC})
#' differ between levels of a categorical grouping column in a \code{GRanges}
#' object.
#'
#' @details
#' The test used depends on \code{method} and the number of group levels:
#'
#' \tabular{lll}{
#'   Groups \tab \code{method = "wilcoxon"} \tab \code{method = "t.test"} \cr
#'   2 \tab Wilcoxon rank-sum (or signed-rank if \code{paired = TRUE}) \tab
#'       Welch two-sample \eqn{t}-test (or paired \eqn{t}-test) \cr
#'   \eqn{\ge} 3 \tab Kruskal-Wallis \tab One-way ANOVA (Welch) \cr
#' }
#'
#' When there are three or more groups and \code{posthoc = TRUE}, pairwise
#' follow-up tests are run (Wilcoxon or Welch \eqn{t}-test) with \code{p.adjust}
#' multiple-testing correction.
#'
#' @param gr A \code{GRanges} object.
#' @param target Character vector of metadata column names to test (e.g.
#'   \code{c("Tm", "GC")}). All must be numeric.
#' @param method Statistical test family. One of:
#'   \itemize{
#'     \item \code{"wilcoxon"}: rank-based tests (Wilcoxon / Kruskal-Wallis).
#'     \item \code{"t.test"}: Gaussian-based tests (\eqn{t}-test / ANOVA).
#'   }
#' @param group Character. Name of a metadata column in \code{gr} defining groups.
#' @param min_n_per_group Minimum number of non-missing observations required
#'   per group. Default: \code{2}.
#' @param paired Logical. Only for exactly two groups. Default: \code{FALSE}.
#' @param alternative Character. Alternative hypothesis for two-group tests.
#'   One of \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#' @param posthoc Logical. For \eqn{\ge} 3 groups, run pairwise follow-up tests
#'   with multiple-testing correction. Default: \code{TRUE}.
#' @param p.adjust.method Multiple-testing correction for post-hoc pairwise
#'   comparisons. Passed to \code{\link[stats]{pairwise.wilcox.test}} /
#'   \code{\link[stats]{pairwise.t.test}}. One of:
#'   \code{"holm"} (default), \code{"hochberg"}, \code{"hommel"},
#'   \code{"bonferroni"}, \code{"BH"} (Benjamini-Hochberg FDR),
#'   \code{"BY"} (Benjamini-Yekutieli), \code{"fdr"} (alias of \code{"BH"}),
#'   or \code{"none"}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{results}}{Data frame with one row per \code{target}: omnibus
#'       test name, statistic, \code{p.value}, and group information.}
#'     \item{\code{summary}}{Per-group descriptive statistics
#'       (\code{n}, \code{mean}, \code{sd}, \code{median}).}
#'     \item{\code{pairwise}}{Pairwise comparisons when \code{posthoc = TRUE} and
#'       there are \eqn{\ge} 3 groups; otherwise \code{NULL}.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' set.seed(1)
#' gr <- GRanges(
#'   seqnames = rep("chr1", 90),
#'   ranges = IRanges(start = sort(sample(1e6, 90)), width = 50),
#'   Tm = c(rnorm(30, 74, 2), rnorm(30, 70, 2), rnorm(30, 66, 2)),
#'   GC = runif(90, 40, 60)
#' )
#' gr$group <- rep(c("high", "mid", "low"), each = 30)
#'
#' # Two groups
#' gr2 <- gr[gr$group %in% c("high", "low")]
#' compare_groups(gr2, target = "Tm", group = "group", posthoc = FALSE)
#'
#' # Three or more groups (Kruskal-Wallis + pairwise Wilcoxon)
#' compare_groups(gr, target = c("Tm", "GC"), group = "group",
#'                      method = "wilcoxon")
#'
#' # Three or more groups (one-way ANOVA + pairwise t-tests)
#' compare_groups(gr, target = "Tm", group = "group", method = "t.test")
#'
#' # Post-hoc with Benjamini-Hochberg FDR control
#' compare_groups(gr, target = "Tm", group = "group",
#'                      p.adjust.method = "BH")
#' }
#'
#' @importFrom GenomicRanges mcols
#' @importFrom stats wilcox.test kruskal.test t.test oneway.test pairwise.wilcox.test
#'   pairwise.t.test p.adjust
#' @encoding UTF-8
#' @author Junhui Li
#' @export

compare_groups <- function(
    gr,
    target = "Tm",
    method = c("wilcoxon", "t.test"),
    group,
    min_n_per_group = 2L,
    paired = FALSE,
    alternative = c("two.sided", "less", "greater"),
    posthoc = TRUE,
    p.adjust.method = c(
      "holm", "hochberg", "hommel", "bonferroni",
      "BH", "BY", "fdr", "none"
    )
) {
  if (!inherits(gr, "GRanges")) {
    stop("'gr' must be a GRanges object.")
  }

  method <- match.arg(method)
  alternative <- match.arg(alternative)
  p.adjust.method <- match.arg(p.adjust.method)
  meta_cols <- names(GenomicRanges::mcols(gr))

  if (missing(group) || length(group) != 1L || !group %in% meta_cols) {
    stop(
      "'group' must be a single metadata column name in gr.\n",
      "Available columns: ", paste(meta_cols, collapse = ", ")
    )
  }

  target <- unique(as.character(target))
  if (length(target) == 0L) {
    stop("'target' must name at least one metadata column.")
  }
  missing_tgt <- setdiff(target, meta_cols)
  if (length(missing_tgt) > 0L) {
    stop(
      "Column(s) not found in gr: ", paste(missing_tgt, collapse = ", "),
      "\nAvailable columns: ", paste(meta_cols, collapse = ", ")
    )
  }

  grp_raw <- GenomicRanges::mcols(gr)[[group]]
  if (all(is.na(grp_raw))) {
    stop("All values in group column '", group, "' are NA.")
  }
  grp <- droplevels(factor(grp_raw))
  n_groups <- nlevels(grp)
  if (n_groups < 2L) {
    stop(
      "Group column '", group, "' must have at least two non-empty levels ",
      "(found: ", n_groups, ")."
    )
  }
  if (paired && n_groups >= 3L) {
    stop("'paired = TRUE' is only supported when there are exactly two groups.")
  }
  run_posthoc <- posthoc && n_groups >= 3L

  .group_counts <- function(g) {
    table(droplevels(g))
  }

  .validate_counts <- function(counts) {
    if (any(counts < min_n_per_group)) {
      stop(
        "Each group needs at least ", min_n_per_group,
        " non-missing value(s). Counts: ",
        paste(names(counts), counts, sep = "=", collapse = ", ")
      )
    }
  }

  .run_omnibus <- function(x, g) {
    k <- nlevels(g)

    if (method == "wilcoxon") {
      if (k == 2L) {
        if (paired) {
          split_x <- split(x, g)
          fit <- stats::wilcox.test(
            split_x[[1]], split_x[[2]],
            paired = TRUE,
            alternative = alternative,
            exact = FALSE
          )
          test_name <- "Wilcoxon signed-rank (paired)"
        } else {
          fit <- stats::wilcox.test(
            x ~ g,
            alternative = alternative,
            exact = FALSE
          )
          test_name <- "Wilcoxon rank-sum"
        }
      } else {
        fit <- stats::kruskal.test(x ~ g)
        test_name <- "Kruskal-Wallis"
      }
    } else if (k == 2L) {
      if (paired) {
        split_x <- split(x, g)
        fit <- stats::t.test(
          split_x[[1]], split_x[[2]],
          paired = TRUE,
          alternative = alternative
        )
        test_name <- "Paired t-test"
      } else {
        fit <- stats::t.test(
          x ~ g,
          alternative = alternative,
          var.equal = FALSE
        )
        test_name <- "Welch two-sample t-test"
      }
    } else {
      fit <- stats::oneway.test(x ~ g, var.equal = FALSE)
      test_name <- "One-way ANOVA (Welch)"
    }

    list(
      test = test_name,
      statistic = unname(fit$statistic[1]),
      p.value = fit$p.value,
      parameter = if (!is.null(fit$parameter)) unname(fit$parameter[1]) else NA_real_
    )
  }

  .run_posthoc <- function(x, g) {
    if (method == "wilcoxon") {
      pw <- stats::pairwise.wilcox.test(
        x, g,
        p.adjust.method = p.adjust.method,
        exact = FALSE
      )
      test_name <- "Pairwise Wilcoxon"
    } else {
      pw <- stats::pairwise.t.test(
        x, g,
        p.adjust.method = p.adjust.method,
        pool.sd = FALSE
      )
      test_name <- "Pairwise Welch t-test"
    }

    pmat <- pw$p.value
    if (all(is.na(pmat))) {
      return(NULL)
    }

    pairs <- which(!is.na(pmat), arr.ind = TRUE)
    data.frame(
      group1 = rownames(pmat)[pairs[, 1]],
      group2 = colnames(pmat)[pairs[, 2]],
      p.value = pmat[pairs],
      p.adjust.method = p.adjust.method,
      test = test_name,
      stringsAsFactors = FALSE
    )
  }

  .target_summary <- function(x, g) {
    ok <- !is.na(x) & !is.na(g)
    x <- as.numeric(x[ok])
    g <- droplevels(factor(g[ok]))
    do.call(
      rbind,
      lapply(levels(g), function(lv) {
        xi <- x[g == lv]
        data.frame(
          group = lv,
          n = length(xi),
          mean = mean(xi),
          sd = stats::sd(xi),
          median = stats::median(xi),
          stringsAsFactors = FALSE
        )
      })
    )
  }

  result_rows <- vector("list", length(target))
  summary_rows <- vector("list", length(target))
  pairwise_rows <- vector("list", length(target))

  for (i in seq_along(target)) {
    col <- target[i]
    x <- GenomicRanges::mcols(gr)[[col]]
    if (!is.numeric(x)) {
      stop("Column '", col, "' must be numeric.")
    }

    ok <- !is.na(x) & !is.na(grp)
    x_use <- as.numeric(x[ok])
    g_use <- droplevels(factor(grp[ok]))
    counts <- .group_counts(g_use)

    if (length(x_use) == 0L) {
      stop("No non-missing observations for target '", col, "'.")
    }
    .validate_counts(counts)

    tst <- .run_omnibus(x_use, g_use)

    result_rows[[i]] <- data.frame(
      target = col,
      group = group,
      method = method,
      test = tst$test,
      statistic = tst$statistic,
      df = tst$parameter,
      p.value = tst$p.value,
      n_groups = nlevels(g_use),
      n_total = length(x_use),
      group_levels = paste(names(counts), collapse = " | "),
      group_n = paste(names(counts), counts, sep = "=", collapse = ", "),
      stringsAsFactors = FALSE
    )

    sm <- .target_summary(x, grp)
    sm$target <- col
    summary_rows[[i]] <- sm

    if (run_posthoc) {
      pw <- .run_posthoc(x_use, g_use)
      if (!is.null(pw)) {
        pw$target <- col
        pairwise_rows[[i]] <- pw
      }
    }
  }

  pairwise_out <- if (length(pairwise_rows) > 0L && any(!vapply(pairwise_rows, is.null, logical(1)))) {
    do.call(rbind, pairwise_rows[!vapply(pairwise_rows, is.null, logical(1))])
  } else {
    NULL
  }

  list(
    results = do.call(rbind, result_rows),
    summary = do.call(rbind, summary_rows),
    pairwise = pairwise_out
  )
}
