# ===========================================================================
# Overlap-length weighting in integrate_granges()
#
# The reference case is a coverage track, because it has a defined truth: the
# mean per-base depth over a window is a single number, and any aggregation
# that claims to summarise the window can be checked against it.
#
# Layout, all on one contig:
#
#   window 1   [  1 .. 200]      window 2   [401 .. 600]
#   depth      [  1 .. 190] = 5
#              [191 .. 400] = 200
#                                           [401 .. 450] = 60
#
# Window 1 straddles a sharp transition; window 2 is three quarters empty.
# ===========================================================================

suppressPackageStartupMessages(library(GenomicRanges))

win <- GenomicRanges::GRanges(
  "chr1", IRanges::IRanges(start = c(1, 401), width = 200), Tm = c(70, 72))
cov <- GenomicRanges::GRanges(
  "chr1", IRanges::IRanges(start = c(1, 191, 401), end = c(190, 400, 450)),
  cov = c(5, 200, 60))

run <- function(...) TmCalculator::integrate_granges(win, cov, ...)

test_that("the default is unchanged, to the digit", {
  # Weighting is opt-in precisely so that existing analyses do not move.
  expect_equal(run(strategy = "overlap")$cov, c(102.5, 60))
})

test_that("unweighted aggregation misstates a window that straddles a step", {
  # 5 and 200 enter with equal weight although one covers 190 bp and the
  # other 10 bp, which is a sevenfold overstatement of the true 14.75.
  expect_equal(run(strategy = "overlap")$cov[1], 102.5)
})

test_that("overlap weighting recovers the true mean depth", {
  # (5 * 190 + 200 * 10) / 200
  expect_equal(run(strategy = "overlap", weight = "overlap")$cov[1], 14.75)
})

test_that("min_overlap cannot fix it, it only reverses the bias", {
  # Excluding the 10 bp feature leaves 5 against a true 14.75: a threefold
  # understatement where the default gave a sevenfold overstatement.
  expect_equal(run(strategy = "overlap", min_overlap = 20L)$cov[1], 5)
})

test_that("weighting alone does not describe a partly covered window", {
  # Window 2 holds one feature, so weighting changes nothing; 60 is the mean
  # over the covered 50 bp, not over the 200 bp window, and only the coverage
  # column distinguishes the two readings.
  w <- run(strategy = "overlap", weight = "overlap")
  expect_equal(w$cov[2], 60)

  cf <- run(strategy = "overlap", weight = "overlap",
            report_coverage = TRUE)$covered_frac
  expect_equal(cf, c(1, 0.25))
  expect_equal(w$cov[2] * cf[2], 15)          # mean over the whole window
})

test_that("covered_frac never exceeds one when features overlap", {
  # Without reducing the features first this would be 1.5 and would not be a
  # fraction at all.
  dup <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(1, 1, 101), end = c(100, 100, 200)),
    cov = c(1, 1, 1))
  cf <- TmCalculator::integrate_granges(win[1], dup, strategy = "overlap",
                                        report_coverage = TRUE)$covered_frac
  expect_equal(cf, 1)
})

test_that("a weighted mean is refused rather than faked", {
  expect_error(run(strategy = "overlap", weight = "overlap", agg_fun = max),
               "cannot honour a different")
  expect_error(run(strategy = "nearest", weight = "overlap"),
               "does not apply to strategy")
})

test_that("character columns are joined under either weighting", {
  lab <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(1, 191), end = c(190, 400)),
    cov = c(5, 200), name = c("a", "b"))
  a <- TmCalculator::integrate_granges(win[1], lab, strategy = "overlap")
  b <- TmCalculator::integrate_granges(win[1], lab, strategy = "overlap",
                                       weight = "overlap")
  expect_equal(a$name, "a,b")
  expect_equal(b$name, "a,b")
  expect_equal(b$cov, 14.75)
})

test_that("weighting is available in the window and bin strategies", {
  expect_silent(run(strategy = "window", weight = "overlap",
                    window_size = 50L))
  expect_s4_class(
    TmCalculator::integrate_granges(win, cov, strategy = "bin",
                                    weight = "overlap", bin_size = 1000),
    "GRanges")
})
