test_that("test_granges_targets compares two groups", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1:20, width = 10),
    Tm = c(rnorm(10, 75, 1), rnorm(10, 65, 1)),
    group = rep(c("A", "B"), each = 10)
  )

  out <- test_granges_targets(
    gr, target = "Tm", group = "group",
    method = "wilcoxon", posthoc = FALSE
  )
  expect_type(out, "list")
  expect_equal(nrow(out$results), 1L)
  expect_equal(out$results$test, "Wilcoxon rank-sum")
  expect_true(out$results$p.value < 0.05)
  expect_null(out$pairwise)
})

test_that("test_granges_targets uses Kruskal-Wallis for three groups", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1:30, width = 10),
    Tm = c(rnorm(10, 80, 1), rnorm(10, 75, 1), rnorm(10, 70, 1)),
    group = rep(c("high", "mid", "low"), each = 10)
  )

  out <- test_granges_targets(
    gr, target = "Tm", group = "group",
    method = "wilcoxon", posthoc = TRUE
  )
  expect_equal(out$results$test, "Kruskal-Wallis")
  expect_equal(out$results$n_groups, 3L)
  expect_false(is.null(out$pairwise))
  expect_equal(out$pairwise$test[1], "Pairwise Wilcoxon")
  expect_equal(nrow(out$pairwise), 3L)
})

test_that("test_granges_targets uses ANOVA for three groups with t.test", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1:30, width = 10),
    Tm = c(rnorm(10, 80, 1), rnorm(10, 75, 1), rnorm(10, 70, 1)),
    group = rep(c("high", "mid", "low"), each = 10)
  )

  out <- test_granges_targets(
    gr, target = "Tm", group = "group",
    method = "t.test", posthoc = TRUE
  )
  expect_equal(out$results$test, "One-way ANOVA (Welch)")
  expect_equal(out$pairwise$test[1], "Pairwise Welch t-test")
})

test_that("test_granges_targets runs for multiple targets", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1:12, width = 10),
    Tm = c(rnorm(6, 80, 1), rnorm(6, 70, 1)),
    GC = c(rnorm(6, 50, 1), rnorm(6, 40, 1)),
    group = rep(c("mut", "wt"), each = 6)
  )

  out <- test_granges_targets(
    gr, target = c("Tm", "GC"), group = "group",
    method = "t.test", posthoc = FALSE
  )
  expect_equal(nrow(out$results), 2L)
  expect_equal(sort(out$results$target), c("GC", "Tm"))
})
