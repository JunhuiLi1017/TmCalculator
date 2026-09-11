# ===========================================================================
# Length-of-validity guard on the Wallace rule
#
# The 2 + 4 rule has no length term, so its output is unbounded in sequence
# length: a 200 bp window comes back at roughly 600 degrees C. Nothing in the
# arithmetic signals that the number is meaningless, which is why the guard
# has to be explicit. These tests pin the boundary and, just as importantly,
# pin that the returned values themselves are untouched.
# ===========================================================================

mk <- function(n, base = "AT") {
  TmCalculator::to_genomic_ranges(
    substr(strrep(base, ceiling(n / nchar(base))), 1L, n))
}

test_that("a sequence inside the calibrated range is silent", {
  expect_silent(TmCalculator::tm_wallace(mk(20L)))
})

test_that("the boundary is exclusive at 30 nt", {
  # 30 passes, 31 warns. Stated as two assertions rather than one so that a
  # future change to `limit` fails here rather than drifting unnoticed.
  expect_silent(TmCalculator::tm_wallace(mk(30L)))
  expect_warning(TmCalculator::tm_wallace(mk(31L)), "exceed 30 nt")
})

test_that("the warning names the alternative the reviewer asked for", {
  expect_warning(TmCalculator::tm_wallace(mk(200L)), "tm_nn\\(\\)")
})

test_that("the warning counts offenders rather than listing them", {
  # Mixed input: two long, one short. The message must report 2 of 3 and the
  # longest, because the genome-wide case has millions of windows.
  gr <- TmCalculator::to_genomic_ranges(
    c(strrep("AT", 10L), strrep("AT", 30L), strrep("AT", 50L)))
  expect_warning(TmCalculator::tm_wallace(gr),
                 "2 of 3 sequences exceed 30 nt \\(longest 100 nt\\)")
})

test_that("the guard warns but does not alter a single value", {
  # The rule stays available to users who knowingly apply it out of range.
  gr <- mk(200L)
  out <- suppressWarnings(TmCalculator::tm_wallace(gr))
  # 200 bp of alternating AT: 0 GC, so 2 * 200.
  expect_equal(out$gr$Tm, 400)
  expect_equal(out$gr$GC, 0)
})

test_that("the warning fires once per call, not once per chunk", {
  # A per-chunk check would emit one warning per block and would be lost or
  # reordered on a parallel backend.
  gr <- TmCalculator::to_genomic_ranges(rep(strrep("AT", 40L), 500L))
  w <- testthat::capture_warnings(TmCalculator::tm_wallace(gr))
  expect_length(w, 1L)
})
