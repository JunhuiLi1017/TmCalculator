test_that("sequence box parses plain and FASTA input", {
  plain <- .tcs_parse_sequence_box("ATCG\nGGCC")
  expect_length(plain, 2L)
  expect_equal(unname(plain[1]), "ATCG")

  fa <- .tcs_parse_sequence_box(">p1\nATCG\nTTAA\n>p2\nGGCC")
  expect_length(fa, 2L)
  expect_equal(names(fa), c("p1", "p2"))
  expect_equal(unname(fa["p1"]), "ATCGTTAA")
})

test_that("NULL-coalescing operator works", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
  expect_equal(character(0) %||% "x", "x")
})

test_that("app UI assembles without error", {
  skip_if_not_installed("bslib")
  skip_if_not_installed("shiny")
  ui <- app_ui()
  expect_s3_class(ui, "shiny.tag.list")
})

test_that("end-to-end Tm calculation on pasted sequences", {
  skip_if_not_installed("TmCalculator")
  out <- TmCalculator::tm_calculate(
    c("ATCGTGCGTAGCAGTACGATCAGTAG", "GGGCCCATATATGCGC"),
    method = "tm_nn")
  expect_true(!is.null(out$gr))
  expect_true("Tm" %in% names(GenomicRanges::mcols(out$gr)))
})

test_that("tm_gc honours custom A/B/C/D coefficients", {
  skip_if_not_installed("TmCalculator")
  custom <- TmCalculator::tm_calculate(
    "ATCGTGCGTAGCAGTACGATCAGTAG", method = "tm_gc",
    userset = c(81.5, 0.41, 600, 1))
  expect_true("Tm" %in% names(GenomicRanges::mcols(custom$gr)))
})
