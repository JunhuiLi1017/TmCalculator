test_that(".gc_vec reproduces gc() on unambiguous sequences", {
  set.seed(11)
  seqs <- vapply(seq_len(200), function(i)
    paste0(sample(c("A", "C", "G", "T"), sample(8:300, 1), TRUE), collapse = ""),
    character(1))
  expect_equal(TmCalculator:::.gc_vec(seqs),
               vapply(seqs, function(s) TmCalculator::gc(s), numeric(1),
                      USE.NAMES = FALSE))
})

test_that(".gc_vec reproduces gc() with ambiguous codes apportioned", {
  set.seed(12)
  alph <- c("A", "C", "G", "T", "S", "W", "M", "K", "R", "Y",
            "V", "H", "D", "B", "N")
  seqs <- vapply(seq_len(300), function(i)
    paste0(sample(alph, sample(8:200, 1), TRUE), collapse = ""), character(1))
  for (amb in c(FALSE, TRUE)) {
    expect_equal(
      TmCalculator:::.gc_vec(seqs, ambiguous = amb),
      vapply(seqs, function(s) TmCalculator::gc(s, ambiguous = amb),
             numeric(1), USE.NAMES = FALSE),
      info = paste("ambiguous =", amb))
  }
})

test_that(".gc_vec handles the degenerate cases gc() special-cases", {
  # A sequence with no countable base gives NA in both implementations; an
  # all-N sequence is the case that reaches the ngc + nat == 0 branch.
  expect_true(is.na(TmCalculator:::.gc_vec("NNNNNNNN")))
  expect_true(is.na(TmCalculator::gc("NNNNNNNN")))

  expect_true(is.na(TmCalculator:::.gc_vec(NA_character_)))

  # Mixed vector: NA must not shift the positions of the other results.
  s <- c("GGCC", NA_character_, "ATAT")
  out <- TmCalculator:::.gc_vec(s)
  expect_equal(out[1], 100)
  expect_true(is.na(out[2]))
  expect_equal(out[3], 0)

  # Inosine is in gc()'s accepted alphabet but is rejected by DNAStringSet,
  # which is one reason the counting is done in cpp_base_counts() rather than
  # through Biostrings::letterFrequency().
  expect_silent(TmCalculator:::.gc_vec("AICGT"))
})

test_that("tm_gc is unchanged by the vectorised chunk worker", {
  skip_if_not_installed("BiocParallel")
  set.seed(13)
  seqs <- vapply(seq_len(150), function(i)
    paste0(sample(c("A", "C", "G", "T"), 200, TRUE), collapse = ""), character(1))
  gr <- TmCalculator::to_genomic_ranges(seqs)

  # Reference: the calculation as it was performed per sequence, with the
  # coefficients read from the same table rather than hard-coded, so the test
  # cannot drift from the shipped parameters.
  variant <- "Schildkraut1965"
  co  <- TmCalculator:::get_table("GC_VARTAB")[variant, ]
  A   <- as.numeric(co$A); B <- as.numeric(co$B); Cc <- as.numeric(co$C)
  sm  <- co$salt_correct

  ref <- vapply(seqs, function(s) {
    tm <- A + B * TmCalculator::gc(s) - Cc / nchar(s)
    if (!is.na(sm)) {
      tm <- tm + TmCalculator::salt_correct(
        Na = 50, K = 0, Tris = 0, Mg = 0, dNTPs = 0,
        method = sm, input_seq = s)
    }
    tm
  }, numeric(1), USE.NAMES = FALSE)

  got <- suppressWarnings(
    TmCalculator::tm_calculate(gr, method = "tm_gc", variant = variant, Na = 50))
  expect_equal(as.numeric(got$gr$Tm), ref, tolerance = 1e-8)
})

test_that("every GC path in the package agrees on one definition", {
  # gc(), tm_gc(), tm_wallace() and tm_nn() must return the same GC for the
  # same sequence. Before unification tm_nn used (G+C)/length while the others
  # used (G+C)/(A+C+G+T), and coor_to_genomic_ranges() wrote a 0-1 fraction.
  seqs <- c("ACGTACGTAC", "GGGGCCCCAA", "ATATATATAT")
  gr   <- TmCalculator::to_genomic_ranges(seqs)

  ref <- vapply(seqs, function(s) TmCalculator::gc(s), numeric(1),
                USE.NAMES = FALSE)

  for (m in c("tm_nn", "tm_gc", "tm_wallace")) {
    got <- suppressWarnings(TmCalculator::tm_calculate(gr, method = m))
    expect_equal(as.numeric(got$gr$GC), ref, info = m)
  }
})

test_that("inosine is excluded from the GC denominator", {
  # I is an accepted but undeterminable base, so it must enter neither the
  # numerator nor the denominator. This is the one case where the removed
  # .GC_fast() definition ((G+C)/length) diverged from gc()
  # ((G+C)/(A+C+G+T)): it would have returned 100 * 4/6 here.
  expect_equal(TmCalculator::gc("GGCCII"), 100)

  # Same four G/C bases, but with two determinable A's instead of the two I's,
  # which do enter the denominator.
  expect_equal(TmCalculator::gc("GGCCAA"), 100 * 4 / 6)

  # N behaves the same way as I.
  expect_equal(TmCalculator::gc("GGCCNN"), 100)
})

test_that("gc() accepts a split sequence and a whole sequence alike", {
  expect_equal(TmCalculator::gc(c("a", "t", "c", "g")),
               TmCalculator::gc("atcg"))
})

test_that("cpp_base_counts exposes N and I so the denominator is inspectable", {
  m <- TmCalculator:::cpp_base_counts(c("ACGTNNII", NA_character_))
  expect_equal(unname(m[1, "N"]), 2L)
  expect_equal(unname(m[1, "I"]), 2L)
  expect_equal(unname(m[1, "len"]), 8L)
  expect_true(is.na(m[2, "len"]))
})

test_that("FASTA reading is unchanged after dropping seqinr", {
  # example1.fasta mixes a DNA sequence with an RNA one (U), and carries
  # header lines of two shapes. readBStringSet() must accept both, preserve
  # case, and name sequences by the first word of the header, exactly as
  # seqinr::read.fasta(as.string = TRUE, forceDNAtolower = FALSE) did.
  f  <- system.file("extdata", "example1.fasta", package = "TmCalculator")
  skip_if(f == "")
  gr <- TmCalculator::fa_to_genomic_ranges(f)

  expect_gt(length(gr), 0)
  seqs <- as.character(gr$sequence)
  expect_true(any(grepl("U", seqs, fixed = TRUE)))   # RNA survived
  expect_false(any(is.na(seqs)))
})
