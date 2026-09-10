# ===========================================================================
# User-supplied nearest-neighbor parameter tables
#
# The contract these tests pin down is narrow but load-bearing: the compiled
# core looks each stack up by key, so a table with a missing key does not
# fail, it silently contributes zero enthalpy and entropy for that stack. Any
# route by which an incomplete table reaches the core is a route to a wrong
# Tm with no error, which is why the resolver rejects rather than repairs.
# ===========================================================================

ref_name <- "DNA_NN_SantaLucia_2004"
ref      <- TmCalculator:::get_table(ref_name)
choices  <- eval(formals(TmCalculator::tm_nn)$nn_table)
resolve  <- function(x) TmCalculator:::.resolve_table(x, "nn_table", choices)

# tm_nn() operates on a GRanges, not on raw strings.
gr_seqs <- TmCalculator::to_genomic_ranges(
  c("ATCGATCGATCGATCGATCG", "GGGCCCAAATTTGGGCCCAA",
    "TTTTTTTTTTAAAAAAAAAA"))

test_that("a built-in name still resolves to the shipped table", {
  r <- resolve(ref_name)
  expect_false(r$user)
  expect_identical(r$name, ref_name)
  expect_identical(r$tbl, ref)
})

test_that("the shipped table, passed back as a matrix, is a no-op", {
  # The strongest available check: feeding a built-in table in through the
  # user path must not move a single digit.
  r <- resolve(ref)
  expect_true(r$user)
  expect_equal(unname(r$tbl), unname(ref[rownames(r$tbl), ]))

  a <- TmCalculator::tm_nn(gr_seqs, nn_table = ref_name)
  b <- TmCalculator::tm_nn(gr_seqs, nn_table = ref)
  expect_equal(a$gr$Tm, b$gr$Tm)
})

test_that("row order is canonicalised, so shuffling has no effect", {
  set.seed(1)
  shuffled <- ref[sample(nrow(ref)), , drop = FALSE]
  expect_identical(rownames(resolve(shuffled)$tbl), rownames(ref))
  expect_equal(TmCalculator::tm_nn(gr_seqs, nn_table = shuffled)$gr$Tm,
               TmCalculator::tm_nn(gr_seqs, nn_table = ref)$gr$Tm)
})

test_that("a missing key is an error rather than a silent zero", {
  # Without this check the calculation proceeds and every sequence containing
  # an AA/TT step is wrong by that stack's contribution.
  incomplete <- ref[setdiff(rownames(ref), "AA/TT"), , drop = FALSE]
  expect_error(resolve(incomplete), "missing 1 key")
})

test_that("extra keys are retained after the canonical block", {
  # This is the modified-base case: a 5mC set adds stacks, it does not
  # replace the canonical ones.
  extra <- rbind(ref, "MG/CG" = c(-9.1, -24.0))
  r <- resolve(extra)
  expect_identical(rownames(r$tbl), c(rownames(ref), "MG/CG"))
  expect_equal(unname(r$tbl["MG/CG", ]), c(-9.1, -24.0))
})

test_that("structural defects are rejected", {
  expect_error(resolve(list(a = 1)), "matrix / data.frame")

  chr <- ref; storage.mode(chr) <- "character"
  expect_error(resolve(chr), "must be numeric")

  one_col <- ref[, 1, drop = FALSE]
  expect_error(resolve(one_col), "at least two columns")

  noname <- ref; rownames(noname) <- NULL
  expect_error(resolve(noname), "row names")

  dup <- rbind(ref, ref[1, , drop = FALSE])
  expect_error(resolve(dup), "duplicated row names")

  nonfinite <- ref; nonfinite["AA/TT", 1] <- NA_real_
  expect_error(resolve(nonfinite), "non-finite")
})

test_that("a transposed row warns through the reverse-complement check", {
  # A key and its character reversal are the same duplex read from opposite
  # strands. Four rows of the shipped table were once transposed this way.
  bad <- ref
  bad["AC/TG", ] <- bad["AC/TG", ] + c(1, 1)
  expect_warning(resolve(bad), "reverse complement")
})

test_that("the salt_mM attribute survives and still suppresses correction", {
  tbl <- ref
  attr(tbl, "salt_mM") <- 500
  r <- resolve(tbl)
  expect_equal(as.numeric(attr(r$tbl, "salt_mM")), 500)

  # At the fitted concentration the correction is skipped silently; away from
  # it the user is warned that a correction is being stacked on a
  # condition-specific set.
  at_fit <- TmCalculator::tm_nn(gr_seqs, nn_table = tbl, Na = 500)
  expect_false(at_fit$options[["Salt correction applied"]])
  expect_warning(TmCalculator::tm_nn(gr_seqs, nn_table = tbl, Na = 50),
                 "fitted at")
})

test_that("provenance records that the table was user-supplied", {
  out <- TmCalculator::tm_nn(gr_seqs, nn_table = ref)
  expect_match(out$options[["Thermodynamic NN values"]], "user-supplied")
  expect_match(out$options[["Thermodynamic NN values"]], ref_name, fixed = TRUE)
})
