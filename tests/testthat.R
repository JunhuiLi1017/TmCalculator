if (!require(testthat, quietly = TRUE)) {
  stop("testthat is not installed")
}
if (!require(TmCalculator, quietly = TRUE)) {
  stop("TmCalculator is not installed")
}

test_check("TmCalculator")