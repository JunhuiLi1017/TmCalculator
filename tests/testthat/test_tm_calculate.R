Sys.setenv(R_TESTS="")
test_that("test_tm_calculate failed", {
  input_seq <- c("ATGCGATGCGAAGGCGATGGCGTGTAGAATAGATCACATACTGCATAGCTGATC", "ATGCGATGCGCCCGGAGATAGAAGGCGTAGATACAGATCAGTAGCACCTTGAGAC")
  result <- tm_calculate(input_seq)
  expect_equal(round(result$gr$Tm, 5), c(67.06562, 69.64434), tolerance = 0.1)
  
  fasta_file  <- system.file("extdata","BSgenome.Hsapiens.UCSC.hg38.fasta", package = "TmCalculator")
  gr_tm <- tm_calculate(fasta_file)
  expect_equal(
    round(gr_tm$gr$Tm, 5)[c(1, 10, 50, 100, 141)],
    c(76.57807, 72.88846, 70.04297, 67.24338, 75.47229),
    tolerance = 0.1
  )
})

test_that("BPPARAM parallel results match serial", {
  skip_if_not_installed("BiocParallel")

  set.seed(20260831)
  input_seq <- vapply(seq_len(9), function(i) {
    paste(sample(c("A", "C", "G", "T"), 60, replace = TRUE), collapse = "")
  }, character(1))

  # Explicit SerialParam is identical to the default
  serial <- tm_calculate(input_seq, method = "tm_nn")
  serial2 <- tm_calculate(input_seq, method = "tm_nn",
                          BPPARAM = BiocParallel::SerialParam())
  expect_equal(serial$gr$Tm, serial2$gr$Tm)
  expect_equal(serial$gr$GC, serial2$gr$GC)

  # Multi-worker backend (SnowParam works on all platforms) is identical too,
  # including order, for all three methods
  bp <- BiocParallel::SnowParam(workers = 2)

  # PSOCK workers load the INSTALLED package, not the load_all() dev version.
  # Skip (rather than fail) when the installed copy predates the Rcpp core,
  # e.g. during devtools::test() without a prior devtools::install().
  # R CMD check always installs first, so the test runs there.
  workers_ok <- tryCatch({
    r <- BiocParallel::bplapply(1, function(i) {
      exists("cpp_tm_nn_dhds", envir = asNamespace("TmCalculator"))
    }, BPPARAM = bp)
    isTRUE(r[[1]])
  }, error = function(e) FALSE)
  skip_if_not(
    workers_ok,
    "installed TmCalculator lacks the Rcpp core; run devtools::install() first"
  )

  for (m in c("tm_nn", "tm_gc", "tm_wallace")) {
    res_serial <- tm_calculate(input_seq, method = m)
    res_par <- tm_calculate(input_seq, method = m, BPPARAM = bp)
    expect_equal(res_serial$gr$Tm, res_par$gr$Tm)
    expect_equal(res_serial$gr$GC, res_par$gr$GC)
  }
})
