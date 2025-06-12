Sys.setenv(R_TESTS="")
test_that("test_tm_calculate failed", {
  input_seq <- c("ATGCGATGCGAAGGCGATGGCGTGTAGAATAGATCACATACTGCATAGCTGATC", "ATGCGATGCGCCCGGAGATAGAAGGCGTAGATACAGATCAGTAGCACCTTGAGAC")
  result <- tm_calculate(input_seq)
  expect_equal(round(result$tm$Tm$Tm,5), c(67.06562,69.64434))
  
  fasta_file  <- system.file("extdata","BSgenome.Hsapiens.UCSC.hg38.fasta", package = "TmCalculator")
  gr_tm <- tm_calculate(fasta_file)
  expect_equal(round(gr_tm$tm$Tm$Tm,5), c(67.06562,69.64434))
})
