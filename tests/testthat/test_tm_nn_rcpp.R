Sys.setenv(R_TESTS = "")

# Verify the Rcpp-backed chunk worker reproduces the pure-R reference
# implementation exactly across parameter sets, corrections and edge cases.

.compare_chunk <- function(seqs, cseqs, ...) {
  defaults <- list(
    ambiguous = FALSE, shift = 0,
    nn_tbl = TmCalculator:::get_table("DNA_NN_SantaLucia_2004"),
    tmm_tbl = TmCalculator:::get_table("DNA_TMM_Bommarito_2000"),
    imm_tbl = TmCalculator:::get_table("DNA_IMM_Peyret_1999"),
    de_tbl = TmCalculator:::get_table("DNA_DE_Bommarito_2000"),
    end_tbl = matrix(numeric(0), nrow = 0, ncol = 2),
    dnac_high = 25, dnac_low = 25, self_comp = FALSE,
    Na = 50, K = 0, Tris = 0, Mg = 0, dNTPs = 0,
    salt_fn = "Schildkraut2010",
    DMSO = 0, formamide_unit = list(value = 0, unit = "percent"),
    dmso_factor = 0.75, formamide_factor = 0.65
  )
  args <- utils::modifyList(defaults, list(...))
  chunk <- list(sequence = seqs, complement = cseqs)
  cpp <- do.call(TmCalculator:::.tm_nn_chunk, c(list(chunk), args))
  ref <- do.call(TmCalculator:::.tm_nn_chunk_r, c(list(chunk), args))
  expect_equal(cpp$Tm, ref$Tm, tolerance = 1e-10)
  expect_equal(cpp$GC, ref$GC, tolerance = 1e-10)
}

.rand_seqs <- function(n, len, seed) {
  set.seed(seed)
  vapply(seq_len(n), function(i) {
    paste(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
  }, character(1))
}

test_that("Rcpp core matches R core on perfect duplexes", {
  seqs <- .rand_seqs(50, 40, 1)
  cseqs <- complement_fast(seqs)
  .compare_chunk(seqs, cseqs)
  .compare_chunk(seqs, cseqs, salt_fn = "Owczarzy2004")
  .compare_chunk(seqs, cseqs, salt_fn = "Owczarzy2008", Mg = 1.5)
  .compare_chunk(seqs, cseqs, salt_fn = "Owczarzy2008", Mg = 10, Na = 5)
  .compare_chunk(seqs, cseqs, salt_fn = "none")
  .compare_chunk(seqs, cseqs, salt_fn = "Wetmur1991", K = 20, Tris = 10)
})

test_that("Rcpp core matches R core with internal mismatches", {
  seqs <- .rand_seqs(30, 30, 2)
  cseqs <- complement_fast(seqs)
  # introduce a mismatch in the middle of each complement
  substr(cseqs, 15, 15) <- vapply(substr(cseqs, 15, 15), function(b) {
    sample(setdiff(c("A", "C", "G", "T"), b), 1)
  }, character(1))
  .compare_chunk(seqs, cseqs)
})

test_that("Rcpp core matches R core with shift / dangling ends", {
  seqs <- .rand_seqs(20, 25, 3)
  cseqs <- complement_fast(.rand_seqs(20, 27, 4))
  for (sh in c(-2, -1, 0, 1, 2)) {
    .compare_chunk(seqs, cseqs, shift = sh)
  }
})

test_that("Rcpp core matches R core with self-complementary and chem corrections", {
  s <- "GCGCATGCGC"  # palindromic
  cs <- complement_fast(s)
  .compare_chunk(s, cs, self_comp = TRUE)
  seqs <- .rand_seqs(10, 35, 5)
  cseqs <- complement_fast(seqs)
  .compare_chunk(seqs, cseqs, DMSO = 5)
  .compare_chunk(seqs, cseqs,
                 formamide_unit = list(value = 1.25, unit = "molar"))
  .compare_chunk(seqs, cseqs,
                 formamide_unit = list(value = 10, unit = "percent"),
                 formamide_factor = 0.72)
})

test_that("Rcpp core matches R core with inosine and short sequences", {
  seqs <- c("ACGTIACGT", "AI", "ACGT")
  cseqs <- c(complement_fast("ACGTCACGT"), "TC", complement_fast("ACGT"))
  .compare_chunk(seqs, cseqs)
})

test_that("Rcpp core cleans lowercase/ambiguous characters like the R path", {
  # cleaning (uppercase + strip non-ACGTI) now happens inside the workers
  seqs <- c("acgtACGTacgt", "ACGTRYSWacgt", "nnACGTACGTnn")
  cseqs <- c("tgcaTGCAtgca", "TGCAYRSWtgca", "nnTGCATGCAnn")
  .compare_chunk(seqs, cseqs)
  # cleaned-to-short sequences yield NA instead of aborting
  res <- TmCalculator:::.tm_nn_chunk(
    list(sequence = c("RYSW", "ACGTACGT"),
         complement = c("YRSW", "TGCATGCA")),
    ambiguous = FALSE, shift = 0,
    nn_tbl = TmCalculator:::get_table("DNA_NN_SantaLucia_2004"),
    tmm_tbl = TmCalculator:::get_table("DNA_TMM_Bommarito_2000"),
    imm_tbl = TmCalculator:::get_table("DNA_IMM_Peyret_1999"),
    de_tbl = TmCalculator:::get_table("DNA_DE_Bommarito_2000"),
    end_tbl = matrix(numeric(0), nrow = 0, ncol = 2),
    dnac_high = 25, dnac_low = 25, self_comp = FALSE,
    Na = 50, K = 0, Tris = 0, Mg = 0, dNTPs = 0,
    salt_fn = "Schildkraut2010",
    DMSO = 0, formamide_unit = list(value = 0, unit = "percent"),
    dmso_factor = 0.75, formamide_factor = 0.65
  )
  expect_true(is.na(res$Tm[1]) && !is.na(res$Tm[2]))
})

test_that("result$df is computed lazily and matches result$gr", {
  input_seq <- c("ATGCGATGCGAAGGCGATGGCGTGTAGAATAGATCACATACTGCATAGCTGATC")
  result <- tm_calculate(input_seq, method = "tm_nn")
  expect_null(.subset2(result, "df"))          # not stored eagerly
  expect_s3_class(result$df, "data.frame")     # computed on demand
  expect_equal(nrow(result$df), length(result$gr))
})

test_that("Rcpp core matches R core across RNA/hybrid parameter sets", {
  seqs <- .rand_seqs(15, 30, 6)
  cseqs <- complement_fast(seqs)
  for (tab in c("DNA_NN_Breslauer_1986", "RNA_NN_Xia_1998",
                "RNA_DNA_NN_Sugimoto_1995", "DNA_NN_Weber_2015",
                "RNA_DNA_NN_Banerjee_2020", "RNA_NN_Zuber_2022")) {
    .compare_chunk(seqs, cseqs, nn_tbl = TmCalculator:::get_table(tab),
                   end_tbl = TmCalculator:::get_end_table(tab))
  }
})

test_that("tm_calculate end-to-end matches published regression values", {
  input_seq <- c("ATGCGATGCGAAGGCGATGGCGTGTAGAATAGATCACATACTGCATAGCTGATC",
                 "ATGCGATGCGCCCGGAGATAGAAGGCGTAGATACAGATCAGTAGCACCTTGAGAC")
  result <- tm_calculate(input_seq, method = "tm_nn")
  expect_equal(round(result$gr$Tm, 5), c(67.06562, 69.64434), tolerance = 0.1)
})

test_that("Zuber 2022 end effects are applied at both duplex ends", {
  end_tbl <- TmCalculator:::get_end_table("RNA_NN_Zuber_2022")
  expect_gt(nrow(end_tbl), 0)
  expect_equal(nrow(TmCalculator:::get_end_table("RNA_NN_Xia_1998")), 0)

  nn <- TmCalculator:::get_table("RNA_NN_Zuber_2022")
  args <- list(
    ambiguous = FALSE, shift = 0, nn_tbl = nn,
    tmm_tbl = TmCalculator:::get_table("DNA_TMM_Bommarito_2000"),
    imm_tbl = TmCalculator:::get_table("DNA_IMM_Peyret_1999"),
    de_tbl = TmCalculator:::get_table("RNA_DE_Turner_2010"),
    dnac_high = 25, dnac_low = 25, self_comp = FALSE,
    Na = 1000, K = 0, Tris = 0, Mg = 0, dNTPs = 0, salt_fn = "none",
    DMSO = 0, formamide_unit = list(value = 0, unit = "percent"),
    dmso_factor = 0.75, formamide_factor = 0.65)

  ## GC-ended duplex: no end penalty; AU-ended duplex: penalised at both ends
  gc_end <- list(sequence = "GCGCGCGC", complement = complement_fast("GCGCGCGC"))
  au_end <- list(sequence = "AGCGCGCT", complement = complement_fast("AGCGCGCT"))

  with_end <- function(chunk, et)
    do.call(TmCalculator:::.tm_nn_chunk, c(list(chunk), args, list(end_tbl = et)))
  no_end <- matrix(numeric(0), nrow = 0, ncol = 2)

  ## a GC-ended duplex is unaffected by the end table
  expect_equal(with_end(gc_end, end_tbl)$Tm, with_end(gc_end, no_end)$Tm)
  ## an AU-ended duplex is affected, and destabilised (end terms are positive)
  expect_false(isTRUE(all.equal(with_end(au_end, end_tbl)$Tm,
                                with_end(au_end, no_end)$Tm)))
  expect_lt(with_end(au_end, end_tbl)$Tm, with_end(au_end, no_end)$Tm)

  ## C++ and R reference paths agree with the end table active
  ref <- do.call(TmCalculator:::.tm_nn_chunk_r,
                 c(list(au_end), args, list(end_tbl = end_tbl)))
  expect_equal(with_end(au_end, end_tbl)$Tm, ref$Tm, tolerance = 1e-10)
})
