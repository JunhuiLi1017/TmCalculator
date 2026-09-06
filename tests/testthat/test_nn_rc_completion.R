test_that("reverse-complement completion copies from the reversed key", {
  # A key "XY/WZ" is the duplex 5'-XY-3' / 3'-WZ-5'. Read from the other
  # strand the same duplex is written "ZW/YX", the character reversal of the
  # key. Any completed row must therefore equal the row named by its own
  # reversal. This is checked here by re-deriving the relation rather than by
  # restating the lookup table in R/zzz.R, so that a transposition in that
  # table cannot be reproduced in the test that is supposed to catch it.
  rev_key <- function(k)
    vapply(strsplit(k, "", fixed = TRUE),
           function(ch) paste(rev(ch), collapse = ""), character(1))

  completed <- c("TT/AA", "AC/TG", "AG/TC", "TC/AG", "TG/AC", "CC/GG")

  # Every table that goes through .complete_nn_rc(). RNA_DNA_NN_Sugimoto_1995
  # and the RNA/DNA hybrid sets are excluded: they ship with all pairs and are
  # never completed.
  tables <- c(
    "DNA_NN_SantaLucia_2004", "DNA_NN_Breslauer_1986", "DNA_NN_Sugimoto_1996",
    "DNA_NN_Allawi_1998", "RNA_NN_Freier_1986", "RNA_NN_Xia_1998",
    "RNA_NN_Chen_2012", "RNA_NN_Zuber_2022",
    "DNA_NN_Weber_2015", "DNA_NN_Weber_OW04_69", "DNA_NN_Weber_OW04_119",
    "DNA_NN_Weber_OW04_220", "DNA_NN_Weber_OW04_621", "DNA_NN_Weber_OW04_1020",
    "RNA_NN_Weber_VIF_71", "RNA_NN_Weber_VIF_121", "RNA_NN_Weber_VIF_221",
    "RNA_NN_Weber_VIF_621", "RNA_NN_Weber_VIF_1021",
    "RNA_NN_Weber_FIF_71", "RNA_NN_Weber_FIF_121", "RNA_NN_Weber_FIF_221",
    "RNA_NN_Weber_FIF_621", "RNA_NN_Weber_FIF_1021",
    "DNA_NN_Ghosh_2020_PEG200", "RNA_NN_Ghosh_2023_PEG200"
  )

  for (nm in tables) {
    tbl <- TmCalculator:::get_table(nm)
    present <- intersect(completed, rownames(tbl))
    expect_true(length(present) > 0L,
                info = paste(nm, "has no completed rows"))
    for (k in present) {
      src <- rev_key(k)
      expect_true(src %in% rownames(tbl),
                  info = paste(nm, ": source row", src, "absent"))
      expect_equal(unname(tbl[k, ]), unname(tbl[src, ]),
                   info = paste(nm, ":", k, "should equal", src))
    }
  }
})

test_that("all 16 Watson-Crick stacks are present after completion", {
  wc <- c("AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT",
          "CG/GC", "GC/CG", "GG/CC", "TT/AA", "AC/TG", "AG/TC", "TC/AG",
          "TG/AC", "CC/GG")
  tbl <- TmCalculator:::get_table("DNA_NN_SantaLucia_2004")
  expect_true(all(wc %in% rownames(tbl)))
})

test_that("the four transposed rows carry their published values", {
  # Regression test for the transposition present up to 1.10.0, pinned to
  # literal values from SantaLucia & Hicks (2004) so that it fails loudly if
  # the mapping is ever reverted.
  tbl <- TmCalculator:::get_table("DNA_NN_SantaLucia_2004")
  expect_equal(unname(tbl["TG/AC", ]), c(-8.5, -22.7))  # = CA/GT
  expect_equal(unname(tbl["AC/TG", ]), c(-8.4, -22.4))  # = GT/CA
  expect_equal(unname(tbl["AG/TC", ]), c(-7.8, -21.0))  # = CT/GA
  expect_equal(unname(tbl["TC/AG", ]), c(-8.2, -22.2))  # = GA/CT
})
