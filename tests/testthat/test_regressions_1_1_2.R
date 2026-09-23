Sys.setenv(R_TESTS = "")

# ===========================================================================
# Regression tests for every defect fixed in 1.1.2.
#
# One block per defect, numbered as in NEWS.md. Each states what was wrong,
# then asserts the behaviour that must now hold. Blocks 6 to 10 concern the
# stack lookup and are the ones a future refactor is most likely to undo, so
# they are stated twice: once as a frozen value, and once as the property the
# value follows from.
#
# The frozen values all use DNA_NN_Allawi_1998 at dnac_high = 250, dnac_low =
# 0 and salt_method = "none", so that only the stacking sum is compared and
# nothing downstream can mask a change.
# ===========================================================================

flip <- function(s) vapply(strsplit(s, "", fixed = TRUE),
                           function(x) paste(rev(x), collapse = ""), character(1))

SEQ  <- "GCATCGTAGGCTAGCT"
PERF <- chartr("ACGT", "TGCA", SEQ)          # exact complement, 3' to 5'

tm_of <- function(s, cmp, shift = 0, ...) {
  gr <- TmCalculator::to_genomic_ranges(s, complement_seq = cmp)
  r  <- TmCalculator::tm_nn(gr, nn_table = "DNA_NN_Allawi_1998", shift = shift,
                            dnac_high = 250, dnac_low = 0,
                            salt_method = "none", ...)
  as.numeric(GenomicRanges::mcols(r$gr)$Tm)
}

rand_duplex <- function() {
  s <- paste(sample(c("A", "C", "G", "T"), sample(18:24, 1), replace = TRUE),
             collapse = "")
  list(s = s, c = chartr("ACGT", "TGCA", s))
}

mutate_at <- function(cmp, i) {
  ch <- strsplit(cmp, "", fixed = TRUE)[[1]]
  ch[i] <- sample(setdiff(c("A", "C", "G", "T"), ch[i]), 1)
  paste(ch, collapse = "")
}

# ---------------------------------------------------------------------------
# 1. Owczarzy2008 returned NA whenever magnesium dominated.
#
# The correction is piecewise in R = sqrt([Mg2+]free)/[Mon] and the third
# regime, R >= 6, had no branch: the scalar function fell off the end of its
# if/else without assigning a result and the vectorized one returned NA. The
# regime is ordinary for a low-monovalent buffer -- Na = 1 mM with Mg = 5 mM
# gives R = 71 -- so Tm and GC both came back NA, silently. Present since the
# first release; reported as issue #9.
# ---------------------------------------------------------------------------

test_that("Owczarzy2008 is finite in all three regimes", {
  seq20 <- "GCATCGTAGGCTAGCTTGCA"                      # 20 nt, 55% GC
  monovalent <- TmCalculator::salt_correct(Na = 50, Mg = 0,
                                           method = "Owczarzy2008",
                                           input_seq = seq20)
  competing  <- TmCalculator::salt_correct(Na = 50, Mg = 1.5,
                                           method = "Owczarzy2008",
                                           input_seq = seq20)
  divalent   <- TmCalculator::salt_correct(Na = 1, Mg = 5,
                                           method = "Owczarzy2008",
                                           input_seq = seq20)
  expect_equal(monovalent, 1.3200659324611811e-04, tolerance = 1e-9)
  expect_equal(competing,  7.4877849736476802e-05, tolerance = 1e-9)
  expect_equal(divalent,   5.6022561348541469e-05, tolerance = 1e-9)

  res <- TmCalculator::tm_nn(TmCalculator::to_genomic_ranges(seq20),
                             salt_method = "Owczarzy2008", Na = 1, Mg = 5)
  expect_equal(as.numeric(GenomicRanges::mcols(res$gr)$Tm), 63.162,
               tolerance = 1e-3)
  expect_equal(as.numeric(GenomicRanges::mcols(res$gr)$GC), 55)
})

test_that("the scalar and vectorized salt corrections agree", {
  seq20 <- "GCATCGTAGGCTAGCTTGCA"
  for (cond in list(list(Na = 50, Mg = 0), list(Na = 50, Mg = 1.5),
                    list(Na = 1, Mg = 5), list(Na = 1, Mg = 5, dNTPs = 10),
                    list(Na = 0, Mg = 5), list(Na = 0, Mg = 0))) {
    args <- utils::modifyList(
      list(Na = 0, K = 0, Tris = 0, Mg = 0, dNTPs = 0), cond)
    scalar <- do.call(TmCalculator::salt_correct,
                      c(args, list(method = "Owczarzy2008",
                                   input_seq = seq20)))
    vec <- do.call(TmCalculator:::.salt_correct_vec,
                   c(args, list(method = "Owczarzy2008",
                                gc_pct = 55, seq_len = 20)))
    expect_equal(vec, scalar, tolerance = 1e-9)
  }
})

# ---------------------------------------------------------------------------
# 2. Owczarzy2008 with magnesium but no monovalent cation applied nothing.
#
# Mon == 0 fell into the guard that returns a zero correction, which is right
# for the six methods that take log([Mon]) and wrong for this one: in the
# divalent-dominated regime [Mon] drops out of the expression entirely, so the
# correction is defined and a solution with no monovalent cation gives the
# same answer as one with a trace of it. Worse than the NA above, because
# nothing marked the result. Present since the first release.
# ---------------------------------------------------------------------------

test_that("magnesium alone is corrected for, not ignored", {
  seq20 <- "GCATCGTAGGCTAGCTTGCA"
  only_mg <- TmCalculator::salt_correct(Na = 0, Mg = 5,
                                        method = "Owczarzy2008",
                                        input_seq = seq20)
  expect_false(only_mg == 0)
  expect_equal(only_mg,
               TmCalculator::salt_correct(Na = 1, Mg = 5,
                                          method = "Owczarzy2008",
                                          input_seq = seq20),
               tolerance = 1e-9)
  # nothing at all in solution is the one case with nothing to correct for
  expect_equal(TmCalculator::salt_correct(Na = 0, Mg = 0,
                                          method = "Owczarzy2008",
                                          input_seq = seq20), 0)

  # dNTPs bind magnesium with Ka = 3e4 /M, so 5 mM dNTPs against 5 mM Mg
  # leaves 0.39 mM free: the answer must stay finite and move back towards
  # the magnesium-free one rather than diverge
  chelated   <- TmCalculator::salt_correct(Na = 50, Mg = 5, dNTPs = 5,
                                           method = "Owczarzy2008",
                                           input_seq = seq20)
  mg_free    <- TmCalculator::salt_correct(Na = 50, Mg = 0,
                                           method = "Owczarzy2008",
                                           input_seq = seq20)
  unchelated <- TmCalculator::salt_correct(Na = 50, Mg = 5,
                                           method = "Owczarzy2008",
                                           input_seq = seq20)
  expect_true(is.finite(chelated))
  expect_lt(abs(chelated - mg_free), abs(unchelated - mg_free))
})

# ---------------------------------------------------------------------------
# 3. tm_gc accepted the Owczarzy corrections and added them to Tm.
#
# They correct the reciprocal of the melting temperature in kelvin, referenced
# to the same duplex in 1 M Na+, and carry a 1/(2(N-1)) duplex-length term of
# their own; a GC-content formula is on neither footing and already has a
# length term. Adding a quantity of order 1e-5 K^-1 to a Celsius value left
# the Tm essentially uncorrected. Present since the first release, reachable
# only through `userset`.
# ---------------------------------------------------------------------------

test_that("tm_gc refuses the Owczarzy corrections", {
  gr <- TmCalculator::to_genomic_ranges(SEQ)
  coef <- c(81.5, 0.41, 675, 1)
  expect_error(TmCalculator::tm_gc(gr, userset = coef,
                                   salt_method = "Owczarzy2008"),
               "not available for tm_gc")
  expect_error(TmCalculator::tm_gc(gr, userset = coef,
                                   salt_method = "Owczarzy2004"),
               "not available for tm_gc")
  expect_error(TmCalculator::tm_calculate(SEQ, method = "tm_gc",
                                          userset = coef,
                                          salt_method = "Owczarzy2008"),
               "not available for method")
})

# ---------------------------------------------------------------------------
# 4. tm_gc silently overrode salt_method, and neither NULL nor NA worked.
#
# With a built-in `variant` the salt term is part of the published formula, so
# a different one was discarded without a word, while the documented NULL and
# NA spellings both errored. Naming a different method is still overridden --
# the formula has to be the one it is labelled as -- but now says so; NA and
# "none" drop the correction, which is not a substitution and is honoured.
# ---------------------------------------------------------------------------

test_that("tm_gc says when a named salt_method is overridden", {
  gr <- TmCalculator::to_genomic_ranges(SEQ)
  expect_warning(TmCalculator::tm_gc(gr, variant = "vonAhsen2001",
                                     salt_method = "Wetmur1991"),
                 "is ignored")
  expect_warning(TmCalculator::tm_gc(gr, variant = "Chester1993",
                                     salt_method = "Wetmur1991"),
                 "no salt term of its own")
  # the variant's own method, and leaving it alone, are both silent
  expect_warning(TmCalculator::tm_gc(gr, variant = "vonAhsen2001",
                                     salt_method = "SantaLucia1998-1"),
                 regexp = NA)
  expect_warning(TmCalculator::tm_gc(gr, variant = "vonAhsen2001"),
                 regexp = NA)
  # and the published formula is what was actually used
  own <- TmCalculator::tm_gc(gr, variant = "vonAhsen2001", Na = 50)
  expect_equal(own$options[["Salt correction"]], "SantaLucia1998-1")
})

test_that("tm_gc honours NA and 'none' as no correction", {
  gr <- TmCalculator::to_genomic_ranges(SEQ)
  expect_warning(bare <- TmCalculator::tm_gc(gr, variant = "Primer3Plus",
                                             salt_method = NA, Na = 50),
                 regexp = NA)
  expect_true(is.na(bare$options[["Salt correction"]]))
  salted <- TmCalculator::tm_gc(gr, variant = "Primer3Plus", Na = 50)
  expect_equal(GenomicRanges::mcols(salted$gr)$Tm -
                 GenomicRanges::mcols(bare$gr)$Tm,
               16.6 * log10(0.05), tolerance = 1e-9)
  expect_equal(GenomicRanges::mcols(
                 TmCalculator::tm_gc(gr, variant = "Primer3Plus",
                                     salt_method = "none", Na = 50)$gr)$Tm,
               GenomicRanges::mcols(bare$gr)$Tm)
  expect_error(TmCalculator::tm_gc(gr, userset = c(81.5, 0.41, 675, 1),
                                   salt_method = c("Wetmur1991",
                                                   "SantaLucia1996")),
               "single method name")
})

# ---------------------------------------------------------------------------
# 5. GC followed Tm into NA, and NA was returned silently.
#
# Base composition is a property of the sequence, not of the thermodynamic
# model, so it survives a Tm the model cannot produce. And an NA that nobody
# is told about is carried into a mean or a plot unnoticed.
# ---------------------------------------------------------------------------

test_that("GC survives a Tm the model cannot compute, and NA is announced", {
  gr <- TmCalculator::to_genomic_ranges(c("G", SEQ))
  expect_warning(res <- TmCalculator::tm_nn(gr), "Tm is NA for 1 region")
  expect_true(is.na(GenomicRanges::mcols(res$gr)$Tm[1]))
  expect_equal(GenomicRanges::mcols(res$gr)$GC[1], 100)
  expect_false(is.na(GenomicRanges::mcols(res$gr)$Tm[2]))
})

# ---------------------------------------------------------------------------
# 6. A stack stored in the reversed spelling was not found, and contributed
#    zero.
#
# A key "XY/WZ" is 5'-XY-3' over 3'-WZ-5'; the same stack read from the other
# strand is the whole key reversed, "ZW/YX". The published tables store each
# stack in one orientation only (85 of the 87 keys of DNA_IMM_Peyret_1999
# have no reversed twin), and the lookup was a single probe. Exactly one of
# the two stacks flanking any internal mismatch needs the reversed spelling,
# so no internal mismatch was ever scored completely.
#
# The T.C mismatch at position 12 below flanks CT/GC (stored as written) and
# TA/CT (stored only as TC/AT). Dropping the second understated the penalty by
# 45%: 3.67 C instead of 6.62 C.
# ---------------------------------------------------------------------------

test_that("both stacks flanking an internal mismatch are counted", {
  perfect <- tm_of(SEQ, PERF)
  mm      <- tm_of(SEQ, "CGTAGCATCCGCTCGA")
  expect_equal(perfect, 66.5968, tolerance = 1e-4)
  expect_equal(mm,      59.9731, tolerance = 1e-4)
  expect_lt(mm, perfect)
})

test_that("the four calls reported in issue #10 agree with each other", {
  # Verbatim from the report, including its arguments. The first two describe
  # one duplex with a T.C mismatch at position 12, read from either strand;
  # the last two are the same duplex without the mismatch. The report observed
  # 62.927 against 61.633 for the mismatched pair -- the same molecule giving
  # two answers 1.3 C apart -- while the perfect pair already agreed at 66.597.
  f <- function(p, cmp) {
    r <- TmCalculator::tm_nn(
      TmCalculator::to_genomic_ranges(p, complement_seq = cmp),
      nn_table = "DNA_NN_Allawi_1998", dnac_high = 250, dnac_low = 0,
      salt_method = "none", Na = 1000)
    as.numeric(S4Vectors::mcols(r$gr)$Tm)
  }
  mm1 <- f("GCATCGTAGGCTAGCT", "CGTAGCATCCGCTCGA")   # read from strand 1
  mm2 <- f("AGCTCGCCTACGATGC", "TCGATCGGATGCTACG")   # read from strand 2
  ok1 <- f("GCATCGTAGGCTAGCT", "CGTAGCATCCGATCGA")   # perfect, strand 1
  ok2 <- f("AGCTAGCCTACGATGC", "TCGATCGGATGCTACG")   # perfect, strand 2

  expect_equal(mm1, mm2, tolerance = 1e-10)          # was 62.927 vs 61.633
  expect_equal(ok1, ok2, tolerance = 1e-10)
  expect_equal(mm1, 59.9731, tolerance = 1e-4)       # the full penalty
  expect_equal(ok1, 66.5968, tolerance = 1e-4)       # unchanged since 1.1.1
  expect_lt(mm1, ok1)
})

# ---------------------------------------------------------------------------
# 7. The terminal-mismatch key was built in the wrong orientation, at both
#    ends, which both missed real terminal mismatches and invented others.
#
# The TMM tables carry the penultimate pair first and the terminal pair second
# ("AA/TA" is a Watson-Crick pair then a mismatch). The walk built the key
# with the terminal pair first. A genuine terminal mismatch therefore never
# matched. A duplex whose terminal pair IS Watson-Crick but whose next pair is
# not produces a terminal-first key of exactly the shape the table stores, so
# it matched spuriously and was charged a penalty it had not earned.
# ---------------------------------------------------------------------------

TMM_BENT <- local({
  b <- TmCalculator:::get_table("DNA_TMM_Bommarito_2000")
  b[, 1] <- b[, 1] - 5                           # 5 kcal/mol more stable
  b
})

# `SEQ` is GCATCGTAGGCTAGCT, so position 1 pairs G and position 16 pairs T.
# Changing the complement at either end makes a terminal mismatch there; the
# `mm_at_2` and `mm_at_15` cases move the mismatch one position inwards,
# leaving a Watson-Crick terminus, and must not reach the table at all.
# The two key columns are the spellings the table is now probed with: the
# reversal of the first stack at the left end, the last stack as-is at the
# right end.
#
#   case          complement           left key   right key  in table
#   left  G.A     A + PERF[2:]         GA/CG      CT/GA      left
#   left  G.G     G + PERF[2:]         GG/CG      CT/GA      left
#   left  G.T     T + PERF[2:]         GT/CG      CT/GA      left
#   right T.C     PERF[1:15] + C       GC/CG      CT/GC      right
#   right T.G     PERF[1:15] + G       GC/CG      CT/GG      right
#   right T.T     PERF[1:15] + T       GC/CG      CT/GT      right
#   both ends     A + PERF[2:15] + C   GA/CG      CT/GC      both
#   mm_at_2       PERF with 2 -> A     AC/CG      CT/GA      neither
#   mm_at_15      PERF with 15 -> T    GC/CG      CT/TA      neither
#   perfect       PERF                 GC/CG      CT/GA      neither

at <- function(...) {
  x <- strsplit(PERF, "", fixed = TRUE)[[1]]
  for (p in list(...)) x[p[[1]]] <- p[[2]]
  paste(x, collapse = "")
}

TMM_CASES <- list(
  list(name = "left G.A",  cmp = at(list(1, "A")),  tm = 63.7946, hits = TRUE),
  list(name = "left G.G",  cmp = at(list(1, "G")),  tm = 46.1302, hits = TRUE),
  list(name = "left G.T",  cmp = at(list(1, "T")),  tm = 63.5560, hits = TRUE),
  list(name = "right T.C", cmp = at(list(16, "C")), tm = 65.1777, hits = TRUE),
  list(name = "right T.G", cmp = at(list(16, "G")), tm = 65.0678, hits = TRUE),
  list(name = "right T.T", cmp = at(list(16, "T")), tm = 65.3648, hits = TRUE),
  list(name = "both ends", cmp = at(list(1, "A"), list(16, "C")),
                                                    tm = 62.2927, hits = TRUE),
  list(name = "mm_at_2",   cmp = at(list(2, "A")),  tm = 55.1589, hits = FALSE),
  list(name = "mm_at_15",  cmp = at(list(15, "T")), tm = 55.4735, hits = FALSE),
  list(name = "perfect",   cmp = PERF,              tm = 66.5968, hits = FALSE)
)

test_that("the terminal mismatch table is reached exactly when it should be", {
  for (k in TMM_CASES) {
    plain <- tm_of(SEQ, k$cmp)
    moved <- tm_of(SEQ, k$cmp, tmm_table = TMM_BENT)
    expect_equal(plain, k$tm, tolerance = 1e-4, info = k$name)
    if (k$hits) {
      # a genuine terminal mismatch: perturbing the table must move the answer
      expect_false(isTRUE(all.equal(plain, moved)), info = k$name)
    } else {
      # a Watson-Crick terminus: the table must not be consulted at all
      expect_equal(plain, moved, info = k$name)
    }
  }
})

test_that("every terminal mismatch reads the same from either strand", {
  for (k in TMM_CASES) {
    expect_equal(tm_of(SEQ, k$cmp),
                 tm_of(flip(k$cmp), flip(SEQ)),
                 tolerance = 1e-10, info = k$name)
  }
})

test_that("the terminal mismatch identity changes the answer", {
  # Before the orientation was corrected the left-hand terminal stack matched
  # nothing in either table, so G.A, G.G and G.T at position 1 all returned
  # the same 62.1434 -- the stack contributed nothing whatever the mismatch
  # was. The three must now differ from each other and from the perfect duplex.
  left <- vapply(TMM_CASES[1:3], function(k) tm_of(SEQ, k$cmp), numeric(1))
  expect_length(unique(round(left, 6)), 3L)
  expect_true(all(left < tm_of(SEQ, PERF)))
})

test_that("a perfect duplex never reaches the terminal mismatch table", {
  set.seed(11)
  for (i in 1:20) {
    d <- rand_duplex()
    expect_equal(tm_of(d$s, d$c), tm_of(d$s, d$c, tmm_table = TMM_BENT))
  }
})

# ---------------------------------------------------------------------------
# 8. A stack present in both nn_table and imm_table was counted twice.
#
# A stack has one delta_H and delta_S. The two tables were consulted
# independently and both added. The overlap is the G.U wobble stacks of the
# RNA sets, whose spellings also occur in the default DNA_IMM_Peyret_1999
# meaning a DNA G.T mismatch, so such a stack was charged an RNA wobble
# parameter plus an unrelated DNA mismatch parameter. The nearest-neighbor
# set now wins. No DNA set overlaps the mismatch table in either orientation.
# ---------------------------------------------------------------------------

test_that("no DNA stack is defined by both the nn and the mismatch table", {
  # the premise of "DNA is unaffected": if this ever stops holding, the
  # precedence rule starts moving DNA values and must be revisited
  rev_str <- function(k) paste(rev(strsplit(k, "", fixed = TRUE)[[1]]),
                               collapse = "")
  imm <- rownames(TmCalculator:::get_table("DNA_IMM_Peyret_1999"))
  for (nm in grep("^DNA_NN_", names(TmCalculator:::.TM_CONSTANTS),
                  value = TRUE)) {
    k <- rownames(TmCalculator:::get_table(nm))
    k <- grep("^[ACGTI]{2}/[ACGTI]{2}$", k, value = TRUE)
    expect_length(intersect(k, imm), 0L)
    expect_length(intersect(vapply(k, rev_str, character(1L),
                                   USE.NAMES = FALSE), imm), 0L)
  }
})

test_that("a stack the nn table defines is not taken from the mismatch table", {
  # An RNA wobble spelled the same way in both: bending the mismatch table
  # must not move an RNA duplex whose wobble the RNA set already defines.
  cmp <- PERF
  substr(cmp, 6, 6) <- "T"       # position 6 is G, so G.T: a wobble in RNA
  # both stacks flanking it, CG/GT and GT/TA, are in RNA_NN_Chen_2012 and in
  # DNA_IMM_Peyret_1999 alike
  bent <- TmCalculator:::get_table("DNA_IMM_Peyret_1999")
  bent[, 1] <- bent[, 1] - 5
  run <- function(imm) {
    gr <- TmCalculator::to_genomic_ranges(SEQ, complement_seq = cmp)
    r <- TmCalculator::tm_nn(gr, nn_table = "RNA_NN_Chen_2012",
                             imm_table = imm, de_table = "RNA_DE_Turner_2010",
                             dnac_high = 250, dnac_low = 0,
                             salt_method = "none")
    as.numeric(GenomicRanges::mcols(r$gr)$Tm)
  }
  # the RNA set defines these stacks, so the DNA mismatch table must not be
  # consulted for them and bending it must change nothing
  expect_equal(run("DNA_IMM_Peyret_1999"), run(bent))
})

# ---------------------------------------------------------------------------
# 9. init_5T/A was charged for the top strand only.
#
# It is due once for each strand whose 5' end is T: the first base of the
# sequence, and the last base of the complement. Charging only the first made
# Tm depend on which strand was handed over. The row is zero in every shipped
# parameter set, so no shipped value ever moved -- which is why this needs a
# table carrying a non-zero value to show at all.
# ---------------------------------------------------------------------------

test_that("the 5'-T penalty is charged per strand", {
  tb <- TmCalculator:::get_table("DNA_NN_Allawi_1998")
  tb["init_5T/A", 1:2] <- c(-3.5, -9.0)
  run <- function(s, cmp) {
    gr <- TmCalculator::to_genomic_ranges(s, complement_seq = cmp)
    r  <- TmCalculator::tm_nn(gr, nn_table = tb, dnac_high = 250,
                              dnac_low = 0, salt_method = "none")
    as.numeric(GenomicRanges::mcols(r$gr)$Tm)
  }
  s   <- "TCATCGTAGGCTAGCA"                    # 5'-T on both strands
  cmp <- chartr("ACGT", "TGCA", s)
  expect_equal(run(s, cmp), run(flip(cmp), flip(s)), tolerance = 1e-10)

  tb0 <- tb; tb0["init_5T/A", 1:2] <- c(0, 0)
  gr <- TmCalculator::to_genomic_ranges(s, complement_seq = cmp)
  zeroed <- TmCalculator::tm_nn(gr, nn_table = tb0, dnac_high = 250,
                                dnac_low = 0, salt_method = "none")
  expect_false(isTRUE(all.equal(
    run(s, cmp), as.numeric(GenomicRanges::mcols(zeroed$gr)$Tm))))
})

test_that("init_5T/A is zero in every shipped parameter set", {
  # the premise of "no shipped value moved" in blocks 9 above
  nms <- grep("_NN_", names(TmCalculator:::.TM_CONSTANTS), value = TRUE)
  vals <- do.call(rbind, lapply(nms, function(n) {
    tb <- TmCalculator:::get_table(n)
    if ("init_5T/A" %in% rownames(tb)) tb["init_5T/A", 1:2] else NULL
  }))
  expect_gt(nrow(vals), 25)
  expect_true(all(vals == 0))
})

# ---------------------------------------------------------------------------
# 10. The property that blocks 6 to 9 all violated.
#
# A duplex read from either strand is the same molecule, so its Tm must not
# depend on which strand was handed over as `sequence`:
#
#     5'-G C A T-3'  sequence        5'-A T G C-3'  = rev(complement)
#     3'-C G T A-5'  complement      3'-T A C G-5'  = rev(sequence)
#
# This needs no reference implementation and no frozen value. It would also
# have caught the table transposition fixed in 1.1.0. Before 1.1.2 the two
# readings differed by up to about 3 C for an internal mismatch and 4 C for a
# terminal one.
# ---------------------------------------------------------------------------

both_ways <- function(s, cmp, shift = 0) {
  c(tm_of(s, cmp, shift), tm_of(flip(cmp), flip(s), -shift))
}

test_that("Tm does not depend on which strand is passed", {
  cases <- list(
    perfect = function(d) d$c,
    internal_mismatch = function(d) mutate_at(d$c, sample(4:(nchar(d$s) - 3), 1)),
    terminal_mismatch = function(d) mutate_at(d$c, sample(c(1L, nchar(d$s)), 1)),
    beside_terminus   = function(d) mutate_at(d$c, sample(c(2L, nchar(d$s) - 1L), 1)),
    two_mismatches    = function(d) mutate_at(mutate_at(d$c, 5L), nchar(d$s) - 4L)
  )
  set.seed(7)
  for (nm in names(cases)) {
    for (i in 1:25) {
      d <- rand_duplex()
      x <- both_ways(d$s, cases[[nm]](d))
      expect_equal(x[1], x[2], tolerance = 1e-10,
                   info = paste(nm, "iteration", i))
    }
  }
})

test_that("Tm does not depend on which strand is passed, with a dangling end", {
  # the complement is one base short at its 3' end, so the sequence overhangs
  # at the left; read from the other strand the overhang is at the right and
  # the padding comes from the length difference rather than from `shift`
  set.seed(8)
  for (i in 1:25) {
    d <- rand_duplex()
    short <- substring(d$c, 2)
    expect_equal(tm_of(d$s, short, shift = -1),
                 tm_of(flip(short), flip(d$s), shift = 0),
                 tolerance = 1e-10)
  }
})

test_that("perfect duplexes are untouched by all of the above", {
  # the values the manuscript's benchmarks and the cross-tool comparison rest
  # on: every stack of a perfect duplex is Watson-Crick, and the NN tables
  # have carried both orientations of those since 1.1.0
  expect_equal(tm_of(SEQ, PERF), 66.5968, tolerance = 1e-4)
  expect_equal(tm_of("AGCTAGCCTACGATGC", "TCGATCGGATGCTACG"), 66.5968,
               tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# 11. complement_seq did not say which direction it wanted, and
#     generate_complement's documentation named the two directions the wrong
#     way round.
#
# `complement_seq` wants the plain complement, aligned base for base and so
# written 3' to 5'. The reverse complement is the same strand written 5' to
# 3', which is what Biostrings::reverseComplement() gives; passing it pairs
# every position against the wrong base and returns a large negative Tm rather
# than an error.
# ---------------------------------------------------------------------------

test_that("a reverse complement passed as complement_seq is called out", {
  expect_warning(TmCalculator::to_genomic_ranges(SEQ,
                                                 complement_seq = flip(PERF)),
                 "looks like a reverse complement")
  # the correct form, and a genuinely mismatched duplex, are both silent
  expect_warning(TmCalculator::to_genomic_ranges(SEQ, complement_seq = PERF),
                 regexp = NA)
  expect_warning(
    TmCalculator::to_genomic_ranges(SEQ, complement_seq = "CGTAGCATCCGCTCGA"),
    regexp = NA)
})

test_that("generate_complement's two directions are what the docs claim", {
  # reverse = FALSE pairs position by position; reverse = TRUE is its reversal
  expect_equal(TmCalculator::generate_complement("ATGCG"), "TACGC")
  expect_equal(TmCalculator::generate_complement("ATGCG", reverse = TRUE),
               "CGCAT")
  expect_equal(TmCalculator::generate_complement(SEQ), PERF)
  expect_equal(TmCalculator::generate_complement(SEQ, reverse = TRUE),
               flip(PERF))
})
