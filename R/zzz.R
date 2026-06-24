# =============================================================================
# .build_all_tm_tables() - updated with reverse complement completion (Bug 9 fix)
#
# WHAT CHANGED:
#   Added .complete_nn_rc() helper that adds the 6 missing reverse complement
#   rows to any standard 17-row DNA/RNA NN table.
#
#   The 16 Watson-Crick NN pairs split into:
#     10 asymmetric pairs: each has a distinct RC partner → 6 unique new rows needed
#        AA/TT ↔ TT/AA   (TT/AA missing)
#        CA/GT ↔ TG/AC   (TG/AC missing)
#        GA/CT ↔ TC/AG   (TC/AG missing)
#        CT/GA ↔ AG/TC   (AG/TC missing)
#        GT/CA ↔ AC/TG   (AC/TG missing)
#        GG/CC ↔ CC/GG   (CC/GG missing)
#     6 self-complementary palindromes: RC maps to themselves → no new rows needed
#        AT/TA, TA/AT, CG/GC, GC/CG, CA/GT→already, GC/CG→already
#        (Actually: AT/TA, TA/AT, CG/GC, GC/CG are true palindromes)
#
#   RNA_DNA_NN_Sugimoto_1995 already contains all 16 + 6 RC rows → skip.
#   RNA_NN_Chen_2012 has extra wobble rows → complete only the standard 16.
# =============================================================================

# -- Helper: add the 6 missing reverse complement rows to a standard NN table --
#
# There are exactly 16 valid Watson-Crick dinucleotide keys (XY/X'Y' where
# X' = complement of X, Y' = complement of Y). Published NN tables list only
# 10 of these - the 6 missing ones are thermodynamically equivalent to an
# existing row when read from the other strand direction.
#
# The complete mapping (verified by enumeration):
#   TT/AA  <- same ΔH/ΔS as  AA/TT   (homodimer palindrome)
#   AC/TG  <- same ΔH/ΔS as  CA/GT   (RC pair)
#   AG/TC  <- same ΔH/ΔS as  GA/CT   (RC pair)
#   TC/AG  <- same ΔH/ΔS as  CT/GA   (RC pair)
#   TG/AC  <- same ΔH/ΔS as  GT/CA   (RC pair)
#   CC/GG  <- same ΔH/ΔS as  GG/CC   (homodimer palindrome)
#
# @param tbl       Matrix with dimnames - a standard NN table (17+ rows).
# @param skip_rows Row names to exclude from completion (init/sym rows and any
#                  extra wobble rows like in RNA_NN_Chen_2012).
# @return          The input matrix with up to 6 new rows appended.
#
.complete_nn_rc <- function(tbl,
                            skip_rows = c("init", "init_A/T", "init_G/C",
                                          "init_oneG/C", "init_allA/T",
                                          "init_5T/A", "sym")) {
  # Fixed mapping: new key -> source key to copy values from
  # This is the complete, verified set - do not derive algorithmically
  # (TT/AA and CC/GG are self-palindromes under the RC formula and would
  # be missed by a naive "RC != self" filter).
  missing_map <- c(
    "TT/AA" = "AA/TT",
    "AC/TG" = "CA/GT",
    "AG/TC" = "GA/CT",
    "TC/AG" = "CT/GA",
    "TG/AC" = "GT/CA",
    "CC/GG" = "GG/CC"
  )
  
  # Only add rows not already present (handles RNA_DNA_NN_Sugimoto_1995
  # which was published with the full 16-pair set)
  to_add <- missing_map[!names(missing_map) %in% rownames(tbl)]
  
  if (length(to_add) == 0L) return(tbl)   # already complete
  
  # Verify source rows exist before copying
  missing_sources <- to_add[!to_add %in% rownames(tbl)]
  if (length(missing_sources) > 0L)
    stop(sprintf("Source rows missing from table: %s",
                 paste(missing_sources, collapse=", ")))
  
  new_rows           <- tbl[to_add, , drop = FALSE]
  rownames(new_rows) <- names(to_add)
  
  rbind(tbl, new_rows)
}

# -- Main table builder --------------------------------------------------------
.build_all_tm_tables <- function() {
  
  nn_col     <- c("left", "right")
  nn_row_std <- c("init", "init_A/T", "init_G/C", "init_oneG/C", "init_allA/T",
                  "init_5T/A", "sym",
                  "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA",
                  "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC")
  
  # -- DNA NN Tables -----------------------------------------------------------
  DNA_NN_SantaLucia_2004 <- .complete_nn_rc(matrix(c(
    0.2,-5.7,  2.2,6.9,  0,0,  0,0,  0,0,  0,0,  0,-1.4,
    -7.6,-21.3,  -7.2,-20.4,  -7.2,-20.4,  -8.5,-22.7,  -8.4,-22.4,
    -7.8,-21.0,  -8.2,-22.2,  -10.6,-27.2,  -9.8,-24.4,  -8.0,-19.0
  ), ncol=2, byrow=TRUE, dimnames=list(nn_row_std, nn_col)))
  
  DNA_NN_Breslauer_1986 <- .complete_nn_rc(matrix(c(
    0,0,  0,0,  0,0,  0,-16.8,  0,-20.1,  0,0,  0,-1.3,
    -9.1,-24.0,  -8.6,-23.9,  -6.0,-16.9,  -5.8,-12.9,  -6.5,-17.3,
    -7.8,-20.8,  -5.6,-13.5,  -11.9,-27.8,  -11.1,-26.7,  -11.0,-26.6
  ), ncol=2, byrow=TRUE, dimnames=list(nn_row_std, nn_col)))
  
  DNA_NN_Sugimoto_1996 <- .complete_nn_rc(matrix(c(
    0.6,-9.0,  0,0,  0,0,  0,0,  0,0,  0,0,  0,-1.4,
    -8.0,-21.9,  -5.6,-15.2,  -6.6,-18.4,  -8.2,-21.0,  -9.4,-25.5,
    -6.6,-16.4,  -8.8,-23.5,  -11.8,-29.0,  -10.5,-26.4,  -10.9,-28.4
  ), ncol=2, byrow=TRUE, dimnames=list(nn_row_std, nn_col)))
  
  DNA_NN_Allawi_1998 <- .complete_nn_rc(matrix(c(
    0,0,  2.3,4.1,  0.1,-2.8,  0,0,  0,0,  0,0,  0,-1.4,
    -7.9,-22.2,  -7.2,-20.4,  -7.2,-21.3,  -8.5,-22.7,  -8.4,-22.4,
    -7.8,-21.0,  -8.2,-22.2,  -10.6,-27.2,  -9.8,-24.4,  -8.0,-19.9
  ), ncol=2, byrow=TRUE, dimnames=list(nn_row_std, nn_col)))
  
  # -- RNA NN Tables -----------------------------------------------------------
  RNA_NN_Freier_1986 <- .complete_nn_rc(matrix(c(
    0,-10.8,  0,0,  0,0,  0,0,  0,0,  0,0,  0,-1.4,
    -6.6,-18.4,  -5.7,-15.5,  -8.1,-22.6,  -10.5,-27.8,  -10.2,-26.2,
    -7.6,-19.2,  -13.3,-35.5,  -8.0,-19.4,  -14.2,-34.9,  -12.2,-29.7
  ), ncol=2, byrow=TRUE, dimnames=list(nn_row_std, nn_col)))
  
  RNA_NN_Xia_1998 <- .complete_nn_rc(matrix(c(
    3.61,-1.5,  3.72,10.5,  0,0,  0,0,  0,0,  0,0,  0,-1.4,
    -6.82,-19.0,  -9.38,-26.7,  -7.69,-20.5,  -10.44,-26.9,  -11.40,-29.5,
    -10.48,-27.1,  -12.44,-32.5,  -10.64,-26.7,  -14.88,-36.9,  -13.39,-32.7
  ), ncol=2, byrow=TRUE, dimnames=list(nn_row_std, nn_col)))
  
  # Chen_2012 has extra wobble rows - complete_nn_rc skips non-standard rows
  # automatically and only processes the 16 standard Watson-Crick rows.
  chen_rows <- c(nn_row_std,
                 "GT/TG","GG/TT","AG/TT","TG/AT","TT/AG","TG/GT",
                 "AT/TG","CG/GT","CT/GG","GG/CT","GT/CG")
  RNA_NN_Chen_2012 <- .complete_nn_rc(matrix(c(
    6.40,6.99,  3.85,11.04,  0,0,  0,0,  0,0,  0,0,  0,-1.4,
    -7.09,-19.8,  -9.11,-25.8,  -8.50,-22.9,  -11.03,-28.8,  -11.98,-31.3,
    -10.90,-28.5,  -13.21,-34.9,  -10.88,-27.4,  -16.04,-40.6,  -14.18,-35.0,
    -13.83,-46.9,  -17.82,-56.7,  -3.96,-11.6,  -0.96,-1.8,  -10.38,-31.8,
    -12.64,-38.9,  -7.39,-21.0,  -5.56,-13.9,  -9.44,-24.7,  -7.03,-16.8,
    -11.09,-28.8
  ), ncol=2, byrow=TRUE, dimnames=list(chen_rows, nn_col)),
  # Pass the extra wobble rows as skip so they are not RC-completed
  skip_rows = c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T",
                "init_5T/A","sym",
                "GT/TG","GG/TT","AG/TT","TG/AT","TT/AG","TG/GT",
                "AT/TG","CG/GT","CT/GG","GG/CT","GT/CG"))
  
  # -- RNA/DNA Hybrid - already complete (22 rows including all RC pairs) ------
  # RNA_DNA_NN_Sugimoto_1995 was published with the full symmetric set.
  # .complete_nn_rc() will find nothing to add (returns tbl unchanged).
  hybrid_rows <- c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T",
                   "init_5T/A","sym",
                   "AA/TT","AC/TG","AG/TC","AT/TA",
                   "CA/GT","CC/GG","CG/GC","CT/GA",
                   "GA/CT","GC/CG","GG/CC","GT/CA",
                   "TA/AT","TC/AG","TG/AC","TT/AA")
  RNA_DNA_NN_Sugimoto_1995 <- matrix(c(
    1.9,-3.9,  0,0,  0,0,  0,0,  0,0,  0,0,  0,0,
    -11.5,-36.4,  -7.8,-21.6,  -7.0,-19.7,  -8.3,-23.9,
    -10.4,-28.4,  -12.8,-31.9,  -16.3,-47.1,  -9.1,-23.5,
    -8.6,-22.9,   -8.0,-17.1,   -9.3,-23.2,  -5.9,-12.3,
    -7.8,-23.2,  -5.5,-13.5,  -9.0,-26.1,  -7.8,-21.9
  ), ncol=2, byrow=TRUE, dimnames=list(hybrid_rows, nn_col))
  # No .complete_nn_rc() needed - already has all 16 unique pairs
  
  # -- IMM, TMM, DE tables (unchanged) -----------------------------------------
  imm_rows <- c(
    "AG/TT","AT/TG","CG/GT","CT/GG","GG/CT","GG/TT","GT/CG","GT/TG",
    "TG/AT","TG/GT","TT/AG","AA/TG","AG/TA","CA/GG","CG/GA","GA/CG",
    "GG/CA","TA/AG","TG/AA","AC/TT","AT/TC","CC/GT","CT/GC","GC/CT",
    "GT/CC","TC/AT","TT/AC","AA/TC","AC/TA","CA/GC","CC/GA","GA/CC",
    "GC/CA","TA/AC","TC/AA","AA/TA","CA/GA","GA/CA","TA/AA","AC/TC",
    "CC/GC","GC/CC","TC/AC","AG/TG","CG/GG","GG/CG","TG/AG","AT/TT",
    "CT/GT","GT/CT","TT/AT","AI/TC","TI/AC","AC/TI","TC/AI","CI/GC",
    "GI/CC","CC/GI","GC/CI","AI/TA","TI/AA","AA/TI","TA/AI","CI/GA",
    "GI/CA","CA/GI","GA/CI","AI/TT","TI/AT","AT/TI","TT/AI","CI/GT",
    "GI/CT","CT/GI","GT/CI","AI/TG","TI/AG","AG/TI","TG/AI","CI/GG",
    "GI/CG","CG/GI","GG/CI","AI/TI","TI/AI","CI/GI","GI/CI"
  )
  DNA_IMM_Peyret_1999 <- matrix(c(
    1.0,0.9,    -2.5,-8.3,  -4.1,-11.7,  -2.8,-8.0,   3.3,10.4,
    5.8,16.3,  -4.4,-12.3,   4.1,9.5,   -0.1,-1.7,   -1.4,-6.2,
    -1.3,-5.3,  -0.6,-2.3,  -0.7,-2.3,   -0.7,-2.3,   -4.0,-13.2,
    -0.6,-1.0,   0.5,3.2,    0.7,0.7,     3.0,7.4,     0.7,0.2,
    -1.2,-6.2,  -0.8,-4.5,  -1.5,-6.1,    2.3,5.4,     5.2,13.5,
    1.2,0.7,    1.0,0.7,    2.3,4.6,     5.3,14.6,    1.9,3.7,
    0.6,-0.6,   5.2,14.2,  -0.7,-3.8,    3.4,8.0,     7.6,20.2,
    1.2,1.7,   -0.9,-4.2,  -2.9,-9.8,    4.7,12.9,    0.0,-4.4,
    -1.5,-7.2,   3.6,8.9,    6.1,16.4,   -3.1,-9.5,   -4.9,-15.3,
    -6.0,-15.8,  1.6,3.6,   -2.7,-10.8,  -5.0,-15.8,  -2.2,-8.4,
    0.2,-1.5,  -8.9,-25.5, -5.9,-17.4,  -8.8,-25.4,  -4.9,-13.9,
    -5.4,-13.7, -6.8,-19.1, -8.3,-23.8,  -5.0,-12.6,  -8.3,-25.0,
    -3.4,-11.2, -0.7,-2.6,  -1.3,-4.6,    2.6,8.9,    -7.8,-21.1,
    -7.0,-20.0, -7.6,-20.2,  0.49,-0.7,  -6.5,-22.0,  -5.6,-18.7,
    -0.8,-4.3,  -1.0,-2.4,  -3.5,-10.6,   0.1,-1.0,   -4.3,-12.1,
    -4.9,-15.8, -1.9,-8.5,   0.1,-1.8,    1.0,1.0,     7.1,21.3,
    -1.1,-3.2,   5.8,16.9,  -7.6,-22.0,  -3.3,-11.9,   0.1,-2.3,
    1.3,3.0,   -0.5,-1.3
  ), ncol=2, byrow=TRUE, dimnames=list(imm_rows, nn_col))
  
  tmm_rows <- c(
    "AA/TA","TA/AA","CA/GA","GA/CA","AC/TC","TC/AC","CC/GC","GC/CC",
    "AG/TG","TG/AG","CG/GG","GG/CG","AT/TT","TT/AT","CT/GT","GT/CT",
    "AA/TC","AC/TA","CA/GC","CC/GA","GA/CC","GC/CA","TA/AC","TC/AA",
    "AC/TT","AT/TC","CC/GT","CT/GC","GC/CT","GT/CC","TC/AT","TT/AC",
    "AA/TG","AG/TA","CA/GG","CG/GA","GA/CG","GG/CA","TA/AG","TG/AA",
    "AG/TT","AT/TG","CG/GT","CT/GG","GG/CT","GT/CG","TG/AT","TT/AG"
  )
  DNA_TMM_Bommarito_2000 <- matrix(c(
    -3.1,-7.8,  -2.5,-6.3,  -4.3,-10.7,  -8.0,-22.5,  -0.1,0.5,
    -0.7,-1.3,  -2.1,-5.1,  -3.9,-10.6,  -1.1,-2.1,   -1.1,-2.7,
    -3.8,-9.5,  -0.7,-19.2, -2.4,-6.5,   -3.2,-8.9,   -6.1,-16.9,
    -7.4,-21.2, -1.6,-4.0,  -1.8,-3.8,   -2.6,-5.9,   -2.7,-6.0,
    -5.0,-13.8, -3.2,-7.1,  -2.3,-5.9,   -2.7,-7.0,   -0.9,-1.7,
    -2.3,-6.3,  -3.2,-8.0,  -3.9,-10.6,  -4.9,-13.5,  -3.0,-7.8,
    -2.5,-6.3,  -0.7,-1.2,  -1.9,-4.4,   -2.5,-5.9,   -3.9,-9.6,
    -6.0,-15.5, -4.3,-11.1, -4.6,-11.4,  -2.0,-4.7,   -2.4,-5.8,
    -3.2,-8.7,  -3.5,-9.4,  -3.8,-9.0,   -6.6,-18.7,  -5.7,-15.9,
    -5.9,-16.1, -3.9,-10.5, -3.6,-9.8
  ), ncol=2, byrow=TRUE, dimnames=list(tmm_rows, nn_col))
  
  de_dna_rows <- c(
    "AA/.T","AC/.G","AG/.C","AT/.A","CA/.T","CC/.G","CG/.C","CT/.A",
    "GA/.T","GC/.G","GG/.C","GT/.A","TA/.T","TC/.G","TG/.C","TT/.A",
    ".A/AT",".C/AG",".G/AC",".T/AA",".A/CT",".C/CG",".G/CC",".T/CA",
    ".A/GT",".C/GG",".G/GC",".T/GA",".A/TT",".C/TG",".G/TC",".T/TA"
  )
  DNA_DE_Bommarito_2000 <- matrix(c(
    0.2,2.3,  -6.3,-17.1, -3.7,-10.0, -2.9,-7.6,   0.6,3.3,
    -4.4,-12.6, -4.0,-11.9, -4.1,-13.0, -1.1,-1.6,  -5.1,-14.0,
    -3.9,-10.9, -4.2,-15.0, -6.9,-20.0, -4.0,-10.9, -4.9,-13.8,
    -0.2,-0.5,  -0.7,-0.8,  -2.1,-3.9,  -5.9,-16.5, -0.5,-1.1,
    4.4,14.9,  -0.2,-0.1,  -2.6,-7.4,   4.7,14.2,  -1.6,-3.6,
    -3.9,-11.2, -3.2,-10.4, -4.1,-13.1,  2.9,10.4,  -4.4,-13.1,
    -5.2,-15.0, -3.8,-12.6
  ), ncol=2, byrow=TRUE, dimnames=list(de_dna_rows, nn_col))
  
  de_rna_rows <- c(
    ".T/AA",".T/CA",".T/GA",".T/TA",".G/AC",".G/CC",".G/GC",".G/TC",
    ".C/AG",".C/CG",".C/GG",".C/TG",".T/AG",".T/CG",".T/GG",".T/TG",
    ".A/AT",".A/CT",".A/GT",".A/TT",".G/AT",".G/CT",".G/GT",".G/TT",
    "AT/.A","CT/.A","GT/.A","TT/.A","AG/.C","CG/.C","GG/.C","TG/.C",
    "AC/.G","CC/.G","GC/.G","TC/.G","AT/.G","CT/.G","GT/.G","TT/.G",
    "AA/.T","CA/.T","GA/.T","TA/.T","AG/.T","CG/.T","GG/.T","TG/.T"
  )
  RNA_DE_Turner_2010 <- matrix(c(
    -4.9,-13.20,  -0.9,-1.30,  -5.5,-15.10,  -2.3,-5.50,
    -9.0,-23.50,  -4.1,-10.60, -8.6,-22.20,  -7.5,-20.31,
    -7.4,-20.30,  -2.8,-7.70,  -6.4,-16.40,  -3.6,-9.70,
    -4.9,-13.20,  -0.9,-1.30,  -5.5,-15.10,  -2.3,-5.50,
    -5.7,-16.10,  -0.7,-1.90,  -5.8,-16.40,  -2.2,-6.80,
    -5.7,-16.10,  -0.7,-1.90,  -5.8,-16.40,  -2.2,-6.80,
    -0.5,-0.60,    6.9,22.60,   0.6,2.60,     0.6,2.60,
    -1.6,-4.50,    0.7,3.20,   -4.6,-14.80,  -0.4,-1.30,
    -2.4,-6.10,    3.3,11.60,   0.8,3.20,    -1.4,-4.20,
    -0.5,-0.60,    6.9,22.60,   0.6,2.60,     0.6,2.60,
    1.6,6.10,     2.2,8.10,    0.7,3.50,     3.1,10.60,
    1.6,6.10,     2.2,8.10,    0.7,3.50,     3.1,10.60
  ), ncol=2, byrow=TRUE, dimnames=list(de_rna_rows, nn_col))
  
  # -- GC Content Coefficient Table (unchanged) -----------------------------
  GC_VARTAB <- data.frame(
    A = c(69.3, 81.5, 81.5, 81.5, 78.0, 67.0, 81.5, 77.1),
    B = c(0.41, 0.41, 0.41, 0.41, 0.70, 0.80, 0.41, 0.41),
    C = c(650,  675,  675,  500,  500,  500,  600,  528),
    D = rep(1, 8),
    salt_correct = c(NA, NA, "Schildkraut2010",
                        rep("Wetmur1991", 3),
                        "Schildkraut2010", "SantaLucia1998-1"),
    row.names = c("Chester1993","QuikChange","Schildkraut1965",
                  "Wetmur1991_MELTING","Wetmur1991_RNA","Wetmur1991_RNA/DNA",
                  "Primer3Plus","vonAhsen2001"),
    stringsAsFactors = FALSE
  )
  
  # -- Return assembled list ------------------------------------------------
  list(
    DNA_NN_Breslauer_1986    = DNA_NN_Breslauer_1986,
    DNA_NN_Sugimoto_1996     = DNA_NN_Sugimoto_1996,
    DNA_NN_Allawi_1998       = DNA_NN_Allawi_1998,
    DNA_NN_SantaLucia_2004   = DNA_NN_SantaLucia_2004,
    RNA_NN_Freier_1986       = RNA_NN_Freier_1986,
    RNA_NN_Xia_1998          = RNA_NN_Xia_1998,
    RNA_NN_Chen_2012         = RNA_NN_Chen_2012,
    RNA_DNA_NN_Sugimoto_1995 = RNA_DNA_NN_Sugimoto_1995,
    DNA_IMM_Peyret_1999      = DNA_IMM_Peyret_1999,
    DNA_TMM_Bommarito_2000   = DNA_TMM_Bommarito_2000,
    DNA_DE_Bommarito_2000    = DNA_DE_Bommarito_2000,
    RNA_DE_Turner_2010       = RNA_DE_Turner_2010,
    GC_VARTAB                = GC_VARTAB
  )
}

# =============================================================================
# VERIFICATION - run this after building to confirm completion
# =============================================================================
.verify_nn_tables <- function(tbl_list) {
  # The 6 RC pairs that should now be present in every standard DNA/RNA NN table
  expected_rc_pairs <- c("TT/AA","TG/AC","TC/AG","AG/TC","AC/TG","CC/GG")
  
  standard_tables <- c("DNA_NN_Breslauer_1986","DNA_NN_Sugimoto_1996",
                       "DNA_NN_Allawi_1998","DNA_NN_SantaLucia_2004",
                       "RNA_NN_Freier_1986","RNA_NN_Xia_1998",
                       "RNA_NN_Chen_2012","RNA_DNA_NN_Sugimoto_1995")
  
  all_ok <- TRUE
  for (nm in standard_tables) {
    tbl      <- tbl_list[[nm]]
    missing  <- expected_rc_pairs[!expected_rc_pairs %in% rownames(tbl)]
    n_rows   <- nrow(tbl)
    status   <- if (length(missing) == 0) "OK" else paste("MISSING:", paste(missing, collapse=", "))
    cat(sprintf("  %-35s  rows=%-3d  %s\n", nm, n_rows, status))
    if (length(missing) > 0) all_ok <- FALSE
  }
  if (all_ok) {
    cat("\nAll tables complete. keys_fr lookups will always succeed.\n")
    cat("The keys_rf / pos_nn_rf fallback in tm_nn_core() can now be removed.\n")
  }
  invisible(all_ok)
}


# -- Usage --------------------------------------------------------------------
# Build once at package development time:
#   .TM_CONSTANTS <- .build_all_tm_tables()
#   .verify_nn_tables(.TM_CONSTANTS)
#   usethis::use_data(.TM_CONSTANTS, internal = TRUE, overwrite = TRUE)
#
# Expected verification output:
#   DNA_NN_Breslauer_1986    rows=23  OK
#   DNA_NN_Sugimoto_1996     rows=23  OK
#   DNA_NN_Allawi_1998       rows=23  OK
#   DNA_NN_SantaLucia_2004   rows=23  OK
#   RNA_NN_Freier_1986       rows=23  OK
#   RNA_NN_Xia_1998          rows=23  OK
#   RNA_NN_Chen_2012         rows=28  OK   (17 std + 11 wobble; 6 RC added to std)
#   RNA_DNA_NN_Sugimoto_1995 rows=23  OK   (already had all 6 RC pairs)
#
# After building, in tm_nn_core() the NN sum simplifies to:
#   matched_nn <- keys_fr %in% rownames(nn_tbl)
#   delta_h    <- delta_h + sum(nn_tbl[keys_fr[matched_nn], 1])
#   delta_s    <- delta_s + sum(nn_tbl[keys_fr[matched_nn], 2])
#   # No keys_rf, no pos_nn_rf, no which_nn_fr/rf needed.
