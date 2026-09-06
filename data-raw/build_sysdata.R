#!/usr/bin/env Rscript
# ===========================================================================
# build_sysdata.R -- regenerate R/sysdata.rda from R/zzz.R
#
# .TM_CONSTANTS (every nearest-neighbor and GC parameter table) is built once
# at development time and stored in R/sysdata.rda, not rebuilt when the
# package loads. Editing .build_all_tm_tables() in R/zzz.R therefore has NO
# EFFECT until this script is run: the installed package keeps serving the
# tables that were baked in when sysdata.rda was last written.
#
# Run from the package root:
#   Rscript data-raw/build_sysdata.R
#
# The script prints every table row whose value changed, so that the effect
# of an edit is visible rather than silent. Any other object already in
# sysdata.rda is preserved.
# ===========================================================================

stopifnot(file.exists("R/zzz.R"), file.exists("DESCRIPTION"))

## -- 1. Keep whatever else lives in sysdata.rda ----------------------------
old <- new.env(parent = emptyenv())
if (file.exists("R/sysdata.rda")) {
  load("R/sysdata.rda", envir = old)
  message("objects currently in R/sysdata.rda: ",
          paste(ls(old, all.names = TRUE), collapse = ", "))
}
previous <- if (exists(".TM_CONSTANTS", envir = old, inherits = FALSE))
  get(".TM_CONSTANTS", envir = old) else NULL

## -- 2. Rebuild from source ------------------------------------------------
# sys.source() defines the builders without attaching the package, so the
# tables come from the file on disk rather than from an installed copy.
env <- new.env(parent = globalenv())
sys.source("R/zzz.R", envir = env)
fresh <- env$.build_all_tm_tables()

if (is.function(env$.verify_nn_tables)) {
  message("\n-- .verify_nn_tables() --")
  env$.verify_nn_tables(fresh)
}

## -- 3. Independent check of the reverse-complement rows -------------------
# .verify_nn_tables() only counts rows, which is why a transposition in the
# completion map survived it. A key "XY/WZ" and its character reversal are
# the same duplex read from opposite strands, so the two rows must carry
# identical values wherever both are present.
rev_key <- function(k)
  vapply(strsplit(k, "", fixed = TRUE),
         function(ch) paste(rev(ch), collapse = ""), character(1))

bad <- character(0)
for (nm in names(fresh)) {
  tbl <- fresh[[nm]]
  if (is.null(dim(tbl)) || is.null(rownames(tbl))) next
  keys <- grep("^[ACGU T]{2}/[ACGU T]{2}$", rownames(tbl), value = TRUE)
  for (k in keys) {
    src <- rev_key(k)
    if (src %in% rownames(tbl) && !isTRUE(all.equal(unname(tbl[k, ]),
                                                    unname(tbl[src, ]))))
      bad <- c(bad, sprintf("%s: %s != %s", nm, k, src))
  }
}
if (length(bad)) {
  message("\n-- reverse-complement rows that DISAGREE --")
  message(paste(bad, collapse = "\n"))
  message("\nSome parameter sets are fitted without imposing strand symmetry, ",
          "so a disagreement is not automatically an error. Confirm against ",
          "the source .par file or publication before accepting it.")
} else {
  message("\nAll reverse-complement rows agree with their reversed key.")
}

## -- 4. Report what changed ------------------------------------------------
if (!is.null(previous)) {
  message("\n-- rows whose values changed --")
  n_changed <- 0L
  for (nm in union(names(previous), names(fresh))) {
    a <- previous[[nm]]; b <- fresh[[nm]]
    if (is.null(a) || is.null(b) || is.null(dim(a)) || is.null(dim(b))) next
    shared <- intersect(rownames(a), rownames(b))
    for (r in shared) {
      if (!isTRUE(all.equal(unname(a[r, ]), unname(b[r, ])))) {
        message(sprintf("  %-26s %-8s  %s  ->  %s", nm, r,
                        paste(format(unname(a[r, ]), width = 8), collapse = " "),
                        paste(format(unname(b[r, ]), width = 8), collapse = " ")))
        n_changed <- n_changed + 1L
      }
    }
    for (r in setdiff(rownames(b), rownames(a)))
      message(sprintf("  %-26s %-8s  ADDED", nm, r))
    for (r in setdiff(rownames(a), rownames(b)))
      message(sprintf("  %-26s %-8s  REMOVED", nm, r))
  }
  if (n_changed == 0L) message("  (none)")
  message("\n", n_changed, " row(s) changed. Every Tm computed with the ",
          "previous sysdata.rda is affected.")
}

## -- 5. Write --------------------------------------------------------------
assign(".TM_CONSTANTS", fresh, envir = old)
save(list = ls(old, all.names = TRUE), envir = old,
     file = "R/sysdata.rda", version = 2, compress = "xz")
message("\nWritten R/sysdata.rda (", length(fresh), " tables). ",
        "Reinstall the package before running the tests.")
