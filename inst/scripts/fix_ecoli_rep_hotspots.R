#!/usr/bin/env Rscript
# fix_ecoli_rep_hotspots.R
#
# Standardise the ecoli_rep_hotspots dataset:
#   1. Every data.frame gets a `chr` column with value "U00096.3"
#   2. Existing NC_000913.3 values are replaced with U00096.3
#   3. `seqnames` columns are renamed to `chr`
#
# Run once from the package root:
#   Rscript inst/scripts/fix_ecoli_rep_hotspots.R

load("data/ecoli_rep_hotspots.rda")

chr_name <- "U00096.3"

for (nm in names(ecoli_rep_hotspots)) {
  df <- ecoli_rep_hotspots[[nm]]

  # rename seqnames -> chr if present
  if ("seqnames" %in% colnames(df)) {
    colnames(df)[colnames(df) == "seqnames"] <- "chr"
    cat(sprintf("  %s: renamed 'seqnames' -> 'chr'\n", nm))
  }

  if ("chr" %in% colnames(df)) {
    # replace any existing values (e.g. NC_000913.3) with U00096.3
    df$chr <- chr_name
    cat(sprintf("  %s: set chr = '%s' (%d rows)\n", nm, chr_name, nrow(df)))
  } else {
    # no chr column at all — add it as the first column
    df <- cbind(chr = chr_name, df)
    cat(sprintf("  %s: added chr = '%s' (%d rows)\n", nm, chr_name, nrow(df)))
  }

  ecoli_rep_hotspots[[nm]] <- df
}

# --- save ---
save(ecoli_rep_hotspots, file = "data/ecoli_rep_hotspots.rda", compress = "xz")
cat("\nSaved updated data/ecoli_rep_hotspots.rda\n\n")

# --- verify ---
for (nm in names(ecoli_rep_hotspots)) {
  df <- ecoli_rep_hotspots[[nm]]
  cat(sprintf("  %s: cols = [%s], chr[1] = '%s'\n",
              nm, paste(colnames(df), collapse = ", "), df$chr[1]))
}
