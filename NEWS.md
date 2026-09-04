# TmCalculator 1.10.0

## Breaking changes

* **`GC` is now a percentage everywhere, computed as
  `100 * (G+C)/(A+C+G+T)`.** `coor_to_genomic_ranges()` previously wrote this
  column as a fraction (0–1) computed over the full window width, so the same
  `GC` column meant different things depending on how the object had been
  built, and windows overlapping assembly gaps appeared GC-poor because N sat
  in the denominator. Objects built by `coor_to_genomic_ranges()` will report
  values 100× larger than before, and slightly larger again wherever N is
  present. Values from `tm_gc()`, `tm_wallace()` and `gc()` are unchanged.

* **`tm_nn()` now reports GC on the same definition it already used for salt
  correction.** It previously reported `(G+C)/length` while correcting with
  `(G+C)/(A+C+G+T)`. The two differ only when inosine is present, since `I`
  counts in the length but is not a determinable base; sequences containing N
  are skipped before this point and so were never affected. Tm values are
  unchanged in all cases.

## Performance

* **`tm_gc()` is roughly 260× faster.** On the *E. coli* case study (23,208
  windows of 200 bp) it fell from 51.7 s to 0.198 s, and is now faster than
  `tm_nn()` rather than 77× slower. `tm_wallace()` receives the same fix.
  The cause was that both looped over sequences in R calling `gc()`, which
  split each sequence into a character vector with `seqinr::s2c()` and scanned
  it five times, while `salt_correct()` repeated the same work.

* Base counting now happens once per sequence in compiled code
  (`cpp_base_counts()`), reading string bytes directly. `gc()`, `.gc_vec()`,
  `tm_gc()`, `tm_wallace()`, `tm_nn()` and `salt_correct()` all share this one
  implementation and one definition of GC; the internal `.GC_fast()`, which
  carried a second definition, has been removed.

## Dependencies

* **`seqinr` is no longer required.** It was used for `s2c()`/`c2s()` in
  `gc()`, `tm_gc()`, `tm_wallace()`, `salt_correct()` and
  `generate_complement()`, and for `read.fasta()` in `fa_to_genomic_ranges()`.
  The first group is replaced by the compiled base counter, the second by
  `Biostrings::readBStringSet()`, and `generate_complement()` now uses base R
  `chartr()`. `readBStringSet()` rather than `readDNAStringSet()` is used
  deliberately, since the latter validates against the DNA alphabet and would
  reject the RNA input this package supports. FASTA parsing behaviour is
  unchanged: no alphabet restriction, case preserved, sequences named by the
  first word of the header.

* `rlang` removed from Suggests; it was referenced nowhere.

* `BSgenome.Hsapiens.UCSC.hg38` removed from Suggests. Its examples are inside
  `\dontrun{}` and the hg38 vignette sets `eval = FALSE`, so `R CMD check`
  never loads it, but listing it obliged every check environment to download
  roughly 850 MB. The requirement is still stated in the vignette text.

## Bug fixes

* `gc()` given a character vector of length > 1 treated the elements as
  individual bases, so passing several complete sequences silently returned a
  value computed from the wrong thing. Such input is now concatenated into one
  sequence, consistent with the documented `gc(c("a","t","g","c"))` form.

## Known issues

* `tm_nn()` skips windows containing N and warns, whereas `tm_gc()` and
  `tm_wallace()` retain them and compute GC over the remaining bases, so
  changing method can change which windows are returned rather than only their
  values. See `ROADMAP.md` item 7(c); this is unchanged in this release.

* With `mismatch = TRUE`, `tm_gc()` evaluates `sequence %in% "X"`, which is
  true only when an entire sequence is the single character `"X"`, so the
  mismatch penalty is inert. Preserved verbatim during the performance work so
  that no value changed; see `ROADMAP.md` item 4.
