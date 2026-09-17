# TmCalculator 1.1.0

## New features

* **`tm_profile()` builds a Tm profile from a sequence source and a set of
  regions.** It tiles the requested regions into windows, retrieves each
  window's sequence and computes its Tm, returning one `GRanges`. The
  source is given by name, either an installed `BSgenome` package or the
  path to a FASTA file, because each task opens the source itself: only a
  record name and a coordinate pair cross between processes, and no
  sequence is ever serialized. That is what makes parallelism worth having
  here, where dividing the sequences of a single `tm_calculate()` call is
  not (see Breaking changes).

  `regions` takes chromosome or record names, numbers, `"name:start-end"`
  strings, a mixture of those, or a `GRanges`. On a BSgenome the `chr`
  prefix is added or removed as the genome requires; FASTA record names are
  matched exactly. Supply a `BPPARAM` to spread the tasks over workers, or
  leave it `NULL` to run them in this process with no dependency on
  `BiocParallel`.

  `window = NULL` gives one window per region, which is the form short
  records call for: a FASTA of array probes, primers or synthetic oligos
  returns one Tm per record.

  A vector of sequences is also accepted, and is staged as a temporary
  FASTA file so that the workers read it rather than receive it. That is
  not a detour: it moves window construction, complement generation and
  result assembly into the worker as well, where sending the sequences
  over a socket would parallelise the inner loop alone. The result is
  keyed by input position so that a sequence dropped for containing `N`
  can still be identified, with the caller's names returned in a `name`
  column. For a handful of sequences the staging and the worker start-up
  cost more than the calculation; `tm_calculate()` remains the right call
  there.

## Bug fix affecting all nearest-neighbor Tm values

* **Four of the six reverse-complement rows added to every nearest-neighbor
  table were transposed, so `tm_nn()` returned incorrect Tm values for any
  sequence containing an AC, AG, TC or TG step.** A key `XY/WZ` denotes the
  duplex 5'-XY-3' / 3'-WZ-5'; read from the other strand the same stack is
  written as the character reversal of the key. `.complete_nn_rc()` instead
  gave `AC/TG` and `TG/AC` each other's parameters, and likewise `AG/TC` and
  `TC/AG`.

  Every table built through that helper was affected, which is all of the
  DNA and RNA sets including `DNA_NN_Breslauer_1986` and
  `DNA_NN_SantaLucia_2004`. `RNA_DNA_NN_Sugimoto_1995` and the RNA/DNA hybrid
  sets ship with all sixteen pairs and were never completed, so they are
  unaffected.

  The size of the error depends on how far apart the transposed rows are in a
  given table and on how often the affected steps occur, so it is
  sequence-dependent rather than a constant offset. It is small for
  SantaLucia 2004 (the transposed pairs differ by 0.1 and 0.4 kcal/mol) and
  considerably larger for Breslauer 1986 (0.7 and 2.2 kcal/mol). **Any Tm
  computed with a previous release should be recomputed.**

  Found by recovering dH and dS from `Tm_NN` in Biopython and from MELTING 5
  and comparing them with this package: those two agree with each other on
  the sequence-dependent part of the sum, and this package did not.
  `tests/testthat/test_nn_rc_completion.R` now re-derives the mapping from
  the reversal rule rather than restating it, and pins the four affected rows
  to their published values.

## Breaking changes

* **`gc()` has been renamed `gc_content()` and is no longer exported under its
  old name.** The exported `gc()` masked `base::gc()` for every user of the
  package, so attaching it printed a masking warning and any subsequent call
  to the garbage collector needed a `base::` prefix, including inside this
  package's own benchmark scripts.

  No deprecated alias is kept, because keeping one would preserve exactly the
  masking the rename is meant to remove. Existing calls fail loudly rather
  than silently: `gc("ACGT")` now reaches `base::gc()`, whose first argument
  is `verbose`, and errors instead of returning a plausible number.

* **`GC` is now a percentage everywhere, computed as
  `100 * (G+C)/(A+C+G+T)`.** `coor_to_genomic_ranges()` previously wrote this
  column as a fraction (0–1) computed over the full window width, so the same
  `GC` column meant different things depending on how the object had been
  built, and windows overlapping assembly gaps appeared GC-poor because N sat
  in the denominator. Objects built by `coor_to_genomic_ranges()` will report
  values 100× larger than before, and slightly larger again wherever N is
  present. Values from `tm_gc()`, `tm_wallace()` and `gc()` are unchanged.

* **The `BPPARAM` argument has been removed from `tm_calculate()`,
  `tm_nn()`, `tm_gc()` and `tm_wallace()`, and `BiocParallel` is no longer
  imported.** With the compiled nearest-neighbor core the per-window loop is
  a minority of a call's runtime; the rest (N filtering, coercion, result
  assembly) runs once in the calling process and cannot be divided. Measured
  on chr1 of GRCh38 (about 1.2 million windows), `SnowParam(5)` never beat
  the serial run: worker start-up and shipping half a gigabyte of sequence
  to the workers cost more than the loop they were dividing. The argument
  therefore offered a slower path and no faster one, and it is gone rather
  than kept as a no-op, so that calls passing it fail loudly.

  Parallelism belongs outside the call, one region per worker, each a
  serial `tm_calculate()`. That pattern needs no support from this package
  and works with any backend; `vignette("hg38_performance_parallel")`
  measures it with `BiocParallel` across a whole genome, which is why
  `BiocParallel` remains in Suggests.

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

* **Fewer S4 operations per call.** Profiling a 100-sequence `tm_nn()` call
  found `validObject()`, `updateObject()` and method dispatch accounting for
  most of the run time, against about 2% in the compiled core. Three sources
  were removed: the two `GRanges` subsets used to drop N-containing regions
  are now taken only when something is actually dropped; the `GC` and `Tm`
  metadata columns are written in a single `mcols<-` assignment instead of
  two `$<-` calls, each of which replaced and revalidated the whole metadata
  table; and the metadata table is extracted once and reused. `tm_gc()`
  receives the same treatment for its two column writes.

  This is a fixed saving per call rather than per sequence, so it does not
  change genome-scale timings. It matters when the functions are called
  repeatedly on short sequences. Results are unchanged.

  Measured on a 100-sequence input of 25 bp oligonucleotides, the two changes
  in this section together reduced the cost of one `tm_calculate()` call from
  25.1 ms to 11.6 ms, a factor of 2.2 (`BPPARAM` default 25.1 -> 21.4 ms; the
  S4 reductions 21.4 -> 11.6 ms). Profiling after the change shows no single
  remaining hotspot: the residual cost is S4 method dispatch distributed
  across the Biostrings, S4Vectors and GenomicRanges accessors, and removing
  it would require keeping the hot path out of S4 entirely.

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

* **`BiocParallel` moved from Imports to Suggests.** It was used only by the
  within-call `BPPARAM` path removed above (see Breaking changes); the
  parallel vignette and the benchmark scripts still use it.

* **`BSgenome` moved from Imports to Suggests, cutting load time by about
  two thirds.** Attaching it pulls in rtracklayer, Rsamtools,
  GenomicAlignments and their dependencies: measured on its own it took 7.2 s
  to load, against 6.5 s for all of TmCalculator and roughly 2 s for the rest
  of the Bioconductor packages combined. Every user paid that whether or not
  they touched a genome.

  Almost nothing used it. `available.genomes` was imported and never called.
  `organism()` and `provider()` are reached only in the second fallback of
  `.resolve_pkg_name()`, when a `BSgenome` object carries no `Package`
  metadata and has a class name shorter than six characters, and they are now
  guarded by `requireNamespace()`. Sequence extraction goes through
  `Biostrings::getSeq()`, whose method for `BSgenome` objects is registered
  when the genome package itself is loaded, which
  `coor_to_genomic_ranges()` already does on demand.

  Genome-wide workflows are unaffected: a genome package such as
  `BSgenome.Ecoli.NCBI.ASM584v2` depends on BSgenome, so it is present
  whenever it is needed.

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
