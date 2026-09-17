# TmCalculator 1.1.0

## New features

* **`tm_calculate()` now takes a source and a set of regions, and can spread
  the work over processes.** It accepted four kinds of input already; it can
  now be told which part of that input to use, how finely to tile it, and how
  many workers to divide it among:

  ```r
  tm_calculate("BSgenome.Hsapiens.UCSC.hg38", window = 200, slide = 200,
               BPPARAM = SnowParam(workers = 5))          # a whole genome
  tm_calculate("contigs.fa.gz", regions = c("contig_7", "contig_9:1-50000"))
  tm_calculate(oligos)                                    # unchanged
  ```

  `regions` means the same thing for every source, because the identifier
  before the colon is resolved against whatever names the source itself
  offers and falls back to position. So `"chr1"` is a chromosome in a
  BSgenome, a record in a FASTA file and a `seqname` in a `GRanges`, and
  `"1:1-200"` is the first 200 bases of the first sequence in an unnamed
  character vector.

  Parallelism divides the work **by region and never by the sequences of one
  region**. Each task opens the source itself, so only a name and a
  coordinate pair cross between processes; sequences supplied directly are
  staged to a temporary FASTA for the same reason. That is the arrangement
  that pays: dividing the sequences of a single call parallelises the inner
  loop alone and leaves window construction, sequence retrieval and result
  assembly in the calling process, which measured slower than not dividing
  them at all.

  A call that passes sequences and nothing else behaves exactly as before,
  and takes the same short path through the function.

* **`tm_profile()` is deprecated**, having been merged into
  `tm_calculate()`. It still works, warns once per session, and returns a
  bare `GRanges` as it always did; `tm_calculate()` returns a `TmCalculator`
  object, so the profile is `$gr`.

* **FASTA input keeps its record names.** `tm_calculate("reads.fa")`
  previously labelled every window `chr1`, because a record name that is not
  in `chr:start-end` form fell through to the default; the record name now
  becomes the `seqname`. Tm values are unaffected.

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

* **`BPPARAM` no longer divides the sequences of one call, and is gone from
  `tm_nn()`, `tm_gc()` and `tm_wallace()` entirely.** `tm_calculate()` keeps
  a `BPPARAM`, but it now means something else: the work is divided by
  region, with each task opening the source for itself, and never by
  splitting the sequences of one region among workers. With the compiled nearest-neighbor core the per-window loop is
  a minority of a call's runtime; the rest (N filtering, coercion, result
  assembly) runs once in the calling process and cannot be divided. Measured
  on chr1 of GRCh38 (about 1.2 million windows), `SnowParam(5)` never beat
  the serial run: worker start-up and shipping half a gigabyte of sequence
  to the workers cost more than the loop they were dividing. The argument
  therefore offered a slower path and no faster one, and it is gone rather
  than kept as a no-op, so that calls passing it fail loudly.

  Parallelism is now the business of `tm_calculate()` itself, one region per
  worker, which is why `BiocParallel` stays in Imports rather than moving to
  Suggests: dividing the work is part of what the function does, not an
  optional extra a user assembles around it.
  `vignette("hg38_performance_parallel")` measures it across a whole genome.

* **`tm_nn()` now reports GC on the same definition it already used for salt
  correction.** It previously reported `(G+C)/length` while correcting with
  `(G+C)/(A+C+G+T)`. The two differ only when inosine is present, since `I`
  counts in the length but is not a determinable base; sequences containing N
  are skipped before this point and so were never affected. Tm values are
  unchanged in all cases.

## Performance

* **`coor_to_genomic_ranges(method = "preload_chr")` now loads the span its
  windows cover, not the chromosome they sit on.** The preload path called
  `genome[[chr]]`, which decompresses the whole chromosome whatever was
  asked for, so a request for 200 kb of chr21 read all 46.7 Mb of it. That
  cost was paid once per task, and a segmented genome-scale run is many
  tasks per chromosome: five 50 Mb segments of chr1 each decompressed its
  full 249 Mb, five times over, to cover 249 Mb once. The span is now read
  in one `getSeq()` call and window coordinates are shifted onto it.

  Windows scattered along a whole chromosome still span it, so the case the
  preload path was written for is unchanged. A dense run over part of a
  chromosome, which is what `regions` and `segment_size` produce, reads what
  it uses. Sequences and Tm values are identical either way;
  `inst/scripts/test_tm_calculate_merged.R` pins them against `getSeq()` at
  the same coordinates rather than against another run of this code, since a
  mistake in the shift would move the whole profile consistently and agree
  with itself.


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

* **`BiocParallel` stays in Imports.** It briefly moved to Suggests when the
  within-call `BPPARAM` path was removed (see Breaking changes), on the
  reasoning that nothing in the package used it any more. Merging
  `tm_profile()` into `tm_calculate()` put it back: `.tm_run()` dispatches
  the tasks of a genome-scale call with `bplapply()`, so dividing the work is
  something the package does rather than something a user assembles around
  it, and the dependency is not optional.

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
