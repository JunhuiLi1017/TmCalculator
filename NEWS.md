# TmCalculator 1.1.2

## Bug fix affecting every mismatched duplex

Perfect duplexes are **not affected**: all of them return exactly the value
they did in 1.1.1, so the vignettes, the shipped benchmarks and the
cross-tool comparison are unchanged. What changes is `tm_nn()` on a duplex
whose `complement` is not the exact complement of its `sequence`, which is
reached only by supplying `complement_seq` or a non-zero `shift`.

This is a **regression, not an original defect**. From the first CRAN
release in 2022 the stacking loop was an explicit
`if key / else if reversed key / else if other table / else stop()` chain:
it retried the reversal, took one table or the other, and refused to guess at
a stack no table defined. All three properties were lost in a single commit
on 2026-05-26, when that per-base loop was vectorised, and 1.1.0 and 1.1.1
shipped without them. This release restores the first two and records the
third in `ROADMAP.md` item 11. The terminal-key orientation was wrong in both
forms and is corrected here for the first time.

* **A stack and its character reversal are the same stack, and the lookup did
  not know it.** A key `"XY/WZ"` means 5'-XY-3' paired with 3'-WZ-5'; read
  from the other strand the same stack is `"ZW/YX"`, the whole key reversed.
  The published tables store each stack in one orientation only -- of the 87
  keys in `DNA_IMM_Peyret_1999`, 85 have no reversed twin, and neither
  `DNA_TMM_Bommarito_2000` (48) nor `DNA_DE_Bommarito_2000` (32) has any --
  and the lookup was a single hash probe with no retry. A miss contributes
  zero by design, so the omission was silent. Two consequences:

  - **Internal mismatches lost one of their two flanking stacks.** Exactly one
    of the two needs the reversed spelling, so no internal mismatch was ever
    scored completely. The penalty came out about 45% too small on the
    example from issue #10: a T·C mismatch in a 16-mer read as 62.93 °C
    against a perfect 66.60 °C, where the full penalty gives 59.97 °C.

  - **The terminal mismatch table was consulted in the wrong orientation,
    which both missed the mismatches it should have caught and invented ones
    it should not.** Its keys carry the terminal pair *second* (`"AA/TA"` is a
    Watson-Crick pair then a mismatch), while the walk built the terminal key
    with the terminal pair first. A genuine terminal mismatch therefore never
    matched -- 0 times across 400 random duplexes carrying one -- and was
    scored with the internal-mismatch parameters instead. Worse, a duplex
    whose terminal pair *is* Watson-Crick but whose next pair is not produces
    a terminal-first key of exactly the shape the table stores, so it matched
    spuriously and collected a terminal-mismatch penalty it had not earned.
    The terminal keys are now built in the orientation the table uses (the
    reversal of the first stack at the left-hand end, the last stack as-is at
    the right-hand end), and that table is excluded from the reversed retry,
    because for it the orientation carries meaning rather than being two
    spellings of one thing.

  Because the two strands were walked asymmetrically, the same molecule gave
  two different answers depending on which strand was passed as `sequence`.
  Over 400 random 18-24-mers read both ways the two answers differed by up to
  about 3 °C for an internal mismatch and about 4 °C for a terminal one;
  perfect duplexes agreed exactly. That asymmetry is now the regression test
  (`test_regressions_1_1_2.R`): it needs no reference implementation, and it
  would also have caught the table transposition fixed in 1.1.0.

  Reported by \@haraldn in
  [#10](https://github.com/JunhuiLi1017/TmCalculator/issues/10).

* **Dangling ends and the Zuber 2022 end-effect table were already correct**
  and are unchanged. `.right_key()` already rewrites the right-hand terminal
  stack into the orientation those tables use. The end-effect table is
  excluded from the reversed retry on purpose: it lists both orientations of
  most of its keys, so a retry there could return a different stack's value.

* **A stack is now taken from one table, not two.** `nn_table` and
  `imm_table` were consulted independently and both added, so a stack present
  in both received two values. That happens for the G·U wobble stacks of
  `RNA_NN_Chen_2012` and `RNA_NN_Zuber_2022`, whose spellings also occur in
  the default `DNA_IMM_Peyret_1999`, where they mean a DNA G·T mismatch: such
  a stack was charged an RNA wobble parameter plus an unrelated DNA mismatch
  parameter. The nearest-neighbor set now wins, being the one chosen for the
  molecule. No DNA parameter set overlaps the mismatch table in either
  orientation, so **no DNA result moves**; RNA duplexes carrying wobbles move
  by up to about 2.4 °C. (The pre-vectorisation loop already took one table
  or the other; the additive form arrived with the vectorisation.)

* **The 5'-T initiation penalty is charged per strand rather than per
  sequence.** `init_5T/A` is due once for each strand whose 5' end is T: the
  first base of the sequence, and the last base of the complement. Only the
  first was charged, which by itself made Tm depend on which strand was
  passed. The term is zero in all 31 shipped parameter sets, so **no shipped
  result moves**; a user-supplied table with a non-zero value used to break
  the strand symmetry and now does not. (This one is not a regression: the
  pre-vectorisation code had its own version of the fault, testing the first
  base for `A` where it meant the last.)

## Strand direction, which is what invites the mistake

* **`generate_complement()`'s `reverse` argument had its two directions
  written the wrong way round.** It described `reverse = FALSE` as giving a
  sequence "in the same direction (5' to 3')" when that is the plain
  complement, which pairs base for base and therefore reads 3' to 5'; and
  `reverse = TRUE`, the reverse complement, as "3' to 5'" when that is the
  strand written 5' to 3'. Both labels are corrected.

* **`to_genomic_ranges(complement_seq =)` now states the direction, and warns
  when a reverse complement is supplied instead.** The argument wants the
  plain complement, aligned base for base with `input_seq`. The reverse
  complement is what `Biostrings::reverseComplement()` and a supplier's order
  form give you, and passing it pairs every position against the wrong base:
  the duplex is read as almost entirely mismatched and the Tm comes back as a
  large negative number rather than an error. The check is cheap and
  unambiguous -- if reversing the supplied complement makes it pair better, it
  was supplied reversed -- and it only reports, so a genuinely mismatched
  duplex still goes through untouched.

## Bug fix affecting Owczarzy2008 with magnesium

* **`salt_method = "Owczarzy2008"` returned `NA` for both `Tm` and `GC`
  whenever magnesium dominated the monovalent cations.** The method is
  piecewise in `R = sqrt([Mg2+]free) / [Mon]` and the third of its
  three regimes, `R >= 6`, was not implemented: the scalar
  `salt_correct()` reached the end of its `if`/`else` without assigning a
  result, and the vectorized path returned `NA` to reproduce that. Ordinary
  PCR-like conditions land there, so a whole class of calls silently produced
  nothing:

  ```r
  tm_nn(to_genomic_ranges("GCATCGTAGGCTAGCTTGCA"),
        salt_method = "Owczarzy2008", Na = 1, Mg = 5)
  #  before: Tm NA, GC NA      now: Tm 63.16, GC 55
  ```

  The missing branch is the published expression with its constants
  unmodified (a = 3.92, b = -0.911, c = 6.26, d = 1.42, e = -48.2, f = 52.5,
  g = 8.31); only the competing regime reparameterises them on `[Mon]`.

* **Owczarzy2008 with magnesium but no monovalent cation applied no
  correction at all.** `Na = 0, K = 0, Tris = 0, Mg = 5` fell into the guard
  that returns zero when the monovalent concentration is zero, which is right
  for the six methods that take the logarithm of it and wrong for this one:
  in the divalent-dominated regime `[Mon]` drops out of the expression
  and the correction is defined. This was worse than the `NA` above, because
  nothing marked the result as untrustworthy. Such calls now return the
  divalent form and the value changes.

## Breaking change in `tm_gc()`

* **`salt_method = "Owczarzy2004"` and `"Owczarzy2008"` are no longer accepted
  by `tm_gc()`, and `tm_calculate(method = "tm_gc", ...)` rejects them too.**
  They are corrections to the reciprocal of the melting temperature in
  kelvin, referenced to the same duplex in 1 M Na+, and they carry a
  `1/(2(N-1))` duplex-length term of their own. The GC-content formulas are
  on neither footing and already have a length term, so the two cannot be
  combined even with the reciprocal arithmetic done correctly. Both remain
  available in `tm_nn()`, which is where they belong.

  They were reachable only through `userset`, where `tm_gc()` added them to
  the Tm additively rather than reciprocally and so applied a shift of order
  1e-5 degrees instead of the intended correction. Any such call was
  therefore already producing an essentially uncorrected Tm.

* **`tm_gc(salt_method =)` now defaults to `NULL` and says when it is
  ignored.** With a built-in `variant` the salt term is part of the published
  formula, and naming a different one was silently overridden. It is still
  overridden, because the formula has to be the one it is labelled as, but
  the call now warns. `NULL` selects the variant's own correction, or
  `"Schildkraut2010"` with `userset`, which is what the function has always
  applied in that case. `NA` and `"none"` drop the correction altogether on
  either path: that is not a substitution, so it is honoured rather than
  warned about, and `$options` reports the result as uncorrected. The
  documentation had promised `NA` since 1.1.0 without the code supporting it.

* **`tm_calculate()` reports the salt correction that was applied.** On the
  profiling route `$options$salt_method` echoed the requested value, so a
  `method = "tm_gc"` run could say `"Schildkraut2010"` while `vonAhsen2001`'s
  own `"SantaLucia1998-1"` term had been used, or name a correction for
  `Chester1993`, which has none. The direct route was already correct.

## Behaviour changes

* **`GC` no longer follows `Tm` into `NA`.** A sequence the thermodynamic
  model cannot evaluate gets `NA` for `Tm` alone; its base composition is a
  property of the sequence, not of the model, and is still reported. `GC` is
  `NA` only where the sequence itself has no countable base, which is what
  `gc_content()` returns for the same input.

* **`tm_nn()` now warns when it returns `NA`.** It used to do so silently, so
  a caller who did not test for `NA` carried it into a mean or a plot without
  ever being told. The two causes are reported separately, because they call
  for different fixes: a sequence the model cannot evaluate is an input
  problem, an undefined salt correction is a condition problem.

# TmCalculator 1.1.1

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

* **`tm_profile()` has been removed**, having been merged into
  `tm_calculate()`, which takes the same arguments under the same names. The
  one difference to watch for is the return value: `tm_profile()` gave a bare
  `GRanges`, whereas `tm_calculate()` gives a `TmCalculator` object, so the
  profile is `$gr`.

* **An unnamed sequence is keyed by its position, not by `chr1`.**
  `tm_calculate(c("ACGT...", "GGCC..."))` labelled every window `chr1`,
  naming a chromosome that was not involved and giving every row the same
  key, so nothing in the result said which input a given Tm came from. The
  `seqname` is now `"1"`, `"2"` and so on, which is also what
  `regions = "3:1-20"` already meant for an unnamed vector. Sequences that
  carry names, or names in `chr:start-end` form, are unaffected.

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
