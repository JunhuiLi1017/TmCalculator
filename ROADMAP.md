# TmCalculator roadmap

Candidate items for the next feature release (post-1.1.0), roughly in
order of expected user-visible impact. Background and measurements are in
`vignettes/hg38_performance_parallel.Rmd`; current baseline on the
6-core/16 GB reference machine: chr1 pipeline ~52 s single-process.

Whole-genome figures are deliberately omitted here until
`inst/scripts/bench_parallel_cluster.R` has been run: repeated runs of the
same worker count on the reference laptop differed by nearly 30% (295 s and
379 s for five workers), so no single number from that sweep is quotable.

## Performance

**Items 1-3 form one chain and must be done in order.** The goal is the
*data path*, not compiling more arithmetic: of the ~52 s chr1 pipeline,
window generation is ~5.5 s, sequence retrieval ~15 s and Tm computation
~31 s, and of that 31 s only ~9 s is the arithmetic. The rest is moving
sequence data into a form the core can read. Item 4 is independent of the
chain and can be done at any time.

(An earlier revision of this file asserted that `tm_gc` and `tm_wallace`
were already whole-vector operations and not worth compiling. That was
assumed, not checked, and it is wrong; see item 4.)

1. **Direct DNAStringSet access in the C++ core.** `tm_nn` currently
   coerces sequence and complement columns via `as.character()` before
   calling `cpp_tm_nn_dhds()` — the largest single item in the ~21 s of
   serial preprocessing per chr1. Reading XStringSet raw bytes directly
   in C++ (`LinkingTo: S4Vectors, IRanges, XVector, Biostrings`, as
   DECIPHER does) eliminates the coercion and lets C++ also do the
   N-scan. Expected: preprocessing 21 s -> under 5 s.

2. **Lazy complement.** `coor_to_genomic_ranges()` materializes a
   complement column (~250 MB per chr1, one full pass) that the NN core
   could derive per base on the fly for perfect duplexes. Make the
   column optional: generate only when the user supplies
   `complement_seq` or mismatch analysis is requested. Saves ~5-8 s and
   ~250 MB per chr1. API note: `complement` metadata column becomes
   conditional — document in NEWS.

3. **Shared-memory threading in the core (RcppParallel).** *Blocked by
   item 1 — this is a hard prerequisite, not a preference.* RcppParallel
   worker threads may not touch R objects or call the R API, and
   `cpp_tm_nn_dhds()` currently takes a `CharacterVector` and converts
   each element with `Rcpp::as<std::string>` inside the loop. Until the
   hot loop reads raw bytes and holds no R object, it cannot legally run
   off the main thread.

   Once unblocked: the Tm loop is embarrassingly parallel and already
   thread-safe (read-only hash tables, per-sequence stack state), so
   threads should cut the ~9 s loop to ~2 s with zero serialization —
   no 0.5 GB of sequence copied per worker, and none of the
   copy-on-write page duplication that makes `MulticoreParam()`
   counterproductive today. A single process would then saturate all
   cores, which is what makes `BPPARAM` worth supplying for a single
   chromosome rather than only across chromosomes. Avoid OpenMP (absent
   from Apple clang); RcppParallel/TBB is the portable route.

   With items 1-3 combined, chr1 should drop from ~52 s to ~15 s
   single-process. Note where that leaves the ceiling: sequence
   retrieval (~15 s, inside BSgenome/Biostrings rather than this
   package) would then be essentially the entire remaining cost, so
   further gains require attacking I/O, not arithmetic. Do not promise
   more than this without measuring.

4. **Base counting outside the per-sequence loop.** *`tm_gc` done;
   `tm_wallace` outstanding.*

   `tm_gc` and `tm_wallace` looped over sequences in R calling `gc()`, which
   split each sequence with `s2c()` and scanned it five times, while
   `salt_correct()` did the same again. On the E. coli case study (23,208
   windows of 200 bp, installed build) that made the composition-based
   method **77x slower** than the compiled nearest-neighbor path, which at
   hg38 scale would have meant roughly 9 hours against 7 minutes.

   | method  | before  | after   | per window (after) |
   |---------|---------|---------|--------------------|
   | `tm_nn` | 0.668 s | 0.677 s | 29 us               |
   | `tm_gc` | 51.72 s | 1.754 s | 76 us               |

   Fixed by adding `.gc_vec()` (whole-vector GC via `gsub(fixed = TRUE)`,
   mirroring `gc()` including ambiguity apportioning) and rewriting
   `.tm_gc_chunk()` to compute GC, lengths and salt correction for the whole
   chunk at once, reusing the existing `.salt_correct_vec()`. Equivalence is
   covered by `tests/testthat/test_gc_vec.R`. `Biostrings::letterFrequency()`
   was not used because `gc()` accepts `I`, which `DNAStringSet()` rejects.

   Remaining:
   - `tm_wallace` still loops and still calls `s2c()` twice per sequence.
     Same fix, lower stakes because it warns above 30 bp.
   - The counting was subsequently moved into C++ (`src/gc_core.cpp`,
     `cpp_base_counts()`), which reads the string bytes directly through
     `CHAR(STRING_ELT(...))` and tallies all IUPAC codes in one pass per
     sequence. The apportioning arithmetic stays in R so the published
     formula remains the readable one. Re-measure to confirm the gain over
     the `gsub` version before quoting a number.
   - `mismatch = TRUE` in `.tm_gc_chunk()` computes `seqs %in% "X"`, which is
     TRUE only when an entire sequence is the single character `"X"`, so the
     mismatch penalty is inert. Preserved verbatim during vectorisation so
     that no published value changed. Fix separately and note in NEWS,
     since correcting it *will* alter output.
   - `gc()` branches on `length(input_seq) == 1` and treats a longer vector
     as an already-split character vector, so passing several sequences at
     once returns a number computed from the wrong thing rather than raising
     an error. Make the contract explicit when renaming to `gc_content()`
     (item 6).

5. **Guard against oversubscription once threading lands.** With
   in-core threads *and* `BPPARAM`, `SnowParam(5)` running an 8-thread
   core requests 40 cores. Decide and document the composition rule —
   most likely: threads default to 1 when a multi-worker BPPARAM is
   supplied, overridable explicitly — and add a test that the two
   mechanisms never multiply silently.

## API hygiene

6. **Rename `gc()` to `gc_content()`.** The exported `gc()` masks
   `base::gc()` for every user of the package. Keep `gc()` as a
   soft-deprecated alias (`.Deprecated`) for one or two releases.

7. **Unify what `GC` means.** Three separate issues share this heading; only
   the first is clear-cut.

   **(a) Units, unambiguous bug.** `coor_to_genomic_ranges()` writes the `GC`
   metadata column as a fraction (0-1) while `tm_nn`/`tm_gc`/`tm_wallace`
   write a percentage (0-100). Same column name, different scale depending
   on how the object was built. Pick percent (it matches the Tm literature),
   convert on ingestion, note in NEWS.

   **(b) Definition, narrow in practice.** `tm_nn` reports
   `100 * nGC / len` while `salt_correct()` and `tm_gc` use
   `100 * nGC / (A+C+G+T)`; `tm_nn` computes both internally
   (`R/tm_nn.R:591-594`) and reports one while correcting with the other.
   An earlier version of this note claimed the two diverge on N-containing
   windows. They do not: `tm_nn` drops any sequence containing N at
   `R/tm_nn.R:416` before the calculation, and cleaning reduces the rest to
   A/C/G/T/I, so ambiguous codes never reach it either. The only reachable
   divergence is **inosine**, which counts in `len` but not in A+C+G+T.
   Real, but minor; do not over-engineer it.

   **(c) N handling differs between methods, and this is the user-visible
   one.** `tm_nn` skips windows containing N and warns;
   `tm_gc`/`tm_wallace` keep them and compute GC over the non-N bases:

   ```r
   s <- c("ACGTACGTNN", "ACGTACGTAC")
   tm_calculate(to_genomic_ranges(s), method = "tm_nn")$gr$GC   # 50      (1 row)
   tm_calculate(to_genomic_ranges(s), method = "tm_gc")$gr$GC   # 50 50   (2 rows)
   ```

   Changing method therefore changes which windows exist, not merely their
   values. At hg38 scale that means method comparisons silently operate on
   different window sets. Decide whether the composition-based methods
   should adopt the same skip rule, or whether `tm_nn` should return NA
   rather than dropping rows; either is defensible, but the current
   asymmetry is not documented anywhere.

   Whatever is decided for (a) and (b), route every path through
   `cpp_base_counts()`, which already returns the counts needed for any of
   these conventions, and delete `.GC_fast()` (now reachable only from the R
   reference implementation at `R/tm_nn.R:882`).

## Documentation and interoperability

8. **Verify and document `BatchtoolsParam()` on an HPC scheduler.**
   `.bp_map_chunks()` calls `BiocParallel::bplapply()` generically, so any
   `BiocParallelParam` is accepted by construction — including
   `DoparParam()` (bridging to any registered foreach backend, and through
   `doFuture` to the future framework) and `BatchtoolsParam()` (SLURM, LSF,
   SGE, Torque). None of this is tested: `tests/testthat/test_tm_calculate.R`
   covers only `SerialParam` and `SnowParam`, so cluster support is currently
   an inference from the generic call, not a verified feature, and is
   deliberately not claimed in the manuscript.

   Concretely: while running `inst/scripts/bench_parallel_cluster.R` on the
   LSF cluster, try `BatchtoolsParam(cluster = "lsf")` once. If it works,
   document it in `?tm_calculate` and add a skipped-unless-available test.
   The practical value is that users could dispatch to a scheduler from R
   without writing submission scripts. Note the likely obstacle before
   spending time on it: batchtools workers are separate R sessions on other
   nodes, so the package and the BSgenome must be installed there, and the
   within-call chunking would ship sequence data across the network — which
   is the pattern already shown to be slower than region-level parallelism.
   Expect it to be useful for dispatching whole regions, not for `BPPARAM`
   inside a single call.

## Done in 1.1.0 (for reference)

- BPPARAM (BiocParallel) in tm_calculate/tm_nn/tm_gc/tm_wallace.
- Rcpp nearest-neighbor core (~10x serial speedup); sequence cleaning
  moved into the core.
- Lazy `result$df` via `$.TmCalculator`.
- `preload_chr` extraction rewritten with `extractAt()`; N-end scan made
  in-memory (chr21's 5 Mb leading gap: 30 s -> ~1 s).
- hg38 performance vignette + `inst/scripts/benchmark_hg38.R`.
