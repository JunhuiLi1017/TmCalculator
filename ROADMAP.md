# TmCalculator roadmap

Candidate items for the next feature release (post-1.1.0), roughly in
order of expected user-visible impact. Background and measurements are in
`vignettes/hg38_performance_parallel.Rmd`; current baseline on the
6-core/16 GB reference machine: chr1 pipeline ~52 s, whole genome 295 s
with 5 workers.

## Performance

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

3. **Shared-memory threading in the core (RcppParallel).** The compiled
   Tm loop is embarrassingly parallel and thread-safe (read-only hash
   tables, per-sequence stack state). RcppParallel threads would cut the
   ~9 s loop to ~2 s with zero serialization and no BPPARAM needed,
   making a single process saturate all cores. Avoid OpenMP (absent
   from Apple clang); RcppParallel/TBB is the portable route. With items
   1-3 combined, chr1 should drop from ~52 s to ~15 s single-process.

## Thermodynamic models

6. **Zuber et al. (2022) RNA end effects** (NAR 50:5251,
   doi:10.1093/nar/gkac261; raised by a manuscript reviewer). Their
   improved terminal terms depend on the penultimate base pair ("AU end
   on AU" vs "AU end on CG", plus GU-end variants), which the current
   core's initiation system cannot express. Requires a small core
   extension: an optional end-stack table looked up on the terminal
   dinucleotide keys and added without trimming (unlike tmm). The
   companion Banerjee et al. (2020) hybrid set (gkaa572) is already
   implemented as RNA_DNA_NN_Banerjee_2020 in 1.10.0.

## API hygiene

4. **Rename `gc()` to `gc_content()`.** The exported `gc()` masks
   `base::gc()` for every user of the package. Keep `gc()` as a
   soft-deprecated alias (`.Deprecated`) for one or two releases.

5. **Unify GC units.** `coor_to_genomic_ranges()` writes the `GC`
   metadata column as a fraction (0-1); `tm_nn`/`tm_gc`/`tm_wallace`
   compute it as a percentage (0-100). The same column name currently
   means different units depending on the input path. Pick one (percent
   matches the Tm literature), convert on ingestion, note in NEWS.

## Done in 1.1.0 (for reference)

- BPPARAM (BiocParallel) in tm_calculate/tm_nn/tm_gc/tm_wallace.
- Rcpp nearest-neighbor core (~10x serial speedup); sequence cleaning
  moved into the core.
- Lazy `result$df` via `$.TmCalculator`.
- `preload_chr` extraction rewritten with `extractAt()`; N-end scan made
  in-memory (chr21's 5 Mb leading gap: 30 s -> ~1 s).
- hg38 performance vignette + `inst/scripts/benchmark_hg38.R`.
