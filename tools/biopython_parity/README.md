# Biopython parity harness

Cross-checks `TmCalculator::tm_nn()` against
`Bio.SeqUtils.MeltingTemp.Tm_NN()` across every parameter both implementations
share. Development tool, not part of the package; excluded from the build by
`.Rbuildignore`.

```sh
pip install biopython
Rscript tools/biopython_parity/compare_biopython.R            # ~26k cases
Rscript tools/biopython_parity/compare_biopython.R --quick    # smoke test
```

Writes `out/cases.csv`, `out/biopython.csv` and `out/parity_report.csv`, prints
a per-suite summary, and **exits non-zero** if any case that must have matched
did not. Run it after touching `zzz.R`, `tm_nn.R`, `tm_nn_core.cpp` or
`salt_correction.R`.

## Axes covered

| axis | values |
|---|---|
| `nn_table` | the 8 sets Biopython also ships: Breslauer 1986, Sugimoto 1996, Allawi 1998, SantaLucia 2004, Freier 1986, Xia 1998, Chen 2012, Sugimoto 1995 hybrid |
| `salt_method` | all 8, `none` and `saltcorr` 1–7 |
| ionic | Na only; high Na; Na+K+Tris; +Mg; Mg with dNTP chelation; Mg-dominated |
| strand conc. | symmetric 25/25, asymmetric 250/0 and 50/10 |
| `self_comp` | both, the latter on a palindrome |
| duplex | perfect, read from either strand, terminal mismatch at either or both ends, 5' and 3' dangling ends reached by `shift` and by unequal lengths, one and two internal mismatches, and a stack neither table covers |

Sets the package ships but Biopython does not (Weber, Ferreira, Zuber,
Banerjee, Ghosh) are out of scope here; they are covered by
`tests/testthat/test_regressions_1_1_2.R`.

## What counts as a pass

**Suite A (perfect duplexes) and suite C (internal mismatches only) must agree
exactly.** Nothing in either model's treatment of these is supposed to differ,
so a non-exact row is a real defect on one side. These are what the exit status
is based on.

**Suite B (terminal mismatches and dangling ends) is expected to differ**, and
the difference is deliberate. Biopython consumes the terminal mismatch or
dangling end and then indexes all four initiation terms on the *original*
sequence anyway — its own comment says "Now everything 'unusual' at the ends is
handled and removed and we can look at the initiation", but the code that
follows reads `seq`, not the trimmed duplex. TmCalculator indexes them on what
is left, because SantaLucia & Hicks (2004) define the terminal penalty as a
property of the closing base pair, and a mismatch is not a base pair. The
divergence reaches roughly 0.3–0.5 °C for one mismatched end on a 16-mer, and
about 4 °C with `DNA_NN_Breslauer_1986`, whose `init_allA/T` row is the only
one in the package that differs from `init_oneG/C`. Block 15 of
`tests/testthat/test_regressions_1_1_2.R` pins this convention.

**Suite D (uncovered stack): both sides refusing is agreement.** Biopython
raises under `strict=True`; TmCalculator returns `NA` with a warning. One
answering while the other refuses is a failure.

## Files

- `compare_biopython.R` — builds the grid, runs TmCalculator, drives the Python
  worker, joins and classifies
- `biopython_tm.py` — worker: case CSV in, `id,bio_tm,bio_error` out
