# TmCalculator <a href="https://github.com/JunhuiLi1017/TmCalculator"><img src="man/figures/logo.png" align="right" height="138" /></a>
v1.0.9

Genome-wide nucleic acid melting temperature (Tm) profiling and multi-omics
integration. Results are returned as `GRanges` objects, so Tm can be used
directly as a quantitative genomic feature alongside ATAC-seq, RNA-seq,
ChIP-seq and other assays.

## 1. install

  `install.packages("TmCalculator")`

  install dev version from github

  `pak::pkg_install("JunhuiLi1017/TmCalculator@dev")`

## 2. usage and examples

Please see the vignetts for the details.

```r
library(TmCalculator)

seqs <- to_genomic_ranges("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
tm_calculate(seqs, method = "tm_nn", nn_table = "DNA_NN_SantaLucia_2004", Na = 50)
```

## 3. thermodynamic parameter sets

Twenty-seven nearest-neighbor parameter sets are available, in two families.

**Reference-salt sets** were fitted at a single reference sodium concentration.
Other conditions are reached through the `salt_method` correction formulas.

| Duplex  | Sets |
|---------|------|
| DNA/DNA | `DNA_NN_Breslauer_1986`, `DNA_NN_Sugimoto_1996`, `DNA_NN_Allawi_1998`, `DNA_NN_SantaLucia_2004` (default) |
| RNA/RNA | `RNA_NN_Freier_1986`, `RNA_NN_Xia_1998`, `RNA_NN_Chen_2012` |
| RNA/DNA | `RNA_DNA_NN_Sugimoto_1995` |

**Condition-specific sets** were fitted directly at the sodium concentration
shown, by melting-temperature optimization. They are intended to *replace*
salt correction rather than be corrected. When the requested `Na` matches the
concentration a set was fitted at, salt correction is skipped automatically;
when it does not, the correction is applied with a warning.

| Duplex  | Sets | Fitted at |
|---------|------|-----------|
| DNA/DNA | `DNA_NN_Weber_2015` | 1020 mM |
| DNA/DNA | `DNA_NN_Weber_OW04_69` / `_119` / `_220` / `_621` / `_1020` | 69–1020 mM |
| RNA/RNA | `RNA_NN_Weber_VIF_71` / `_121` / `_221` / `_621` / `_1021` | 71–1021 mM |
| RNA/RNA | `RNA_NN_Weber_FIF_71` / `_121` / `_221` / `_621` / `_1021` | 71–1021 mM |
| RNA/DNA | `RNA_DNA_NN_Weber_2019_FT`, `RNA_DNA_NN_Weber_2019_VH` | 1000 mM |
| RNA/DNA | `RNA_DNA_NN_Weber_2019_LS` | 100 mM |

For RNA, the VIF (variable initiation factors) sets gave better cross-validation
than FIF. For RNA/DNA hybrids at high salt, `..._FT` was the best-performing set
in the source study.

```r
# Fitted at 100 mM, so no salt correction is applied on top of it
res <- tm_calculate(seqs, method = "tm_nn",
                    nn_table = "RNA_DNA_NN_Weber_2019_LS", Na = 100)
res$options[["Salt correction applied"]]                    # FALSE
res$options[["Parameter set fitted at [Na+] (mM)"]]         # 100
```

Pick the set whose fitted salt is closest to your experimental condition rather
than correcting a distant one. See `?tm_nn` for the full list and citations.

## 4. launch an R shiny application

using R function
  `TmCalculatorShiny::TmCalculator_shiny()`

## 5. citation

If you use the melting-temperature-optimized parameter sets, please also cite
the source studies:

- Weber G (2015) *Bioinformatics* **31**:871–877.
  <doi:10.1093/bioinformatics/btu751>
- Ferreira I, Jolley EA, Znosko BM, Weber G (2019) *Chemical Physics*
  **521**:69–76. <doi:10.1016/j.chemphys.2019.01.016>
- Basilio Barbosa V, de Oliveira Martins E, Weber G (2019) *Biophysical
  Chemistry* **251**:106189. <doi:10.1016/j.bpc.2019.106189>