# TmCalculatorShiny

An interactive Shiny front end for the
[TmCalculator](https://github.com/JunhuiLi1017/TmCalculator) package.

## Install & run

```r
# install.packages("devtools")
devtools::install_local("TmCalculatorShiny")   # or load_all() from source
TmCalculatorShiny::runTmCalculatorShiny()
```

## Panels

- **Sequence Tm** — calculate melting temperature from one of four inputs:
  pasted sequences, genomic coordinate strings
  (`chr:start-end:strand:pkg_name:region_id`), an uploaded FASTA, or
  genome-wide tiling of an installed `BSgenome.*` package (chromosome, window,
  slide, start/end, strand, N-handling). Supports `tm_nn`, `tm_gc` (incl.
  custom **A/B/C/D** coefficients via `userset`), and `tm_wallace`. Every
  calculation is saved to a shared store.
- **Visualize** — linear, heatmap, and karyotype views of any stored result,
  plus multi-omics integration: upload feature ranges, run `integrate_granges`,
  and render multi-track linear / circular / karyotype genome plots.
- **Compare groups** — attach grouping metadata to a stored result (quantile
  bins, threshold split, an existing column, or pasted labels) and run
  `compare_groups`.
- **About** — overview and methods.
