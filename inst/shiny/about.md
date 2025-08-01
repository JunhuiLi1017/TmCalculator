# TmCalculator

A Calculator for Melting Temperature of Nucleic Acid Sequences.

## Features

### Input Options
- Multiple input methods:
  - Direct sequence input
  - FASTA file upload
  - Genomic coordinates
- Support for both single and multiple sequences
- Support for ambiguous bases

### Calculation Methods
- Nearest Neighbor (tm_nn)
- GC Content (tm_gc)
- Wallace Rule (tm_wallace)

### Adjustable Parameters
- Salt conditions:
  - Na+ concentration
  - K+ concentration
  - Tris buffer
  - Mg2+ concentration
  - dNTPs concentration
- Chemical modifications:
  - DMSO percentage
  - Formamide (percent or molar)
- Multiple salt correction methods available

### Visualization Options
- Karyotype Plot:
  - Customizable chromosome colors and shapes
  - Adjustable point sizes and text scaling
  - Flexible layout options
- Heatmap Plot:
  - Karyogram and faceted views
  - Customizable color palettes
  - Zoom functionality for detailed analysis
- Genome Tracks:
  - Interactive ideogram display
  - Customizable color schemes
  - Region-specific zooming

## How to Use

1. Select your preferred calculation method
2. Choose your input type (direct, FASTA, or genomic coordinates)
3. Enter your sequence(s) or upload file(s)
4. Adjust calculation parameters as needed
5. Click "Calculate Tm" to start the process
6. View and interact with the results
7. Choose your preferred visualization type
8. Download results and plots for further analysis

## Citation

If you use Tm Calculator in your research, please cite:

```
Junhui Li (2025). TmCalculator: A Calculator for Melting Temperature of Nucleic Acid Sequences. R package version 1.1.0, https://CRAN.R-project.org/package=TmCalculator
```

## Contact

For questions, suggestions, or bug reports, please visit our [GitHub repository](https://github.com/JunhuiLi1017/TmCalculator/issues) or open an issue.
