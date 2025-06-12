# Tm Calculator

A web-based tool for calculating melting temperatures (Tm) of nucleic acid sequences with advanced visualization capabilities.

## Features

### Input Options
- Multiple input methods:
  - Direct sequence input
  - FASTA file upload
  - Genomic coordinates
- Support for both single and multiple sequences
- Support for ambiguous bases
- Real-time progress tracking during calculations

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

### Output Features
- Real-time calculation progress
- Interactive visualization
- Downloadable results and plots
- Multiple export formats

## How to Use

1. Select your preferred calculation method
2. Choose your input type (direct, FASTA, or genomic coordinates)
3. Enter your sequence(s) or upload file(s)
4. Adjust calculation parameters as needed
5. Click "Calculate Tm" to start the process
6. Monitor the calculation progress
7. View and interact with the results
8. Choose your preferred visualization type
9. Download results and plots for further analysis

## Citation

If you use Tm Calculator in your research, please cite:

```
Li, J., & [Your Institution]. (2024). TmCalculator: A web-based tool for calculating and visualizing melting temperatures of nucleic acid sequences. [Journal Name], [Volume], [Pages]. https://github.com/JunhuiLi1017/TmCalculator
```

## Contact

For questions, suggestions, or bug reports, please visit our [GitHub repository](https://github.com/JunhuiLi1017/TmCalculator/issues) or open an issue.

## About the Calculations

The Tm Calculator uses the nearest-neighbor thermodynamic parameters to calculate melting temperatures. The calculations take into account:

- Base composition
- Salt effects
- Chemical modifications
- Ambiguous bases