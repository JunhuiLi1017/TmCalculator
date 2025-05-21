# TmCalculatorShiny

A Shiny web application that provides a user-friendly interface for the TmCalculator package.

Version: 1.0.0 (2024-03-19)

## Features

- Three different methods for calculating melting temperature (Tm):
  - GC-based Method
  - Nearest Neighbor Method
  - Wallace Method
- Support for various parameters:
  - Salt concentrations (Na+, K+, Tris, Mg2+, dNTPs)
  - Chemical modifiers (DMSO, formamide)
  - Ambiguous bases
  - Self-complementary sequences
  - Different hybridization types

## Installation

1. Make sure you have R installed (version 4.0.0 or higher recommended)
2. Install required packages:

```r
install.packages(c("shiny", "shinythemes", "bslib", "DT"))
```

3. Install the TmCalculator package:

```r
install.packages("TmCalculator")
```

4. Clone this repository or download the source code

## Usage

1. Open R or RStudio
2. Set your working directory to the TmCalculatorShiny folder
3. Run the app:

```r
shiny::runApp()
```

## Interface

The application provides three main tabs:

1. **GC-based Method**: Uses empirical formulas based on GC content
2. **Nearest Neighbor Method**: Uses thermodynamic parameters for nearest neighbor interactions
3. **Wallace Method**: Uses a simple formula based on GC content and sequence length

Each method has its own set of parameters that can be adjusted through the interface.

## Help

The Help tab provides information about:
- Available methods
- Input parameters
- References

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

Junhui Li (junhui.li@umassmed.edu) 