# Tm Calculator

A web-based tool for calculating melting temperatures (Tm) of nucleic acid sequences.

## Features

- Calculate Tm using both primer and template sequences
- Support for ambiguous bases
- Adjustable salt conditions:
  - Na+ concentration
  - K+ concentration
  - Tris buffer
  - Mg2+ concentration
  - dNTPs concentration
- Chemical modifications:
  - DMSO percentage
  - Formamide (percent or molar)
- Download results for further analysis

## How to Use

1. Enter your primer sequence (5' to 3')
2. Enter your template sequence (5' to 3')
3. Adjust salt conditions if needed
4. Add chemical modifications if required
5. Click "Calculate Tm" to get results
6. Download results if needed

## About the Calculations

The Tm Calculator uses the nearest-neighbor thermodynamic parameters to calculate melting temperatures. The calculations take into account:

- Base composition
- Salt effects
- Chemical modifications
- Ambiguous bases

## References

For more information about the thermodynamic parameters and calculation methods, please refer to:

- SantaLucia J. (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95:1460-1465
- Owczarzy R. et al. (2008) Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations. Biochemistry 47:5336-5353 