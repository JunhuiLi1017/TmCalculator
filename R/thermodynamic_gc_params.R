#' Thermodynamic parameters for GC-based Tm calculation methods
#' 
#' A data frame containing coefficients and parameters for different GC-based Tm calculation methods.
#' Each row represents a different method with its specific coefficients (A, B, C, D) and salt correction method.
#' 
#' @format A data frame with 8 rows and 5 columns:
#' \describe{
#'   \item{A}{Intercept coefficient}
#'   \item{B}{GC content coefficient}
#'   \item{C}{Length correction coefficient}
#'   \item{D}{Mismatch coefficient}
#'   \item{salt_correction}{Associated salt correction method}
#' }
#' 
#' @details
#' The methods included are:
#' - Chester1993: Tm = 69.3 + 0.41(Percentage_GC) - 650/N
#' - QuikChange: Tm = 81.5 + 0.41(Percentage_GC) - 675/N - Percentage_mismatch
#' - Schildkraut1965: Tm = 81.5 + 0.41(Percentage_GC) - 675/N + 16.6 x log[Na+]
#' - Wetmur1991_MELTING: Tm = 81.5 + 0.41(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#' - Wetmur1991_RNA: Tm = 78 + 0.7(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#' - Wetmur1991_RNA/DNA: Tm = 67 + 0.8(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
#' - Primer3Plus: Tm = 81.5 + 0.41(Percentage_GC) - 600/N + 16.6 x log[Na+]
#' - vonAhsen2001: Tm = 77.1 + 0.41(Percentage_GC) - 528/N + 11.7 x log[Na+]
#' 
"thermodynamic_gc_params"

thermodynamic_gc_params <- data.frame(
  A = c(69.3, 81.5, 81.5, 81.5, 78.0, 67.0, 81.5, 77.1),
  B = c(0.41, 0.41, 0.41, 0.41, 0.70, 0.80, 0.41, 0.41),
  C = c(650, 675, 675, 500, 500, 500, 600, 528),
  D = rep(1, 8),
  salt_correction = c(NA, NA, "Schildkraut2010",
                     rep("Wetmur1991", 3), "Schildkraut2010", "SantaLucia1998-1")
)
rownames(thermodynamic_gc_params) <- c(
  "Chester1993", "QuikChange", "Schildkraut1965",
  "Wetmur1991_MELTING", "Wetmur1991_RNA", "Wetmur1991_RNA/DNA",
  "Primer3Plus", "vonAhsen2001"
) 

usethis::use_data(thermodynamic_gc_params, overwrite = TRUE) 