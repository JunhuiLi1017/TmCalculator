#' Prints melting temperature from a \code{TmCalculator} object
#'
#' \code{print.TmCalculator} prints to console the melting temperature value from an object of
#' class \code{TmCalculator}.
#'
#' @param x An object of class \code{TmCalculator}.
#' @param ... Unused
#'
#' @return The melting temperature value.
#'
#' @export print.TmCalculator
#' @export
print.TmCalculator <- function(x,...){
  nohid <- attr(x, "nonhidden")
  print(x[[nohid]])
}

