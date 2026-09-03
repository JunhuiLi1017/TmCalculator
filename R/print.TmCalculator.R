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
#' @export

print.TmCalculator <- function(x,...){
  nohid <- attr(x, "nonhidden")
  print(x[[nohid]])
}

#' Lazy \code{$df} accessor for \code{TmCalculator} objects
#'
#' \code{result$df} is computed on demand from \code{result$gr} instead of
#' being stored eagerly, so genome-scale results are not converted to a
#' data.frame unless actually requested. All other elements behave as plain
#' list access.
#'
#' @param x A \code{TmCalculator} object.
#' @param name Element name.
#' @return The element, or a data.frame built from \code{x$gr} when
#'   \code{name == "df"} and no stored \code{df} exists.
#' @export
`$.TmCalculator` <- function(x, name) {
  val <- .subset2(x, name)
  if (is.null(val) && identical(name, "df")) {
    gr <- .subset2(x, "gr")
    if (!is.null(gr)) {
      val <- as.data.frame(gr)
    }
  }
  val
}

