# One warning per session rather than one per call: a script that profiles a
# hundred regions in a loop should be told once, not a hundred times.
.TmCalculatorEnv <- new.env(parent = emptyenv())

#' Genome-scale Tm profiling (deprecated)
#'
#' \code{tm_profile()} has been merged into \code{\link{tm_calculate}}, which
#' now accepts the same four sources and the same \code{regions} argument.
#' The two functions had come to take the same inputs and to differ only in
#' whether they spread the work over processes, which is a property of the
#' call rather than a reason for a second function.
#'
#' The translation is mechanical: the first argument keeps its meaning and
#' every other argument keeps its name.
#'
#' \preformatted{
#'   tm_profile(hg38, regions = 1:2, window = 200)
#'   tm_calculate(hg38, regions = 1:2, window = 200)
#' }
#'
#' One difference is worth knowing. \code{tm_profile()} returned a bare
#' \code{GRanges}; \code{tm_calculate()} returns a \code{TmCalculator}
#' object, so the profile is \code{$gr}. This wrapper returns the
#' \code{GRanges}, as it always did, so existing code keeps working until you
#' move it over.
#'
#' @param seq_source,... Passed to \code{\link{tm_calculate}}.
#' @return A \code{GRanges}, as before.
#' @seealso \code{\link{tm_calculate}}.
#' @rdname tm_profile
#' @export
tm_profile <- function(seq_source, ...) {
  .tm_profile_warned <- get0(".tm_profile_warned", envir = .TmCalculatorEnv,
                             ifnotfound = FALSE)
  if (!isTRUE(.tm_profile_warned)) {
    warning("tm_profile() is deprecated and will be removed in a future ",
            "release. It has been merged into tm_calculate(), which takes ",
            "the same arguments;\n  note that tm_calculate() returns a ",
            "TmCalculator object, so the profile is $gr.", call. = FALSE)
    assign(".tm_profile_warned", TRUE, envir = .TmCalculatorEnv)
  }
  tm_calculate(input_seq = seq_source, ...)$gr
}
