#' Obtaining the inverse of marginal empirical cumulative distribution (ECDF)
#'
#' @param x A vector of numbers which is the marginal empirical data.
#' @param p A vector of numbers which is the probability of the simulated data.
#' @param sort.flag A logical value to specify whether to sort the output data.
#' @return The inverse values of \code{p} based on ECDF of \code{x}.
#' @examples
#' ecdf.inv(0:10, c(0.25, 0.75))
#' ecdf.inv(0:10, c(0.25, 0.75), FALSE)
#' @export
#' @importFrom stats approx
#' @importFrom dplyr between
ecdf.inv <- function(x, p, sort.flag = TRUE) {
  if (any(is.na(x)) || any(is.na(p)))
    stop("x or p should not contains NAs.")
  if (!all(between(p, 0, 1)))
    stop("p should be between 0 and 1.")
  p <- if (sort.flag) {
    sort(p)
  } else {
    p
  }
  return(approx(
    seq_along(x),
    sort(x),
    (length(x) - 1) * p + 1,
    "linear",
    yleft = min(x),
    yright = max(x)
  )$y)
}
