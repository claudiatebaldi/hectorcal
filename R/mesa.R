## This thing we're calling "erf" is missing the sqrt(2) factor inside the pnorm
## call that would ordinarily be there if we were calculating the Gaussian error
## function. However, when we compute the mesa function, if we were using the
## _real_ erf the calls would be erf((b-a)/(sqrt(2)*sig)), so dropping the
## sqrt(2) from the definition of erf allows us to drop it to from the arguments
## of the mesa function too, saving us a lot of sqrt(2) factors in our code that
## will only end up canceling out anyway.  It does, however, mean that this
## "erf" is not technically the real erf.
erf <- function(x) {2*pnorm(x)-1}

#' Calculate the mesa function used in defining our likelihoods.
#'
#' The mesa function is designed to be flat over some range and then to fall off
#' more or less sharply as it approaches the edges of that range.  Outside the
#' range, it has Gaussian-like tails.
#'
#' The parameters \code{a} and \code{b} control the width of the flat range.
#' One notable feature of the function is that the probability density at x=a or
#' x=b is always half what it is in the flat region of the function.
#'
#' The \code{sig} parameter controls how rapidly the probability density falls
#' off as you approach the boundary.  Larger values of \code{sig} extend the flat
#' region of the function closer to the boundary, and causes the rolloff to be
#' sharper, and the tails outside the boundary weaker.  Smaller values cause the
#' rolloff to start farther from the boundary, make the rolloff more gradual,
#' and cause the tails to extend farther outside the boundary.
#'
#' The function is vectorized, so it is permissible to pass a vector of x
#' values.  The other parameters will be recycled as vectors in the usual way.
#'
#' @param x Values at which to evaluate the mesa function
#' @param a Parameter controlling the lower boundary of the function.
#' @param b Parameter controlling the upper boundary of the function.
#' @param sig Parameter controlling the steepness of the rolloff at the
#' boundaries.
#' @export
mesa <- function(x, a, b, sig) {(erf((b-x)/sig) - erf((a-x)/sig))/(2*(b-a))}
