#### Functions for producing log-priors that we will be using in the analysis.

#' Create a truncated normal log density function
#'
#' Given the mean, standard deviation, and cutoffs for a truncated normal PDF,
#' produce a function that calculates the normalized log density.
#'
#' @param a Lower bound for the truncated normal distribution
#' @param b Upper bound for the truncated normal distribution
#' @param mu Mean for the untruncated normal distribution the distribution is
#' derived from.
#' @param sig Standard deviation for the untruncated normal distribution the
#' distribution is derived from.
#' @return A function that computes the log density, suitably normalized.
#' @export
mktruncnorm <- function(a, b, mu, sig)
{
    ftmp <- function(x) {dnorm(x, mu, sig)}
    lnfac <- log(integrate(ftmp, a, b)$value)     # normalization factor
    ## create the function to return
    function(x) {
        ifelse(x>=a & x<=b,
               dnorm(x,mu,sig, log=TRUE) - lnfac,
               -Inf)
    }
}

