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



#' @describeIn mktruncnorm  Sample from a truncated normal log density function
#' @param n Number of samples to generate
#' @export
rtruncnorm <- function(n, a, b, mu, sig)
{
    samps <- rnorm(n, mu, sig)
    bad <- samps<a | samps>b
    nbad <- sum(bad)
    while(nbad > 0) {
        samps[bad] <- rnorm(nbad, mu, sig)
        bad <- samps<a | samps>b
        nbad <- sum(bad)
    }

    samps
}

#' @describeIn mktruncnorm CDF for truncated normal distribution
#' @param x Vector of x values
#' @param lower.tail If true, probabilities are P(X<x); if false, probabilities
#' are P(X>x).
#' @param log.p If true, return log probabilities
#' @export
ptruncnorm <- function(x, a, b, mu, sig, lower.tail=TRUE, log.p=FALSE)
{
    if(length(a) == 1 && length(b) ==1 && length(mu)==1 && length(sig)==1) {
        ## Normalization factor will be the same for all x values, so save
        ## ourselves some unnecessary computation by calculating it once.
        ftmp <- function(x) {dnorm(x, mu, sig)}
        normfac <- integrate(ftmp, a, b)$value
    }
    else {
        normfac <- NULL
    }
    mapply(ptnsingle, x, a, b, mu, sig, lower.tail, log.p, normfac)
}

### Helper function for ptruncnorm.  Each of these arguments should be a single
### value (as opposed to a possible vector for the ones in prtruncnorm)
ptnsingle <- function(x, a, b, mu, sig, lower.tail, log.p, normfac=NULL)
{
    val <-
        if(x <= a) {
            0
        }
        else if (x >= b) {
            1
        }
        else {
            ftmp <- function(x) {dnorm(x, mu, sig)}
            if(is.null(normfac)) {
                normfac <- integrate(ftmp, a, b)$value
            }
            integrate(ftmp, a, x)$value / normfac
        }

    if(!lower.tail)
        val <- 1.0 - val

    if(log.p)
        log(val)
    else
        val
}
