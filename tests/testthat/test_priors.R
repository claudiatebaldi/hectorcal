context('Custom probability distribution functions')


test_that('Truncated normal CDF is correct', {
    allvals <- c(-1, 0, 1, 2, 3)
    probs <- ptruncnorm(allvals, 0, 2, 0, 1)
    expect_equal(probs[1:2], c(0,0))
    expect_equal(probs[4:5], c(1,1))
    expect_gt(probs[3], 0)
    expect_lt(probs[3], 1)

    upper.probs <- ptruncnorm(allvals, 0, 2, 0, 1, lower.tail=FALSE)
    expect_equal(upper.probs, 1-probs)

    ## If we truncate from 0 to infinity, then the upper tails for positive x
    ## should be exactly twice the upper tails for the standard normal.
    plusvals <- c(0, 1, 2, 3)
    probs <- ptruncnorm(plusvals, 0, Inf, 0, 1, lower.tail=FALSE)
    ## Check the log probs while we're at it
    lprobs <- ptruncnorm(plusvals, 0, Inf, 0, 1, lower.tail=FALSE, log.p=TRUE)
    expect_equal(lprobs, log(probs))

    chkprobs <- pnorm(plusvals, 0, 1, lower.tail=FALSE)
    expect_equal(probs, 2.0*chkprobs)

    ## Check that changing mu and sigma works as expected
    plusvals <- c(2, 3, 4, 5)
    probs <- ptruncnorm(plusvals, 2, Inf, 2, 2, lower.tail=FALSE)
    chkprobs <- pnorm(plusvals, 2, 2, lower.tail=FALSE)
    expect_equal(probs, 2.0*chkprobs)
})


test_that('Truncated normal PDF is correct', {
    ltruncnorm <- mktruncnorm(0, 5, 0, 2)
    dtruncnorm <- function(x) {exp(ltruncnorm(x))}
    xval <- 2
    pval <- integrate(dtruncnorm, -Inf, xval)$value

    expect_equal(pval, ptruncnorm(xval, 0, 5, 0, 2))
})


test_that('Values returned by rtruncnorm have the right distribution', {
    set.seed(867-5309)
    rvals <- rtruncnorm(1000, 0, 2, 0, 1)
    expect_true(all(rvals >= 0))
    expect_true(all(rvals <= 2))

    ptcdf <- function(x) {ptruncnorm(x, 0, 2, 0, 1)}
    kst <- ks.test(rvals, ptcdf)
    expect_gt(kst$p.value, 0.5)
})
