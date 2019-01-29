
#### Error handler for hector errors
errhandler <- function(e)
{
    message(e)
    NULL
}

#' Construct a log posterior probability function for Bayesian calibration
#'
#' The function constructed by this call will be suitable for passing to either
#' a markov chain sampler or an optimizer.  (Though the optimizer will need to
#' be configured to find a maximum, rather than a minimum.)
#'
#' Note that the \code{use_c_cycle} argument only controls whether carbon cycle
#' parameters are sampled and whether the co2 variable is considered in the
#' likelihood function.  Whether Hector's carbon cycle runs is controled by the
#' input file \code{inifile}.  It is the user's responsibility to supply a
#' compatible input file, as well as comparison data that was run with the
#' desired carbon cycle settings.
#'
#' @section Notes:
#'
#' For some reason we have made the mesa function's \code{sig} parameter
#' settable for co2, but fixed at 0.4 for temperature.
#'
#' @param comp_data Table of esm summary statistics (mean, min, max) by year.
#' @param inifile Input file to initialize hector core.
#' @param years Years to use in the likelihood calculation.  Other years in the
#' comparison data and the Hector output will be ignored.
#' @param smooth_co2 Sig parameter to use in \code{\link{mesa}} for the CO2
#' variable.
#' @param cal_mean If true calibrate to mean; otherwise calibrate to range
#' @param use_c_cycle If true, include carbon cycle parameters; if not, don't.
#' @export
build_mcmc_post <- function(comp_data, inifile, years=seq(2010, 2100, 20),
                            smooth_co2 = 15,
                            cal_mean = TRUE, use_c_cycle=TRUE)
{
    ## indices in the parameter vector for the various parameters.  We have several
    ## combinations of run parameters, so we have to sort them out here.
    iecs <- 1
    iaero <- 2
    ikappa <- 3
    if(use_c_cycle) {
        ibeta <- 4
        iq10 <- 5
        ic0 <- 6
        ilast <- ic0
        compvars <- c(tas=hector::GLOBAL_TEMP(), co2=hector::ATMOSPHERIC_CO2())
    }
    else {
        ibeta <- iq10 <- ic0 <- NA
        ilast <- ikappa
        compvars <- c(tas=hector::GLOBAL_TEMP())
    }

    if(cal_mean) {
        ## We have a sigmat and possibly a sigmaco2 parameter
        isigt <- ilast+1
        if(use_c_cycle)
            isigco2 <- isigt+1
        else
            isigco2 <- NA
    }
    else {
        isigt <- isigco2 <- NA
    }

    #### Filter the comparison data to just what we are supposed to be using
    comp_data <- comp_data[comp_data$variable %in% names(compvars) & comp_data$year %in% years,]
    ## order by year and extract the values
    perm <- order(comp_data$year)
    comp_data <- comp_data[perm,]
    esmtemps <- comp_data[comp_data$variable=='tas', ]
    if(use_c_cycle) {
        esmco2 <- comp_data[comp_data$variable=='co2', ]
    }
    else {
        esmco2 <- NA
    }

    #### Set up the hector core
    hcore <- hector::newcore(inifile, suppresslogging = TRUE)

    ### Some hyperparameters
    ## Prior hyperparameters
    ecsmu   <- 3.0; ecssig   <- 3.0
    aeromu  <- 1.0; aerosig  <- 1.4
    kappamu <- 2.3; kappasig <- 2.0
    betamu  <- 0.3; betasig  <- 0.7
    q10mu   <- 2.0; q10sig   <- 2.0
    c0mu    <- 285; c0sig    <- 14.0
    sigtscale <- 1.0
    sigco2scale <- 10.0
    ## truncated normal functions for constrained params
    betalprior <- mktruncnorm(0, Inf, betamu, betasig)


    #### construct a function to return the log prior
    logprior <- function(p)
    {
        ## Get the normally distributed priors for the parameters that are
        ## always present.
        lp <- sum(dnorm(p[c(iaero, ikappa)],
                        c(aeromu, kappamu),
                        c(aerosig, kappasig),
                        log=TRUE),
                  dlnorm(p[iecs],
                         log(ecsmu),
                         log(ecssig),
                         log = TRUE))
        if(use_c_cycle) {
            lp <- lp + sum(dnorm(p[c(iq10, ic0)],
                                 c(q10mu, c0mu),
                                 c(q10sig, c0sig),
                                 log=TRUE)) +
                           betalprior(p[ibeta])
        }
        if(cal_mean) {
            if(use_c_cycle) {
                sigvals <- p[c(isigt, isigco2)]
                sclvals <- c(sigtscale, sigco2scale)
            }
            else {
                sigvals <- p[isigt]
                sclvals <- sigtscale
            }
            lp <- lp + sum(ifelse(sigvals<0, -Inf, dcauchy(sigvals, 0, sclvals, log=TRUE)))
        }
        ## return lp value
        lp
    }

    #### Construct the likelihood function
    loglik <- function(p)
    {
        ## Protect all of these hector calls.  Any error will result in a
        ## -Inf result
        hdata <- tryCatch({
            ## Set the model parameters
            hector::setvar(hcore, NA, hector::ECS(), p[iecs], 'degC')
            hector::setvar(hcore, NA, hector::AERO_SCALE(), p[iaero], NA)
            hector::setvar(hcore, NA, hector::DIFFUSIVITY(), p[ikappa], 'cm2/s')
            if(use_c_cycle) {
                hector::setvar(hcore, NA, hector::BETA(), p[ibeta], NA)
                hector::setvar(hcore, NA, hector::Q10_RH(), p[iq10], NA)
                hector::setvar(hcore, NA, hector::PREINDUSTRIAL_CO2(), p[ic0], 'ppmv CO2')
            }
            hector::reset(hcore)
            hector::run(hcore, max(years))
            hector::fetchvars(hcore, years, compvars)
        },
        error = errhandler)
        ## Check to see if there was an error. If so, return -Inf.
        if(is.null(hdata)) {
            return(-Inf)
        }

        htemps <- hdata[hdata$variable==hector::GLOBAL_TEMP(), 'value']
        if(cal_mean) {
            ll <- sum(dnorm(htemps, mean=esmtemps$cmean, sd=p[isigt], log=TRUE))
        }
        else {
            ll <- sum(log(mesa(htemps, esmtemps$mina, esmtemps$maxb, 0.4)))
        }
        if(use_c_cycle) {
            hco2 <- hdata[hdata$variable==hector::ATMOSPHERIC_CO2(), 'value']
            if(cal_mean) {
                ll <- ll + sum(dnorm(hco2, esmco2$cmean, sd=p[isigco2], log=TRUE))
            }
            else {
                ll <- ll + sum(log(mesa(hco2, esmco2$mina, esmco2$maxb, smooth_co2)))
            }
        }

        ## return the log of the likelihood
        ll
    }

    #### Finally, return a function that evaluates the prior and the likelihood
    function(p)
    {
        lp <- logprior(p)
        lpost <-
            if(is.finite(lp)) {
                lp + loglik(p)
            }
            else {
                -Inf
            }
        message('p: ', paste(p, collapse=', '), '      lpost: ', lpost)
        lpost
    }
}
