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
#' @param inifiles Vector of input files to initialize hector cores.
#' @param pcs Principal components structure.  If \code{NULL}, work directly
#' with output values.
#' @param smoothing Smoothing parameter for the mesa function, expressed as a
#' fraction of the window size.
#' @param cal_mean If true calibrate to mean; otherwise calibrate to range
#' @param lowcol Column in comparison data to use for the low edge of the mesa
#' function.  Ignored if cal_mean is \code{TRUE}
#' @param hicol Column in comparison data to use for the high edge of the mesa
#' function.  Ignored if cal_mean is \code{TRUE}
#' @param prior_params Named list of alternative values for the numerical
#' parameters in the prior distributions.  Any parameters not mentioned in the
#' list will be set to their default values.
#' @param use_lnorm_ecs If true, use a log-normal prior for climate
#' sensitivity.  If false, use a normal prior.  The false setting is intended
#' primarily for testing the influence of priors on the final result.
#' @param verbose If \code{TRUE}, print diagnostic messages.  Otherwise,
#' diagnostic messages will not be generated.
#' @export
build_mcmc_post <- function(comp_data, inifiles,
                            pcs = NULL,
                            smoothing = 0.1,
                            cal_mean = TRUE,
                            lowcol = 'mina', hicol = 'maxb',
                            prior_params = NULL,
                            use_lnorm_ecs=TRUE,
                            verbose=FALSE)
{
    ## create bindings for NSE vars
    experiment <- variable <- year <- scenario <- NULL

    ## set up the error handler
    if(verbose) {
        errhandler <- function(e) {
            message(conditionMessage(e))
            NULL
        }
    }
    else {
        errhandler <- function(e) {
            NULL
        }
    }

    ## All of the following are valid hector parameters.  We might have some
    ## other parameters, (e.g. sigma values), but those shouldn't be passed to
    ## parameterize_cores, or it will throw an error.
    hector_parms <- c(hector::ECS(), hector::DIFFUSIVITY(),
                       hector::AERO_SCALE(), hector::VOLCANIC_SCALE(),
                       hector::BETA(), hector::Q10(), hector::PREINDUSTRIAL_CO2())


    ## Figure out which parameters we will need to pull for comparison
    if(is.null(pcs)) {
        ## Always compare temperature, compare co2 only if co2 variable is
        ## present in comparison data
        cv <- comp_data$variable
    }
    else {
        cv <- pcs$meta_data$variable
    }
    if('co2' %in% cv) {
        compvars <- c(hector::GLOBAL_TEMP(), hector::ATMOSPHERIC_CO2())
    }
    else {
        compvars <- hector::GLOBAL_TEMP()
    }

    ## Figure out which years we will need for comparison
    if(is.null(pcs)) {
        compyears <- unique(comp_data$year)
        histyears <- compyears[compyears <= 2005]
        futyears <- compyears[compyears >= 2006]
    }
    else {
        histyears <- pcs$meta_data$histyear
        futyears <- pcs$meta_data$year
    }
    histyrmax <- max(histyears)
    futyrmax <- max(futyears)


    ## Figure out which experiments we will be running
    if(is.null(pcs)) {
        expts <- unique(comp_data$experiment)
    }
    else {
        expts <- pcs$meta_data$experiment
    }
    ## We don't need to run the historical experiment separately, unless it's
    ## the only one, since we can get it from any of the other runs.
    histexpts <- grep('[Hh]istorical', expts, value=TRUE)
    if(length(expts) > 1) {
        ## get names of future experiments (any that aren't historical)
        names <- grep('[Hh]istorical', expts, value=TRUE, invert=TRUE)
    }
    else {
        ## With only a single experiment, make sure we use it, even if it's historical.
        names <- expts
    }

    ## If we are operating on output data, order the comparison data by
    ## experiment, variable, and year
    if(is.null(pcs)) {
        comp_data <- dplyr::arrange(comp_data, experiment, variable, year)
    }


    ## Set up the hector core
    hcores <- setup_hector_cores(inifiles, names)

    ## Prior hyperparameters
    ecsmu   <- 3.0; ecssig   <- 3.0
    aeromu  <- 1.0; aerosig  <- 1.4
    kappamu <- 2.3; kappasig <- 2.0
    betamu  <- 0.3; betasig  <- 0.7
    q10mu   <- 2.0; q10sig   <- 2.0
    c0mu    <- 285; c0sig    <- 14.0
    sigtscale <- 1.0
    sigco2scale <- 10.0

    ## Check to see if user has overridden any of the parameters
    if(!is.null(prior_params)) {
        for(param in names(prior_params)) {
            assign(param, prior_params[[param]])
        }
    }

    ## Check to see if user has requested the lognormal prior
    if(use_lnorm_ecs) {
        s_prior <- function(s) {stats::dlnorm(s, log(ecsmu), log(ecssig), log=TRUE)}
    }
    else{
        s_prior <- function(s) {stats::dnorm(s, ecsmu, ecssig, log=TRUE)}
    }

    ## truncated normal functions for constrained params
    betalprior <- mktruncnorm(0, Inf, betamu, betasig)
    q10prior   <- mktruncnorm(0, Inf, q10mu, q10sig)


#### construct a function to return the log prior
    logprior <- function(p)
    {
        ## Get the normally distributed priors for the parameters that are
        ## always present.
        lp <- sum(stats::dnorm(p[c(iaero, ikappa)],
                        c(aeromu, kappamu),
                        c(aerosig, kappasig),
                        log=TRUE),
                  s_prior(p[iecs]))
        if(use_c_cycle) {
            lp <- lp + sum(stats::dnorm(p[ic0],
                                 c0mu,
                                 c0sig,
                                 log=TRUE),
                           betalprior(p[ibeta]),
                           q10prior(p[iq10]))


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
            lp <- lp + sum(ifelse(sigvals<0, -Inf, stats::dcauchy(sigvals, 0, sclvals, log=TRUE)))
        }
        ## return lp value
        lp
    }

#### Construct the likelihood function
    loglik <- function(p) {
        ## Protect all of these hector calls.  Any error will result in a
        ## -Inf result
        hdata <- tryCatch({
            foreach(hcore=hcores, .combine = dplyr::bind_rows) %dopar% {
                ## Set the model parameters
                parm <- p[names(p) %in% hector_parms]
                parameterize_core(parm, hcore)
                ## Run the model and retrieve variables.
                hector::run(hcore, futyrmax)
                hector::fetchvars(hcore, futyears, compvars)
            }
        },
        error = errhandler)
        ## Check to see if there was an error. If so, return -Inf.
        if(is.null(hdata)) {
            return(-Inf)
        }

        ## If we have a historical experiment in the mix, get its data.  It doesn't matter which
        ## core we use for this, so pick the first one arbitrarily.
        if(length(histexpts) > 0) {
            ## This loop should only execute once, but we try to cover our bases
            for(expt in histexpts) {
                hdata <- dplyr::bind_rows(hdata,
                                          hector::fetchvars(hcores[[1]], histyears, compvars, scenario = expt))
            }
        }

        ## The "scenario" column holds the experiment name.  Rename it
        hdata <- dplyr::rename(hdata, experiment=scenario)

        ## Prepare the data for comparison to the hector data
        if(is.null(pcs)) {
            ## We're comparing the output directly, so in this case, all we need to do is arrange the data
            ## in the canonical order
            cdata <- dplyr::arrange(hdata, experiment, variable, year)
            ## This could theoretically fail, if the experiments have different years.  That's not
            ## part of our design, so we'll just throw an error if that happens.
            assert_that(nrow(cdata) == nrow(comp_data))
        }
        else {
            ## We need to project the data onto the principal components.  The projection coefficients
            ## will be used to do the comparison.
            pcproj <- project_climate(hdata, pcs, row_vector = FALSE)
            npc <- comp_data$PC
            cdata <- data.frame(PC=npc, value=pcproj[npc])
        }

        if(cal_mean) {
            ## Add the sigma factors for the normal distribution
            if(is.null(pcs)) {
                sigt <- if('sigt' %in% names(p)) p[['sigt']] else as.numeric(NA)
                sigco2 <- if('sigco2' %in% names(p)) p[['sigco2']] else as.numeric(NA)
                comp_data$sig <- dplyr::if_else(comp_data$variable == 'tas', sigt, sigco2)
                assert_that(!any(is.na(com_data$sig)))   # If a sig parm is missing, that var must not be in the dataset
            }
            else {
                comp_data$sig <- p[['sig']]  # Only one type of sig for pca cal
            }
            ll <- sum(stats::dnorm(cdata$value, mean=compdata$cmean, sd=compdata$sig, log=TRUE))
        }
        else {
            ## Compute the mesa function.  This actually works the same for both output and pc
            sig <- smoothing * (comp_data[[hicol]] - comp_data[[lowcol]])
            ll <- sum(log(mesa(cdata$value, comp_data[[lowcol]], comp_data[[hicol]], sig)))
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
        if(verbose) {
            message('p: ', paste(p, collapse=', '), '      lpost: ', lpost)
        }
        lpost
    }
}
