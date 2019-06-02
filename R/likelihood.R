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
#' @param inifiles Named vector of names of Hector input files for the various
#' experiments to be run. You do not need to provide ini files for historical
#' experiments unless they are the only experiments present.  Extraneous files
#' will be ignored, so it is safe to use a canoncial list of filenames.
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
#' @param hflux_year Year at which to retrieve the ocean heat flux. Default is 2100.
#' @param hflux_expt_regex Regular expression identifying the experiment to pull the heat
#' flux from.  There must be only one experiment in the dataset tha matches.
#' Default is 'rcp85'.
#' @param hflux_smoothing Smoothing parameter for heat flux.  Generally larger than
#' the smoothing parameter for the other variables because heat flux is less
#' well sampled than those variables are.
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
                            hflux_year=2100,
                            hflux_expt_regex='rcp85',
                            hflux_smoothing=0.2,
                            verbose=FALSE)
{
    ## Prior hyperparameters
    default_prior_params <- c(
        ecsmu=3.0, ecssig=log(3.0),
        aeromu=1.0, aerosig=1.4,
        volmu=1.0, volsig=1.4,
        kappamu=2.3, kappasig=2.0,
        betamu=0.3, betasig=0.7,
        q10mu=2.0, q10sig=2.0,
        c0mu=285, c0sig=14.0,
        sigtscale=1.0,
        sigco2scale=10.0,
        sigscale=1.0)

    ## Check to see if user has overridden any of the parameters
    if(!is.null(prior_params)) {
        for(param in names(default_prior_params)) {
            if(! param %in% names(prior_params)) {
                ## Set any unmentioned parameter to its default value
                prior_params[param] <- default_prior_params[param]
            }
        }
    }
    else {
        prior_params <- default_prior_params
    }


    logprior <- make_logprior(prior_params, use_lnorm_ecs)

    loglik <- make_loglikelihood(inifiles, verbose, cal_mean, comp_data,
                                 smoothing, hicol, lowcol, pcs,
                                 hflux_year, hflux_expt_regex, hflux_smoothing)

    ## Return the final posterior function
    make_logpost(logprior, loglik, verbose)
}

#' Internal functions for creating log-posterior functions
#'
#' These functions are called by \code{\link{build_mcmc_post}} to create the
#' log-posterior functions
#' @name likelihoodInternal
#' @keywords internal
NULL

#' @describeIn likelihoodInternal Make a log prior function
#'
#' @param pprior Parameters for prior distributions.  This \emph{must} be a
#' complete set (unlike the input to \code{build_mcmc_post})
#' @param use_lnorm_ecs Flag indicating whether to use the lognormal prior for
#' climate sensitivity
make_logprior <- function(pprior, use_lnorm_ecs)
{
    ## Check to see if user has requested the lognormal prior
    if(use_lnorm_ecs) {
        s_prior <- function(s) {stats::dlnorm(s, log(pprior['ecsmu']), pprior['ecssig'], log=TRUE)}
    }
    else{
        s_prior <- function(s) {stats::dnorm(s, pprior['ecsmu'], pprior['ecssig'], log=TRUE)}
    }

    ## truncated normal functions for constrained params
    betalprior <- mktruncnorm(0, Inf, pprior['betamu'], pprior['betasig'])
    q10prior   <- mktruncnorm(0, Inf, pprior['q10mu'], pprior['q10sig'])

    function(p) {
        sigvals <- p[c('sigt','sigco2', 'sig')] #use these later

        sum(
            ## normal priors
            stats::dnorm(
                p[c(hector::AERO_SCALE(), hector::VOLCANIC_SCALE(),
                    hector::DIFFUSIVITY(), hector::PREINDUSTRIAL_CO2())],
                pprior[c('aeromu', 'volmu', 'kappamu', 'c0mu')],
                pprior[c('aerosig', 'volsig', 'kappasig', 'c0sig')],
                log=TRUE),
            ## custom priors
            s_prior(p[hector::ECS()]),
            betalprior(p[hector::BETA()]),
            q10prior(p[hector::Q10_RH()]),
            ## half-cauchy priors for the sigvals
            ifelse(sigvals<0, -Inf, stats::dcauchy(sigvals, 0,
                                                   pprior[c('sigtscale','sigco2scale', 'sigscale')], log=TRUE)),
            ## This will cause us to ignore any parameters not present
            na.rm = TRUE
            )
    }
}


#' @describeIn likelihoodInternal Make a log likelihood function
#'
#' @param inifiles Named vector of names of ini files for the various
#' experiments. You do not need to provide ini files for historical
#' experiments unless they are the only experiments present.
#' @param verbose Flag indicating whether diagnostic messages should be
#' @param cal_mean Flag indicating whether we are calibrating to the mean. If
#' \code{FALSE}, calibrate to the envelope.
#' @param comp_data Observed comparison data, in summary format (see, for
#' example, \code{\link{esm_comparison}}.
#' @param smoothing Smoothing factor.
#' @param hicol Name of the column with the high edge of the window.
#' @param lowcol Name of the column with the low edge of the window.
#' @param pcs Principal components structure.  If \code{NULL} outputs will be
#' compared directly.
#' @param hflux_year Year at which to evaluate ocean heat flux for comparison.
#' @param hflux_expt_regex Regular expression identifying the experiment to
#' @param hflux_smoothing Smoothing factor for heat flux. (Because it is likely to
#' be less well sampled than other vars)
#' get heat flux from.
make_loglikelihood <- function(inifiles, verbose, cal_mean, comp_data,
                               smoothing, hicol, lowcol, pcs,
                               hflux_year, hflux_expt_regex, hflux_smoothing)
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

    ## Figure out which experiments we will be running
    if(is.null(pcs)) {
        expts <- unique(comp_data$experiment)
    }
    else {
        expts <- pcs$meta_data$experiment
    }

    ## Heat flux is treated specially, since we only use one value.  If there is
    ## a heat flux value in the comparison data, then split it off to be treated
    ## specially
    if(any(comp_data$variable == 'heatflux')) {
        hflux_data <- dplyr::filter(comp_data, variable=='heatflux')
        comp_data <- dplyr::filter(comp_data, variable != 'heatflux')
        expt85 <- grep(hflux_expt_regex, expts, value=TRUE)
        assert_that(length(expt85) == 1)
    }
    else {
        hflux_data <- NULL
        expt85 <- '--none--'
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

    ## Set up the hector core
    hcores <- setup_hector_cores(inifiles[names], names)

    ## If we are operating on output data, order the comparison data by
    ## experiment, variable, and year.
    if(is.null(pcs)) {
        comp_data <- dplyr::arrange(comp_data, experiment, variable, year)
    }
    else if(is.character(comp_data$variable)) {
        ## if the principal components are stored as "PCNN" convert to an
        ## integer index
        comp_data$variable <- as.integer(substring(comp_data$variable, 3))
    }


    ## All of the following are valid hector parameters.  We might have some
    ## other parameters, (e.g. sigma values), but those shouldn't be passed to
    ## parameterize_cores, or it will throw an error.
    hector_parms <- c(hector::ECS(), hector::DIFFUSIVITY(),
                       hector::AERO_SCALE(), hector::VOLCANIC_SCALE(),
                       hector::BETA(), hector::Q10_RH(), hector::PREINDUSTRIAL_CO2())


    ## Construct the likelihood function
    loglik <- function(p) {
        ## Protect all of these hector calls.  Any error will result in a
        ## -Inf result
        parm <- p[names(p) %in% hector_parms]
        hdata <- tryCatch({
            foreach(hcore=hcores, .combine = dplyr::bind_rows) %dopar% {
                ## Set the model parameters
                ## Run the model and retrieve variables.
                parameterize_core(parm, hcore)
                hector::run(hcore, futyrmax)
                expt <- getname(hcore)
                if(expt == expt85) {
                    ## This is the experiment designated to provide heat flux, so
                    ## retrieve it too.
                    hf <- hector::fetchvars(hcore, hflux_year, 'heatflux')
                }
                else {
                    hf <- NULL
                }
                dplyr::bind_rows(hf,
                                 hector::fetchvars(hcore, futyears, compvars))
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
            ## Ugh. When we run with dopar above, we can't be certain that the
            ## cores we have access to here have actually had their parameters
            ## set and been run.  This is really annoying, and we should
            ## probably just create a separate core for the historical runs, but
            ## that has hassles of its own.  For now, explicitly reset the
            ## parameters and rerun the core.
            parameterize_core(parm, hcores[[1]])
            hector::run(hcores[[1]], max(histyears))
            ## This loop should only execute once, but we try to cover our bases
            for(expt in histexpts) {
                hdata <- dplyr::bind_rows(hdata,
                                          hector::fetchvars(hcores[[1]], histyears, compvars, scenario = expt))
            }
        }

        ## The "scenario" column holds the experiment name.  Rename it.  Also
        ## correct the names of the variables.
        hdata <- dplyr::rename(hdata, experiment=scenario) %>%
          dplyr::mutate(variable=hvar2esmvar(variable))

        ## separate off the heatflux data.
        if(!is.null(hflux_data)) {
            cdata <- dplyr::filter(hdata, variable != 'heatflux')
            hdata <- dplyr::filter(hdata, variable == 'heatflux')
        }
        else {
            cdata <- hdata
            hdata <- NULL
        }

        ## Prepare the data for comparison to the hector data
        if(is.null(pcs)) {
            ## We're comparing the output directly, so in this case, all we need to do is arrange the data
            ## in the canonical order
            cdata <- dplyr::arrange(cdata, experiment, variable, year)
            ## This could theoretically fail, if the experiments have different years.  That's not
            ## part of our design, so we'll just throw an error if that happens.
            assert_that(nrow(cdata) == nrow(comp_data))
        }
        else {
            ## We need to project the data onto the principal components.  The projection coefficients
            ## will be used to do the comparison.
            pcproj <- project_climate(cdata, pcs, row_vector = FALSE)
            npc <- comp_data$variable
            cdata <- data.frame(variable=npc, value=pcproj[npc], stringsAsFactors=FALSE)
        }

        if(cal_mean) {
            ## Add the sigma factors for the normal distribution
            if(is.null(pcs)) {
                sigt <- if('sigt' %in% names(p)) p[['sigt']] else NA_real_
                sigco2 <- if('sigco2' %in% names(p)) p[['sigco2']] else NA_real_
                comp_data$sig <- dplyr::if_else(comp_data$variable == 'tas', sigt, sigco2)
                assert_that(!any(is.na(comp_data$sig)))   # If a sig parm is missing, that var must not be in the dataset
            }
            else {
                comp_data$sig <- p[['sig']]  # Only one type of sig for pca cal
            }
            ll <- sum(stats::dnorm(cdata$value, mean=comp_data$cmean,
                                   sd=comp_data$sig, log=TRUE))
            if(!is.null(hflux_data)) {
                assert_that('sighf' %in% names(p))
                sighf <- p[['sighf']]
                ll <- ll + stats::dnorm(hdata$value, mean=hflux_data$cmean,
                                        sd=sighf, log=TRUE)
            }
        }
        else {
            ## Compute the mesa function.  This actually works the same for both output and pc
            sig <- smoothing * (comp_data[[hicol]] - comp_data[[lowcol]])
            ll <- sum(log(mesa(cdata$value, comp_data[[lowcol]], comp_data[[hicol]], sig)))
            if(!is.null(hflux_data)) {
                hfsig <- hflux_smoothing*(hflux_data[[hicol]] - hflux_data[[lowcol]])
                ll <- ll + log(mesa(hdata$value, hflux_data[[lowcol]],
                                    hflux_data[[hicol]], hfsig))
            }
        }
        ## return the log of the likelihood
        ll
    }
}


#' @describeIn likelihoodInternal Combine prior and likelihood into a single function.
#'
#' @param lprior Function that computes the log-prior
#' @param llik Function that computes the log-posterior
#' @param verbose Flag indicating whether diagnostic output should be produced
make_logpost <- function(lprior, llik, verbose)
{
    function(p) {
        lp <- lprior(p)
        lpost <-
            if(is.finite(lp)) {
                lp + llik(p)
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
