

#' Plot spaghetti plots for hector samples
#'
#' Spaghetti plots show the traces of hector outputs for a set of samples from a
#' Bayesian Monte Carlo calculation.
#'
#' @param mcrslt List of metrosamp structures from a Bayesian Monte Carlo run
#' @param nplot Number of traces to plot.
#' @param hcores List of hector cores initialized with the desired emissions.
#' @param pnames Hector names of the parameters in the sample matrix
#' @param vars Variables to plot
#' @param times Vector of times to include in the plot.  Default is 1850-2100.
#' @param alpha Alpha channel (transparency) for the sample traces.
#' @export
spaghetti_plot <- function(mcrslt, nplot, hcores, pnames,
                           vars=c(hector::GLOBAL_TEMP(), hector::ATMOSPHERIC_CO2()),
                           times=1850:2100, alpha=0.3)
{
    allsamps <- do.call(rbind, lapply(mcrslt, function(x) {x$samples}))

    sampout <- run_hector_samples(allsamps, nplot, hcores, times, pnames, vars)

    ## create a mean output dataset
    meanout <-
        dplyr::group_by(sampout, variable, year) %>%
          summarise(value = mean(value))

    ggplot2::ggplot(data=sampout, ggplot2::aes(x=year, y=value)) +
      ggplot2::geom_line(ggplot2::aes(group=scenario), color='blue', alpha=alpha) +
      ggplot2::geom_line(data=meanout, color='red', size=1.25) +
      ggplot2::facet_wrap(facets=~variable, scales='free_y')

}



## Helper function for running the hector samples for various plots
run_hector_samples <- function(allsamps, nsamp, hcores, times, pnames, vars)
{
    ntot <- nrow(allsamps)
    samp_indices <- sample.int(ntot, nsamp)

    if(!is.list(hcores)) {
        ## If a single core was passed, convert it into a list of length one
        hcores <- list(hcores)
    }

    ncore <- length(hcores)
    hout <-
            foreach(icore=seq(1,ncore), .combine=rbind) %dopar% {
                idx <- samp_indices[seq(icore, nsamp, ncore)]
                hcore <- hcores[[icore]]
                foreach(i=seq_along(idx), .combine=rbind) %do% {
                    irow <- idx[i]
                    runid <- i*ncore + icore
                    for(iparm in seq_along(pnames)) {
                        hector::setvar(hcore, NA, pnames[iparm], allsamps[irow,iparm],
                                       hector::getunits(pnames[iparm]))
                    }
                    hector::reset(hcore)
                    hector::run(hcore, max(times))
                    hector::fetchvars(hcore, times, vars,
                                      scenario=as.character(runid))
                }
            }

}
