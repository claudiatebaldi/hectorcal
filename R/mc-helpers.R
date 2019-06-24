## Functions for working with mc output

#' Read in results from Monte Carlo runs
#'
#' This function searches a directory for Monte Carlo outputs and loads them.
#' Individual runs from a single configuration are gathered together into lists
#' and named according to the runid for that configuration with serial number 0.
#'
#' @param dir Filesystem directory in which to look for results.
#' @param filestem Base name used to construct the names of the input files.
#' @return A list of structures.
#' \describe{
#' \item{runstats}{Run statistics such as acceptance rate and effective number of samples.
#' These apply to the merged set of runs for the configuration.}
#' \item{mcobjs}{List of lists of \code{metrosamp} objects.  Each list in the top level
#' corresponds to a configuration; the elements in the second-level list are the objects
#' for the individual runs in the configuration.}
#' \item{codaobjs}{\code{\link[coda]{mcmc.list}} object for each configuration.  These objects
#' represent the merged results of all of the runs in a configuration in a form that can be used
#' with the utilities in the coda package.}
#' }
#' @export
proc_mc_rslts <- function(dir='.', filestem='hectorcal') {
    configs <- expand.grid(pcs=c(TRUE,FALSE), hf=c(TRUE,FALSE), meanf=c(TRUE,FALSE))

    runstats <- list()
    mcobjs <- list()
    codaobjs <- list()

    for (configrow in seq(1,nrow(configs))) {
        conf <- configs[configrow,]
        base_runid <- encode_runid(0, conf$pcs, conf$hf, conf$meanf, FALSE)
        runname <- as.character(base_runid)
        mcobj <- load_mc_output(conf$pcs, conf$hf, conf$meanf, FALSE, dir = dir, filestem = filestem)
        codaobj <- metrosamp2coda(mcobj)

        avgaccept <- mean(sapply(mcobj, function(x){x$accept}))

        runstats[[runname]] <- list(accept=avgaccept, neff=coda::effectiveSize(codaobj),
                                    rhat=coda::gelman.diag(codaobj))
        mcobjs[[runname]] <- mcobj
        codaobjs[[runname]] <- codaobj
    }

    invisible(list(runstats=runstats, mcobjs=mcobjs, codaobjs=codaobjs))
}

#' Make diagnostic plots for Monte Carlo Runs
#'
#' Making these plots is very slow for large runs.  Mostly they are
#' useful for diagnosing pilot runs
#'
#' @param plotruns List of metrosamp objects for runs in a single configuration.
#' @param tofile Flag indicating whether the result should be written to a file
#' @param filestem String giving the stem of the filename to be written.
#' @export
#' @keywords internal
makeplots <- function(plotruns, tofile=TRUE, filestem='hectorcal-run-') {
    for (runname in names(plotruns)) {
        runid0 <- as.integer(runname)
        runid15 <- runid0+15
        filename <- paste0(filestem, runid0, '-', runid15, '.png')
        if(tofile) {
            if(!file.exists(filename)) {
                png(filename, width=1024, height=1024)
                rundefn <- decode_runid(runid0)
                titlestr <- sprintf('runid %d  PC=%s hf=%s mn=%s', runid0,
                                    rundefn$pcsflag, rundefn$hfflag, rundefn$meanflag)
                title(titlestr)
                pltobj <- metrosamp2coda(plotruns[[runname]])
                plot(pltobj)
                dev.off()
            }
        }
        else {
            ## Kind of ugly to repeat this code here.  Replace base plots with appropriate ggplots
            ## so we can create the plot object and then just print it to the appropriate device.
            rundefn <- decode_runid(runid0)
            titlestr <- sprintf('runid %d  PC=%s hf=%s mn=%s', runid0,
                                rundefn$pcsflag, rundefn$hfflag, rundefn$meanflag)
            title(titlestr)
            pltobj <- metrosamp2coda(plotruns[[runname]])
            plot(pltobj)
        }
    }
}


#' Plot pairs plots for the samples in a Monte Carlo configuration.
#'
#' The pairs plots can be used to get some insight into the joint distribution
#' of the parameters.
#'
#' @param mcrslt List of metrosamp objects for runs in a configuration
#' @return The correlation matrix for the samples
#' @export
pairplot <- function(mcrslt) {
    samps <- do.call(rbind, lapply(mcrslt, function(x) {x$samples}))
    nsamp <- nrow(samps)
    idx <- seq(1, nsamp)
    if(nsamp > 500) {
        fac <- ceiling(nsamp/500)
        pltsamps <- samps[idx%%fac == 0, ]
    }
    else {
        pltsamps <- samps
    }
    pairs(pltsamps)
    cor(samps)
}
