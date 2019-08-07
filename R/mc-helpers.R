## Functions for working with mc output

#' Read in results from Monte Carlo runs
#'
#' This function searches a directory for Monte Carlo outputs and loads them.
#' Individual runs from a single configuration are gathered together into lists
#' and named according to the runid for that configuration with serial number 0.
#'
#' @param dir Filesystem directory in which to look for results.
#' @param filestem Base name used to construct the names of the input files.
#' @param codasize Size of the \code{\link[coda]{mcmc.list}} objects returned.  The MCMC
#' results will be thinned to this size for easier plotting and analysis.  If
#' \code{codasize} is specified as \code{NA}, then the output coda objects will not be
#' thinned at all.
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
proc_mc_rslts <- function(dir='.', filestem='hectorcal', codasize=2500)
{
    configs <- expand.grid(pcs=c(TRUE,FALSE), hf=c(TRUE,FALSE), meanf=c(TRUE,FALSE))

    runstats <- list()
    mcobjs <- list()
    codaobjs <- list()

    for (configrow in seq(1,nrow(configs))) {
        conf <- configs[configrow,]
        base_runid <- encode_runid(0, conf$pcs, conf$hf, conf$meanf, FALSE)
        runname <- as.character(base_runid)
        mcobj <- load_mc_output(conf$pcs, conf$hf, conf$meanf, FALSE, dir = dir, filestem = filestem)
        if(length(mcobj) > 0) {
            stat <- proc_mc_helper(mcobj, codasize)
            runstats[[runname]] <- stat[[1]]
            mcobjs[[runname]] <- mcobj
            codaobjs[[runname]] <- stat[[2]]
        }
    }

    invisible(list(runstats=runstats, mcobjs=mcobjs, codaobjs=codaobjs))
}

## mcobj is a list of metrosamp objects loaded by load_mc_output
proc_mc_helper <- function(mcobj, size)
{
    codaobj_full <- metrosamp2coda(mcobj)
    avgaccept <- mean(sapply(mcobj, function(x){x$accept}))

    runstats <- list(accept=avgaccept, neff=coda::effectiveSize(codaobj_full),
                                rhat=coda::gelman.diag(codaobj_full))
    codaobj <- metrosamp2coda(mcobj, size=size)
    list(runstats, codaobj)
}

#' Concatenate a set of base runs and continuations into a single run structure.
#'
#' The inputs to this function are two sets of outputs from \code{\link{proc_mc_rslts}}.  The
#' function will extract the metrosamp objects from the two outputs, concatenate them using
#' the \code{\link[metrosamp]{concat}} function, and produce the \code{runstats} and \code{codaobj}
#' structures erturned by \code{proc_mc_rslts}.  Optionally, write the concatenated results to
#' files so that they can be read in directly in the future.
#'
#' @param bruns List of Monte Carlo results from the base runs.
#' @param cruns List of Monte Carlo results from the continuation runs.
#' @param output_dir Optional output directory to write concatenated objects into.
#' @param name_stem String giving the stem of the filenames for saving concatenated objects.
#' @param codasize Number of samples to keep in the \code{\link[coda]{mc.list}} objects in the output.
#' @return Structure as in \code{\link{proc_mc_rslts}}, for the concatenated runs
#' @export
concat_runs <- function(bruns, cruns, output_dir = NULL, name_stem = 'hectorcal_concat',
                        codasize=2500)
{
    assertthat::assert_that(assertthat::are_equal(names(bruns), names(cruns)))
    runstats <- list()
    mcobjs <- list()
    codaobjs <- list()

    for(run in names(bruns$mcobjs)) {
        if(run %in% names(cruns$mcobjs)) {
            ## concatenate runs that have continuations.
            base <- bruns$mcobjs[[run]]
            contin <- cruns$mcobjs[[run]]
            full <- metrosamp::concat(base, contin)
            stat <- proc_mc_helper(full, codasize)
            runstats[[run]] <- stat[[1]]
            mcobjs[[run]] <- full
            codaobjs[[run]] <- stat[[2]]
        }
        else {
            ## copy any runs that don't have continuations
            runstats[[run]] <- bruns$runstats[[run]]
            mcobjs[[run]] <- bruns$mcobjs[[run]]
            codaobjs[[run]] <- bruns$codaobjs[[run]]
        }
    }

    if(!is.null(output_dir)) {
        repatrn <- "[^-]+-[0-9]+-([0-9]+)-mcrslt.rds"
        for(mcobj in mcobjs) {
            for(run in names(mcobj)) {
                runobj <- mcobj[[run]]
                nsamp <- nrow(runobj$samples)
                match <- regexec(repatrn, run)
                runid <- regmatches(run, match)[[1]][2]
                filename <- file.path(output_dir, paste(name_stem, nsamp, runid, 'mcrslt.rds', sep='-'))
                saveRDS(runobj, filename, compress='xz')
            }
        }
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
    if(nsamp > 1000) {
        fac <- ceiling(nsamp/1000)
        pltsamps <- samps[idx%%fac == 0, ]
    }
    else {
        pltsamps <- samps
    }
    pltsamps <- as.data.frame(pltsamps)
    GGally::ggpairs(data=pltsamps)
}
