#### Misc useful functions

#' Calculate the fraction of the variance of PCA restlts
#'
#' @param pca object returned by \code{\link[stats]{prcomp}}
#' @return Vector of of the portion of variance expalained by the PCs.
#' @seealso \code{\link{plot_varfrac}}
#' @export
calc_variance <- function(pca){
    cumsum(pca$sdev^2)/sum(pca$sdev^2)
}


#' Convert a list of metrosamp structures to a \code{\link[coda]{mcmc.list}}
#' structure.
#'
#' @param mslist A list mof metrosamp structures
#' @export
metrosamp2coda <- function(mslist) {

    coda::mcmc.list(
        lapply(mslist, function(ms) {
                   coda::mcmc(ms$samples)
               })
        )
}


#' Get a data frame of gates for an experiment
#'
#' Find the minimum and maximum attested values by year for a climate variable
#' in the CMIP5 output.  Values are returned at the beginning, end, and median
#' of the requested period.
#'
#' For historical experiments, the final year can be no later than 2005.  For
#' future experiments, the starting year can be no earlier than 2006.
#' The \code{minyear} and \code{maxyear} parameters will be adjusted
#' automatically to these limits.
#'
#' @param expt Experiment to retrieve ('historical', 'rcp26', 'rcp45', etc.)
#' @param var Variable to retrieve ('tas' or 'co2')
#' @param minyear Earliest year to retrieve.
#' @param maxyear Latest year to retrieve.
#' @return Data frame with the matching rows from the
#' \code{\link{esm_comparison}} dataset.
#' @export
get_gates <- function(expt, var, minyear=1861, maxyear=2100) {
    if(grepl('[Hh]istorical', expt)) {
        maxyear <- min(c(2005, maxyear))
    }
    else {
        minyear <- max(c(2006, minyear))
    }

    exptdata <- filter(esm_comparison, experiment==expt, variable==var,
                       year <= maxyear, year >= minyear)
    yrs <- c(min(exptdata$year),
             max(exptdata$year[exptdata$year <= median(exptdata$year)]),
             max(exptdata$year))
    filter(exptdata, year %in% yrs)
}
