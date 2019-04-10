#### Functions for dealing with the output gates

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


#' Find the runids that passed the gates in all experiments
#'
#' @param hdata Data frame of hector ensemble runs, including all of the
#' experiments to use in the filtering
#' @param vars Names of variables to include in the filtering
#' @return Vector of runids that passed all gates in all experiments.
#' @examples
#' ids <- chkgates(rbind(hector_conc_ensemble$rcp26, hector_conc_ensemble$rcp85),
#'                'tas')
#' length(ids)    # == 183
#' @export
chkgates <- function(hdata, vars) {
    variable <- year <- pass <- mina <- maxb <- value <- NULL

    expts <- unique(hdata$experiment)

    gates <- lapply(expts, function(expt) {
        lapply(vars, function(v) {
            get_gates(expt, v)
            }) %>%
            dplyr::bind_rows()}) %>% dplyr::bind_rows()

    gatecmp <- dplyr::left_join(gates, hdata, by=c('variable','year','experiment')) %>%
        dplyr::mutate(pass = mina <= value & value <= maxb)
    gatepass <- dplyr::group_by(gatecmp, runid) %>% dplyr::summarise(allpass = all(pass)) %>%
        dplyr::ungroup()
    gatepass$runid[gatepass$allpass]
}
