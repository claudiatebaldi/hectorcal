
#' Check a dataframe or tibble for required columns, if missing columns throw an error
#'
#' @param input the data frame or the tibble with the column names to check.
#' @param req_cols a vector of the required column name strings
check_columns <- function(input, req_cols){

    missing <- (!req_cols %in% names(input))
    if(any(missing)){
        stop('Missing columns: ', paste(req_cols[missing], collapse = ', '))
    }
}



#' Project a data frame of climate data on to a principal componet basis
#'
#' @param climate_data a data frame of climate data.  The data must contain
#' value, year, experiment, and variable columns.
#' @param principal_components a \code{hectorcal} object that contains results
#' from \code{prcomp} and meta data information.
#' @param row_vector If \code{TRUE}, return the data as a row vector (a 1xN
#' matrix)
#' @return A vector of projection coefficients for the climate data.
#' @importFrom dplyr %>%
#' @importFrom assertthat assert_that
#' @export
project_climate <- function(climate_data, principal_components, row_vector=TRUE)
{

    experiment <- variable <- year <- NULL

    # First check to make sure that the climate data contains the required columns
    check_columns(climate_data, c('value', 'year', 'experiment', 'variable'))

    # Return the output strucutre
    obj <- list()

    ## Check to see that we have all the experiments we need.  We make this
    ## check explicitly because getting all the experiments into a single data
    ## frame probably requires assembling the data from several datasets and is
    ## therefore easy to screw up.  For years and variables we assume that the
    ## datasets are fairly uniform and therefore likely to have all the data
    ## needed by default.
    assert_that(all(principal_components$meta_data$experiment %in% climate_data$experiment))

    ## Filter the data to include just the required data, and put the rows into
    ## canonical order
    cd <- dplyr::filter(climate_data,
                        experiment %in% principal_components$meta_data$experiment,
                        variable %in% principal_components$meta_data$variable)
    ishist <- grepl('[Hh]istorical', cd$experiment)

    ## Handle the filtering by year separately for historical and future and for each variable,
    ## since some of the CMIP inputs have overlapping years for the two types, and each
    ## variable might have different years.
    cd <-
      dplyr::bind_rows(
        lapply(principal_components$meta_data$variable,
               function(var) {
                 dplyr::bind_rows(
                   dplyr::filter(cd, variable==var,
                                 ishist & year %in% principal_components$meta_data$histyear[[var]]),
                   dplyr::filter(cd, variable==var,
                                 !ishist & year %in% principal_components$meta_data$year[[var]]))
               })) %>%
      dplyr::arrange(experiment, variable, year)

    assert_that(nrow(cd) == nrow(principal_components$rotation))

    ## Extract the data, center, and scale
    vp <- (matrix(cd$value, nrow=1) - principal_components$center) / principal_components$scale

    p <- vp %*% principal_components$rotation

    if(row_vector) {
        p
    }
    else {
        as.vector(p)
    }
}


#' Reconstruct climate data from projected climate data and prcinicpal componetns
#'
#' @param projected_climate a vector or row-vector of PC coordinates
#' @param principal_components a \code{hectorcal} object that contains results
#' from \code{prcomp} and meta data information.
#' @param ncomp Number of components to keep in the reconstruction.  Default is
#' to keep them all.
#' @return a data frame of reconstructed cliamte data containing the following
#' columns, year, value, experiment, and variable
#' @importFrom dplyr %>%
#' @export
reconstruct_climate <- function(projected_climate, principal_components, ncomp=NA)
{
    experiment <- variable <- value <- year <- NULL

    if(is.na(ncomp)) {
        ncomp <- ncol(principal_components$rotation)
    }

    if(is.vector(projected_climate)) {
        projected_climate <- matrix(projected_climate, nrow=1)
    }

    p <- projected_climate[ , 1:ncomp, drop=FALSE]
    rinv <- t(principal_components$rotation)[1:ncomp, , drop=FALSE]

    vp <- p %*% rinv

    v <- vp*principal_components$scale + principal_components$center

    # Format information in data frame
    data.frame(value = as.vector(v), col = colnames(rinv)) %>%
        tidyr::separate(col, into = c('experiment', 'variable', 'year'), sep='\\.') %>%
        dplyr::mutate(year = as.integer(year)) %>%
        dplyr::select(year, value, variable, experiment)
}

## This is a helper function for compute_pc below
## Given a vector of years that are meant to apply to all variables, construct
## an explicit named list of year vectors for each variable
##
## Alternatively, if years is already a named list, validate that it has the right
## names and return it if it does.
chk_yrlist <- function(variables, years)
{
  if(is.list(years)) {
    assertthat::assert_that(setequal(names(years), variables))
    years
  }
  else {
    yrlist <- vector('list', length(variables))
    names(yrlist) <- variables
    for(i in seq_along(yrlist)) {
      yrlist[[i]] <- years
    }
    yrlist
  }
}

## Another helper function for compute_pc
## Check that for each variable all years are present
chk_year <- function(df, years, histyears)
{
  spldf <- split(df, df$variable)
  yrchk <- lapply(names(spldf),
                  function(var) {
                    dplyr::group_by(spldf[[var]], experiment) %>%
                      dplyr::summarise(complete =
                                         all(years[[var]] %in% year) |
                                         all(histyears[[var]] %in% year))
                  })
  vargood <- sapply(yrchk, function(df) {all(df$complete)})
  all(vargood)
}

#' Compute principal components from an ensemble of Hector runs
#'
#' Given a collection of Hector runs, filter to the desired experiments,
#' variables, years, and future years and compute the principal components for
#' that combination.
#'
#' This function will work with the package datasets
#' \code{\link{hector_conc_ensemble}} and \code{\link{hector_emiss_ensemble}}
#' as the \code{scenlist} argument.
#'
#' On the ordering of variables when constructing principal components: When
#' flattening a set of runs into a vector, the grouping order shall be
#' experiment, variable, year.  That is, year will vary the most rapidly,
#' experiment the least.  This ordering is assumed in all of the other functions
#' that deal with principal components.
#'
#' @param scenlist List of data frames of hector runs, one for each experiment.
#' The names of the list elements should correspond to the experiments contained
#' in each.
#' @param experiments Vector of strings denoting the experiments to include.
#' Elements of scenlist not matching one of these strings will be dropped.
#' @param variables Vector of Hector output variables to include in the analysis.
#' @param years Vector of years to include in experiments using \emph{future}
#' runs.  Generally these will be the various RCP pathways.  If different variables
#' have different years, then this argument may be a named list of the form
#' \code{list(var1=var1years, var2=var2years, ...)}
#' @param histyears Vector of years to include in experiments using \emph{past}
#' runs.  Any experiment with ``historical'' or ``Historical'' in the name is
#' considered to be a past run.  Different groups of years may be specified for
#' different variables as described for the \code{years} parameter.
#' @param retx Flag indicating whether to retain the projections of the ensemble
#' runs onto the principal components.
#' @importFrom assertthat assert_that
#' @export
compute_pc <- function(scenlist, experiments, variables, years, histyears, retx=FALSE)
{
    variable <- year <- experiment <- runid <- NULL # silence warnings about NSE

    chk_yrlist(variables, years)
    chk_yrlist(variables, histyears)

    ## Drop the scenarios that we won't be using (this allows us to pass the same
    ## list every time)
    if(!is.null(names(scenlist))) {
        ## Use the names as a guide to the experiments
        expts <- names(scenlist)
    }
    else {
        ## get the experiment names from the tables
        expts <- sapply(scenlist, function(d) {
                            expt <- unique(d$experiment)
                            assert_that(length(expt)==1)
                            expt
                        })
    }
    scenlist <- scenlist[expts %in% experiments]
    assert_that(length(scenlist) == length(experiments)) # all the requested
                                        # experiments are represented.

    ## combine into a single dataset
    alldata <- dplyr::bind_rows(scenlist)

    ## Drop any variables and years we will not be using
    ishist <- grepl('[Hh]istorical', alldata$experiment)
    alldata <-
      dplyr::bind_rows(
        lapply(variables,
               function(var) {
                 dplyr::filter(alldata, variable == var,
                               (ishist & year %in% histyears[[var]]) |
                                 (!ishist & year %in% years[[var]]))
               }))

    ## Check that all the requested variables and years are present for all
    ## experiments.  Technically, we should group by runid too, but that will be
    ## really slow, so we take it as given that all the runs were generated by
    ## the same process and will therefore have the same years, variables, etc.
    varchk <- dplyr::group_by(alldata, experiment) %>%
      dplyr::summarise(complete=all(variables %in% variable))
    assert_that(all(varchk$complete))

    yrchk <- chk_year(alldata, years, histyears)
    assert_that(yrchk==TRUE)

    ## Now arrange the data in the order described above and split by runid
    alldata <- dplyr::arrange(alldata, runid, experiment, variable, year)
    alldata_l <- split(alldata, alldata$runid)

    ## Produce a matrix of our output values.  Runs in rows, and data items in
    ## columns.
    data_matrix <- do.call(rbind,
                           lapply(alldata_l, function(d){d$value}))


    ## Add some column names.
    proto <- alldata_l[[1]]
    cnames <- paste(proto$experiment, proto$variable, proto$year, sep='.')
    colnames(data_matrix) <- cnames

    ## Perform the PCA
    pca <- stats::prcomp(data_matrix, center=TRUE, scale.=TRUE, retx=retx)

    ## Add our metadata
    srtyears <- lapply(years, sort)
    srthistyears <- lapply(histyears, sort)
    pca$meta_data <- list(year=srtyears,
                          histyear=srthistyears,
                          variable=sort(variables),
                          experiment=sort(experiments),
                          scenario=sort(unique(alldata$scenario)))

    pca
}



#' Create comparison data for principal components decomposition
#'
#' We have stored the principal components decomposition for each CMIP5 model as
#' package data.  This function summarises that into a format similar to the
#' \code{\link{esm_comparison}} dataset.
#'
#' @param pcdata Input principal components decomposition (probably either
#' \code{\link{cmip_conc_pcproj}} or \code{\link{cmip_emiss_pcproj}}).
#' @param pcmax Maximum PC to use.  PCs higher than this will be dropped.
#' @param pcselect Individual PCs to select.  This is intersected with the list
#' implied by pcmax.
#' @importFrom dplyr %>%
#' @export
summarize_pcdecomp <- function(pcdecomp, pcmax=NULL, pcselect=NULL)
{
    PC <- NULL                          # bindings for NSE vars
    if(!is.null(pcmax)) {
        pcdecomp <- dplyr::filter(pcdecomp, PC <= pcmax)
    }
    if(!is.null(pcselect)) {
        pcdecomp <- dplyr::filter(pcdecomp, PC %in% pcselect)
    }
    dplyr::group_by(pcdecomp, PC) %>%
      dplyr::summarise(mina=min(value), maxb=max(value),
                       a10=quantile(value, 0.1, names=FALSE),
                       b90=quantile(value, 0.9, names=FALSE),
                       cmean=mean(value), cmedian=median(value)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(PC)
}

#' Update PCA structure metadata to the new format
#'
#' The new format for PCA structures has lists instead of vectors for the years and
#' histyears elements.  Use this function when you have a structure saved in the old
#' format and you need to convert it to the new.
#'
#' @param pcastruct Old format PCA structure
#' @return New format PCA structure
#' @export
update_pca_struct <- function(pcastruct)
{
  pcastruct$meta_data$year <- chk_yrlist(pcastruct$meta_data$variable,
                                          pcastruct$meta_data$year)
  pcastruct$meta_data$histyear <- chk_yrlist(pcastruct$meta_data$variable,
                                             pcastruct$meta_data$histyear)
  pcastruct
}
