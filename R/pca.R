
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

#' Determine the elements that are different between two vectors.
#'
#' @param x a vector
#' @param y a second vector
#' @return a vector containing the elements that are unique to the two input vectors
difference_between <- function(x, y){

    c(setdiff(x, y), setdiff(y, x))

}


#' Project a data frame of climate data on to a principal componet basis
#'
#' @param climate_data a data frame of climate data.  The data must contain
#' value, year, experiment, and variable columns.
#' @param principal_components a \code{hectorcal} object that contains results
#' from \code{prcomp} and meta data information.
#' @param row_vector If \code{TRUE}, return the data as a row vector (a 1xN
#' matrix)
#' @return a list of three objects, rotation the PC loadings, scale the values
#' used to scale the input climate data, and center the values used to center
#' the climate data.
#' @importFrom dplyr %>%
#' @importFrom assertthat assert_that
#' @export
project_climate <- function(climate_data, principal_components, row_vector=TRUE)
{

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
                        variable %in% principal_components$meta_data$variable,
                        year %in% principal_components$meta_data$year) %>%
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

