
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
#' @param climate_data an object of climate data, must contain value, year, experiment, and variable columns.
#' It may not contain data for more years, variables, or experiments then in the principal_components argument.
#' @param principal_components a \code{hectorcal} object that contains results from \code{prcomp} and meta data information.
#' @retun a list of three objects, rotation the PC loadings, scale the values used to scale the input climate data, and center the values used to center the climate data.
#' @importFrom dplyr %>%
#' @export
project_climate <- function(climate_data, principal_components){

    # First check to make sure that the climate data contains the required columns
    check_columns(climate_data, c('value', 'year', 'experiment', 'variable'))

    # Return the output strucutre
    obj <- list()

    # Check to see that the climate data being read into the principal componets does not contain data the
    # principal componet strucutre does not have data for.
    missing_years      <- difference_between(climate_data$year, principal_components$meta_data$year)
    missing_variable   <- difference_between(climate_data$variable, principal_components$meta_data$variable)
    missing_experiment <- difference_between(climate_data$experiment, principal_components$meta_data$experiment)

    if(length(missing_years)){stop('climate_data and principal_components are required to have the same years')}
    if(length(missing_variable)){stop('climate_data and principal_components are required to have the same variables')}
    if(length(missing_experiment)){stop('climate_data and principal_components are required to have the same experiments')}

    # Match the cliamte data with the pc center and scales, this ensures that the climate data will be centered and
    # scaled properly in the following step.
    climate_data %>%
        dplyr::mutate(col = paste0('X', year, '_', variable, '_', experiment)) %>%
        dplyr::full_join(tibble::tibble(center = principal_components$center,
                                        col = names(principal_components$center)), by = 'col') %>%
        dplyr::full_join(tibble::tibble(scale = principal_components$scale,
                                        col = names(principal_components$scale)), by = 'col') %>%
        dplyr::full_join(tibble::tibble(col = row.names(principal_components$rotation),
                                        keep = 1), by = 'col')->
        tibble_value_scale_center

    # If there are any NAs throw an error, we expect there to be a 1:1 match.
    if(any(is.na(tibble_value_scale_center))) { stop('There is a mismatch between the information in the climate_data and the principal_components arguments.') }

    # Scale and center the new data
    newdata <- (tibble_value_scale_center$value - tibble_value_scale_center$center) / tibble_value_scale_center$scale
    obj[['rotation']] <- newdata %*% principal_components$rotation

    # Scales
    scale             <- tibble_value_scale_center$scale
    names(scale)      <- tibble_value_scale_center$col
    obj[['scale']]    <- scale

    # Center
    center            <- tibble_value_scale_center$center
    names(center)     <- tibble_value_scale_center$col
    obj[['center']]   <- center

    # Return the object
    obj
}


#' Reconstruct climate data from projected climate data and prcinicpal componetns
#'
#' @param projected_climate an object returned by \code{project_climate}
#' @param principal_components a \code{hectorcal} object that contains results from \code{prcomp} and meta data information.
#' @retun a data frame of reconstructed cliamte data containing the following columns, year, value, experiment, and variable
#' @import dplyr
#' @import tidyr
#' @export
reconstruct_climate <- function(projected_climate, principal_components){

    # Revert the information back to original coordinate system
    uncentered_unscaled_data <- colSums(t(principal_components$rotation) * as.vector(projected_climate$rotation))

    # Recenter and scale the data
    value <- c(uncentered_unscaled_data * principal_components$scale) + principal_components$center

    # Format information in data frame
    data.frame(value = as.vector(value), col = names(value)) %>%
        tidyr::separate(col, into = c('year', 'variable', 'experiment')) %>%
        dplyr::mutate(year = as.integer(gsub('X', '', year))) %>%
        dplyr::select(year, value, variable, experiment)

}


