
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
#' @param PCA_climate an object of climate data, must contain value, year, experiment, and variable columns.
#' It may not contain data for more years, variables, or experiments then in the principal_components argument.
#' @param principal_components a \code{hectorcal} object that contains results from \code{prcomp} and meta data information.
#' @retun a list of three objects \item{rotation}{the PC loadings}, \item{scale}{the values used to scale the input climate data}, and
#' \item{center}{the values used to center the climate data}
#' @export
# DO we want to add the sdev column? Do we want to add the meta data?
project_climate <- function(PCA_climate, principal_components){

    # First check to make sure that the climate data contains the required columns
    check_columns(PCA_climate, c('value', 'year', 'experiment', 'variable'))

    # Check to see that the climate data being read into the principal componets does not contain data the
    # principal componet strucutre does not have data for.
    missing_years      <- setdiff(PCA_climate$year, principal_components$meta_data$year)
    missing_variable   <- setdiff(PCA_climate$variable, principal_components$meta_data$variable)
    missing_experiment <- setdiff(PCA_climate$experiment, principal_components$meta_data$experiment)

    if(length(missing_years)){stop('The principal_components object is missing data for years: ', paste(missing_years, collapse = ', '))}
    if(length(missing_variable)){stop('The principal_components object is missing data for varaible: ', missing_variable)}
    if(length(missing_experiment)){stop('The principal_components object is missing data for experiment: ', missing_experiment)}

    # Make sure that the climate data is not missing experiment or variable data.
    missing_varaible   <- setdiff(principal_components$meta_data$variable, PCA_climate$variable)
    missing_experiment <- setdiff(principal_components$meta_data$experiment, PCA_climate$experiment)

    if(length(missing_variable)){stop('The PCA_climate object is missing data for varaible: ', missing_variable)}
    if(length(missing_experiment)){stop('The PCA_climate object is missing data for experiment: ', missing_experiment)}

    # Center and scale the climate data.
    PCA_climate %>%
        dplyr::mutate(col = paste0('X', year, '_', variable, '_', experiment)) %>%
        dplyr::full_join(tibble::tibble(center = principal_components$center,
                                        col = names(principal_components$center)), by = 'col') %>%
        dplyr::full_join(tibble::tibble(scale = principal_components$scale,
                                        col = names(principal_components$scale)), by = 'col') %>%
        na.omit ->
        inter_df

    inter_df %>%
        dplyr::mutate(value = (value - center) / scale) %>%
        select(col, value) %>%
        distinct  ->
        centered_scaled_values

    # Subset the principal component loadings to match the contents of the cliamte data
    pc_loadings <- principal_components$rotation[row.names( principal_components$rotation) %in% centered_scaled_values$col, ]

    # Multiply the cenered adn scaled values by the
    rotation <- pc_loadings * centered_scaled_values$value

    # Return the output strucutre
    obj <- list()
    # Multiply the cenered and scaled loadings by the value
    obj[['rotation']] <- pc_loadings * centered_scaled_values$value

    # Scales
    scale        <- c(inter_df$scale)
    names(scale) <- inter_df$col
    obj[['scale']]   <- scale

    # Center
    center           <- inter_df$center
    names(center)    <- inter_df$center
    obj[['center']]  <- center

    # Return the object
    obj
}

# Not sure about this funciton...
revert_climate <- function(projected_climate){

    sum_projection <- apply(projected_climate$rotation, 1, sum)
    tibble(value = sum_projection,
           names = names(sum_projection)) %>%
        left_join(tibble(scale = projected_climate$scale,
                         names = names(projected_climate$scale)),
                  by = 'names') %>%
        left_join(tibble(scale = projected_climate$scale,
                         names = names(projected_climate$scale)),
                  by = 'names')

}


