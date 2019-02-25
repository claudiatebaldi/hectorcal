
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
#' @retun a list of three objects, rotation (the PC loadings), scale (the values used to scale the input climate data), and
#' center (the values used to center the climate data).
#' @importFrom dplyr %>%
#' @export
project_climate <- function(PCA_climate, principal_components){

    # First check to make sure that the climate data contains the required columns
    check_columns(PCA_climate, c('value', 'year', 'experiment', 'variable'))

    # Return the output strucutre
    obj <- list()

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
        tibble_value_scale_center


    # Calculate the dot product of the scaled & centered data by the principal componet loadings.
    pc_loadings <- principal_components$rotation[row.names( principal_components$rotation) %in% tibble_value_scale_center$col, ]

    # Scale and center the new data
    n_pc <- ncol(pc_loadings)
    center  <- tibble_value_scale_center$center
    scale   <- tibble_value_scale_center$scale
    newdata <- matrix(rep(tibble_value_scale_center$value, times = n_pc), ncol = length(center))
    #newdata <- scale(matrix(tibble_value_scale_center$value, ncol = length(center)), center, scale)
   obj[['rotation']] <- scale(newdata, center, scale) %*% pc_loadings

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
#' @importFrom dplyr %>%
#' @export
reconstruct_climate <- function(projected_climate, principal_components){

    matrix <- projected_climate$rotation %*% t(principal_components$rotation)
    value  <- matrix * projected_climate$scale + projected_rslt$center

    data.frame(value = as.vector(value), col = colnames(value)) %>%
        tidyr::separate(col, into = c('year', 'variable', 'experiment')) %>%
        dplyr::mutate(year = as.integer(gsub('X', '', year))) %>%
        dplyr::select(year, value, variable, experiment)

}


