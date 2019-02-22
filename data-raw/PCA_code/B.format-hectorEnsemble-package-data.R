# Format Hector results into the correct structure for the PCA and save a external package data.

# 0. Set Up  ----------------------------------------------------------------------
# Required libraries
library(dplyr)
library(tidyr)
library(devtools)

# The size desired size of the ensemble, the number of Hector runs to keep.
ensemble_size <- 1000

# 1. Define Function ------------------------------------------------------------------

# Format and save the Hector results in the standard hectorcal package pca input data format
#
# Args
#   file_flag: a string that is used to describe a set of files in the data-raw/PCA_pic_results directory
#   ensemble_size: the number of Hector runs to include in the package data
# Returns: a matrix of hector runs as rows and years as columns for the scenarios the rda is named after
format_hector_ensemble_package_data <- function(file_flag, ensemble_size){

    # Import files associated with a particular experiment ensemble, such as concentration runs,
    # emission runs, ect and concatenate the ensemble into a single long data frame.
    suppressWarnings(lapply(list.files(file.path('data-raw', 'PCA_pic_results'), file_flag,
                                       full.names = TRUE), readRDS) %>%
                         bind_rows(.) ->
                         rlsts_long)

    # Save a copy of the parameters, this information will be added back to the data at the end.
    rlsts_long %>%
        select(-value, -variable, -year) %>%
        distinct %>%
        na.omit ->
        param_mapping

    # Make sure that there are values for every run for every sceanario.
    rlsts_long %>%
        # Rename the variable to match the cmip notation
        dplyr::mutate(variable = if_else(variable == 'Tgav', 'tas', 'co2')) %>%
        dplyr::mutate(col = paste0(variable, '_', scenario, '_', year)) %>%
        dplyr::select(runid, value, col) %>%
        na.omit %>%
        tidyr::spread(col, value) %>%
        na.omit ->
        rlsts_wide

    # Subset the wide dataframe so that it contains the number of Hector cases for the required ensemble
    rlsts_wide[1:ensemble_size, ] %>%
        # Discard the NAs so now the tibble only contains values for the runs that are consistent across the
        # different RCPs
        na.omit %>%
        tidyr::gather(col, value, -runid) %>%
        # Seperatre the run informatation out from the column name
        tidyr::separate(col, into = c('variable', 'scenario', 'year'), sep = '_') %>%
        # Add an rcp column
        tidyr::separate(scenario, into = c('X', 'experiment'), sep = '-', remove = FALSE) %>%
        dplyr::select(-X) %>%
        # Add the paramter values
        dplyr::inner_join(param_mapping, by = c('runid', 'scenario')) ->
        cleaned_data

    # For each rcp save results as external pacakge data
    cleaned_data %>%
        split(.$scenario) %>%
        lapply(function(input){

            # Subset the data frame to only contain the number of rows of the
            # as the defined ensemble size and name the data frame.
            name <- paste0('PCA_hector_ensemble-', unique(input$scenario))
            assign(name, input)

            # Save the data frame as external package data
            do.call("use_data", list(as.name(name), overwrite = TRUE))

        })

}

# 2. Process the Hector Data -----------------------------------------------------

format_hector_ensemble_package_data('concen-', ensemble_size)
format_hector_ensemble_package_data('emissCC-', ensemble_size)
format_hector_ensemble_package_data('emissConstantC-', ensemble_size)


# End

