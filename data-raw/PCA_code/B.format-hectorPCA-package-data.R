# Format Hector results into the correct structure for the PCA and save a external package data.

# 0. Set Up  ----------------------------------------------------------------------
# Required libraries
library(dplyr)
library(tidyr)
library(devtools)

# The size of the ensemble -- will want to change to be 1000
ensemble_size <- 2

# 1. Define Function ------------------------------------------------------------------

# Format and save the Hector results in the standard hectorcal package pca input data format
#
# Args
#   file_flag: a string that is used to describe a set of files in the data-raw/PCA_pic_results directory
#   ensemble_size: the number of Hector runs to include in the package data
# Returns: a matrix of hector runs as rows and years as columns for the scenarios the rda is named after
format_hector_pca_package_data <- function(file_flag, ensemble_size){

    # Import files associated with a particular experiment ensemble, such as concentration runs,
    # emission runs, ect and concatenate the ensemble into a single long data frame.
    suppressWarnings(lapply(list.files(file.path('data-raw', 'PCA_pic_results'), file_flag,
                                       full.names = TRUE), function(input){readRDS(input)[[1]]}) %>%
                         bind_rows() ->
                         rlsts_long)

    # Save a copy of all the scenario names
    scn_names <- unique(rlsts_long$scenario)

    # Create a matrix of the scenarios with an indicator 0 and 1 that the
    # scenario can be in a combination or not.
    possible_scns <- as_tibble(matrix(rep(c(1, 0), length(scn_names)),
                                      ncol = length(scn_names),
                                      dimnames = list(NULL, scn_names)))

    # Create all list of all the possible scenario combinations, the
    # list of the scenarios will be used to subset the hector ensemble.
    expand.grid(possible_scns) %>%
        mutate(combo = 1:nrow(.)) %>%
        gather(scn, keep, -combo) %>%
        arrange(combo) %>%
        filter(keep == 1) %>%
        split(.$combo) ->
        combination_list

    # For each possible combination of scenarios subset the Hector ensemble and save
    # it as external package data.
    combination_list %>%
        lapply(function(input){

            # Subset the Hector data
            rlsts_long %>%
                filter(scenario %in% input$scn) %>%
                # Name each column after the year, variable, and scenario name
                mutate(col = paste0('X', year, '_', variable, '_', scenario)) %>%
                select(runid, col, value) %>%
                # Format as a wide data frame in prep for the prcomp call
                spread(col, value) %>%
                na.omit() %>%
                select(-runid) ->
                df

            # Check to make sure that the
            if(nrow(df) < ensemble_size) stop('Hector ensemble is too small')

            # Subset the data frame to only contain the number of rows of the
            # as the defined ensemble size and name the data frame.
            name <- paste(input$scn, collapse = '_')
            assign(name, df[1:ensemble_size, ])

            # Save the data frame as external package data
            do.call("use_data", list(as.name(name), overwrite = TRUE))
        })
}

# 2. Process the Hector Data -----------------------------------------------------

format_hector_pca_package_data('concen-', ensemble_size)
format_hector_pca_package_data('emiss-CC-', ensemble_size)

# End

