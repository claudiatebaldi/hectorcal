# Calculate the PC for the ESM data.

# 0. Set Up ------------------------------------------------------------------------
library(hectorcal)
library(tibble)
library(dplyr)
library(tidyr)
library(devtools)

# 1. Define Functions --------------------------------------------------------------

# Perform the PCA on Hector ensemble data to create the Hector PC external data
#
# Args
#   data_list: a list of the PCA_hector_emsble data to use in the PCA, this funciton can handle data from muliplte scenrios
#   combos: optional agument that if set to TRUE will perform the PCA on of the possible different scenario combinations.
# Returns: a prcomp object with a meta data list of infomration about the years, variables, and scenarios that went into the PCA.
get_hector_pcs <- function(data_list, combos = FALSE){

    # Concatenate the list of the data frames into a singl df
    rslts_long <- dplyr::bind_rows(data_list)
    scn_names  <- unique(rslts_long$scenario)

    # Create a matrix of the scenarios with an indicator 0 and 1 that the
    possible_scns <- as_tibble(matrix(rep(c(1, 0), length(scn_names)),
                                      ncol = length(scn_names),
                                      dimnames = list(NULL, scn_names)))

    if(combos){
        # Create all list of all the possible scenario combinations, the
        # list of the scenarios will be used to subset the hector ensemble.
        expand.grid(possible_scns) %>%
            mutate(combo = 1:nrow(.)) %>%
            gather(scenario, keep, -combo) %>%
            arrange(combo) %>%
            filter(keep == 1) %>%
            split(.$combo) ->
            scn_combo_list
    } else {

        scn_combo_list <- list(`1` = distinct(select(rslts_long, scenario, variable)))

    }


    # For each possible combination of scenarios subset the Hector ensemble and save
    # it as external package data.
    scn_combo_list %>%
        lapply(function(input){

            # Construct a list of meta data infomration that is going to be returned with the prcomp
            # object and will be used to compare with the climate data if the PC data contains
            # the correct information to be compared with or not.
            meta_data <- list()
            meta_data$year       <- unique(rslts_long$year)
            meta_data$variable   <- unique(rslts_long$variable)
            meta_data$experiment <- unique(rslts_long$experiment)
            meta_data$scenario   <- unique(input[['scenario']])

            # Subset the Hector data so that it only contains results for the scneario of intrest.
            rslts_long %>%
                filter(scenario %in% input[['scenario']]) %>%
                # Name each column after the year, variable, and scenario name
                mutate(col = paste0('X', year, '_', variable, '_', experiment)) %>%
                select(runid, col, value) %>%
                # Format as a wide data frame in prep for the prcomp call
                spread(col, value) %>%
                na.omit() %>%
                select(-runid) ->
                df

            # PCA
            prcomp_obj <- prcomp(df, center = TRUE, scale. = TRUE, retx = FALSE)

            # Return the prcmop output with the meta data information
            prcomp_obj[['meta_data']] <- meta_data

            # Name the prcomp object to save a external data
            name <- paste0('PC_hector_', paste(unique(input[['scenario']]), collapse = '_'))
            assign(name, prcomp_obj)

            # Save the data frame as external package data
            do.call("use_data", list(as.name(name), overwrite = TRUE))

        })
}


# 2. Calculte Hector ensmbel PCs ---------------------------------------------

# Concentration runs
get_hector_pcs(list(`PCA_hector_ensemble-concen-rcp26`, `PCA_hector_ensemble-concen-rcp45`,
                    `PCA_hector_ensemble-concen-rcp60`, `PCA_hector_ensemble-concen-rcp85`), combos = FALSE)

# Emissions driven runs with varying carbon cycle parameters
get_hector_pcs(list(`PCA_hector_ensemble-emissCC-esmrcp26`, `PCA_hector_ensemble-emissCC-esmrcp45`,
                    `PCA_hector_ensemble-emissCC-esmrcp60`, `PCA_hector_ensemble-emissCC-esmrcp85`), combos = FALSE)

# Emission Driven Runs with constant carbon
get_hector_pcs(list(`PCA_hector_ensemble-emissConstantC-esmrcp26`, `PCA_hector_ensemble-emissConstantC-esmrcp45`,
                    `PCA_hector_ensemble-emissConstantC-esmrcp60`, `PCA_hector_ensemble-emissConstantC-esmrcp85`), combos = TRUE)


