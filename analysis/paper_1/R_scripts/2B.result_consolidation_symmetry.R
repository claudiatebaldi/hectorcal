## This script consolidates all of the symetry run results into a csv file that is easy to
## work with on in the plotting scripts.

# 0. Set Up --------------------------------------------------------------------------------------
# This should be set to the hectorcal/analysis/paper_1 as the working directory
BASE_DIR <- getwd()
OUTPUT_DIR <- file.path(BASE_DIR, 'output')
OUTPUT_DIR_2A <- file.path(OUTPUT_DIR, '2A.symmetry_test_results')

# Load the required libraries
library(dplyr)
library(tidyr)
library(hector)


# 1. Format and Consolidate the Concentration Runs -----------------------------------------------
list.files(OUTPUT_DIR_2A, pattern = 'conc', full.names = TRUE, recursive = TRUE) %>%
    lapply(function(input){
        message(input)

        object <- readRDS(input)
        if( object$convergence != 0){

            message(paste0('problem with: ', input))

        } else {

            # Parse out the part of the file name that will tell us about the calbration
            # protocol that was used. This name will be used to add calibration method
            # and model name to the data frame.
            unique_name <- gsub(x = input, pattern = OUTPUT_DIR_2A, replacement = '')

            # Format the model, method, fitted parmeter values, and optmized value into
            # a data frame.
            data.frame(run = gsub(pattern = '.rds', replacement = '', x = basename(unique_name)),
                       sym_method = gsub(pattern = '/', replacement = '', dirname(unique_name)),
                       stringsAsFactors = FALSE) %>%
                cbind(t(object$par)) %>%
                mutate(optmized_value = object$value,
                       convergence = object$convergence,
                       diff = object$diff)

        }


    }) %>%
    bind_rows() %>%
    # Format the names of the models
    separate(run, c('driven', 'model', 'x'), sep = '_') %>%
    select(-x) ->
    rslts_conc


# 1. Format and Consolidate the Emission Driven Runs ---------------------------------------------
list.files(OUTPUT_DIR_2A, pattern = 'emiss', full.names = TRUE, recursive = TRUE) %>%
    lapply(function(input){
        message(input)
        object <- readRDS(input)

        if( object$convergence != 0){

            message(paste0('problem with: ', input))

        } else {

            # Parse out the part of the file name that will tell us about the calbration
            # protocol that was used. This name will be used to add calibration method
            # and model name to the data frame.
            unique_name <- gsub(x = input, pattern = OUTPUT_DIR_2A, replacement = '')

            # Format the model, method, fitted parmeter values, and optmized value into
            # a data frame.
            data.frame(run = gsub(pattern = '.rds', replacement = '', x = basename(unique_name)),
                       sym_method = gsub(pattern = '/', replacement = '', dirname(unique_name)),
                       stringsAsFactors = FALSE) %>%
                cbind(t(object$par)) %>%
                mutate(optmized_value = object$value,
                       convergence = object$convergence,
                       beta = object$beta)

        }


    }) %>%
    bind_rows() %>%
    # Format the names of the models
    separate(run, c('driven', 'model', 'x'), sep = '_') %>%
    select(-x) %>%
    filter(beta <= 1.5) ->
    rslts_emiss


# Save the results

rslts <- bind_rows(rslts_conc, rslts_emiss)
write.csv(rslts, file = file.path(OUTPUT_DIR, 'sym_results.csv'), row.names = FALSE)



