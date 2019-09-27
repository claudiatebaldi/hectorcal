## Hector calibration paper 1
## Since we observed symmetries about the climate parameters we should
## also check for symmetries for carbon cycle parameters.
## For the emission driven ESM calirbation using temperature, heatflux, and co2 comparison data.
## Calibrate 6/7 Hector parameters at fixed values of beta.
## Working directory should be set to the analysis/paper_1 directory.
## 0. Set Up -----------------------------------------------------------------------------------------------------
# The required packages
library(hector)
library(dplyr)
library(tidyr)
library(hectorcal)

# Directories
BASE_DIR      <- getwd()                                                           # should be hectorcal/analysis/paper_1
OUTPUT_DIR    <- file.path(BASE_DIR, 'output', 'emiss_beta_q10')                   # where to save the restults to
dir.create(OUTPUT_DIR)
BESTGUESS_DIR <- file.path(BASE_DIR, 'output', 'emiss_temp_heatflux_co2', 'final') # where other emission calibration results are saved, these  fits will inform  the initial guesses for the optim routine.dir.create(OUTPUT_DIR)

# Define the models with the comparison data we would like to fit to,
models_to_fit <- c("GFDL-ESM2G", "CESM1-BGC", "CanESM2")
    #c('NorESM1-ME', 'MRI-ESM1', 'MPI-ESM-LR', 'MIROC-ESM', 'GFDL-ESM2G', 'CESM1-BGC', 'CanESM2')

# The sampled valeus of beta.
beta_values   <- seq(from = 1e-5, to = 1.5, length.out = 10)


## 1. Prep Input Data -----------------------------------------------------------------
# Make the tibble that contains info about the Hector cores we are going to create.
tibble::tibble( file = c(system.file('input/hector_rcp85.ini', package = 'hector'),
                         system.file('input/hector_rcp85.ini', package = 'hector')),
                experiment = c('esmHistorical', 'esmrcp85')) %>%
    select(ini_file = file, core_name = experiment, experiment) ->
    ini_files_tib

# Create the list of the Hector cores.
core_list <- mapply(newcore, inifile = ini_files_tib$ini_file, name = ini_files_tib$core_name )

# Split up the the cmip comparison data by model.
cmip_individual %>%
    filter(model %in% models_to_fit & grepl('esm', experiment)) %>%
    split(.$model, drop = TRUE) ->
    comparison_data


## 2. Calibrate at Fixed Beta ---------------------------------------------------------------------

# For each model tune for 6 Hector paramters at sampled values of beta.
mapply(function(data_name, comp_data){

    # Use the final calibtraion results from other emission driven exercises as optims's best guess.
    object     <- readRDS(list.files(BESTGUESS_DIR, data_name, full.names = TRUE))
    assertthat::assert_that(object$optim_rslt$convergence == 0, msg = 'fit did not converge')
    best_guess <- object$optim_rslt$par[names(object$optim_rslt$par) != BETA()]

    # For each kappa value tune the other 6 Hector paramters.
    lapply(beta_values, function(beta){

        file_name <- paste0('emiss_', data_name, '_', signif(beta, 2))
        file = file.path(OUTPUT_DIR, file_name)
        if(!file.exists(file)){


        # Set the beta value.
        lapply(core_list, parameterize_core, params = c('beta' = beta))
        lapply(core_list, reset)

        # Define the function to minimize in optim.
        fn <- make_minimize_function(hector_cores = core_list,
                                     esm_data = comp_data,
                                     normalize = center_scale,
                                     param = best_guess,
                                     cmip_range = NULL,
                                     n = 2,
                                     showMessages = TRUE,
                                     intermediateOutput = FALSE)

        # Minimize the distance between Hector and ESM output.
        rslt           <- stats::optim(par = best_guess, fn = fn, control = list('maxit' = 3000))
        rslt[['beta']] <- beta

        # Save output
        file_name <- paste0('emiss_', data_name, '_', signif(beta, 2), '.rds')
        saveRDS(rslt, file = file.path(OUTPUT_DIR, file_name))

        }
    })

}, data_name = names(comparison_data), comp_data = comparison_data)


## 3. Format Results  ---------------------------------------------------------------------
# Create a summary table of the calibration results.
list.files(OUTPUT_DIR, '.rds', full.names = TRUE) %>%
    lapply(function(input){

       object <- readRDS(input)
       if(object$convergence == 0){

           info <- unlist(strsplit(gsub(pattern = 'emiss_|.rds', replacement = '', basename(input)), split = '_'))

           tibble(model = info[[1]],
                  beta = as.numeric(info[[2]])) %>%
           cbind(matrix(object$par, nrow = 1, dimnames = list(NULL, names(object$par)))) %>%
               mutate(min_value = object$value)


       }

    }) %>%
    bind_rows() %>%
    filter(beta <= 1) ->
    rslts

# Save output
write.csv(rslts, file = file.path(OUTPUT_DIR, 'summary_table.csv'), row.names = FALSE)


