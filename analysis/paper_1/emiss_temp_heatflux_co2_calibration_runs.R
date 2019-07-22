## Hector Calibration Paper 1
## Emission driven single ESM calibration using temperature, heatflux, and co2 comparison data.
## Working directory should be set to the analysis/paper_1 directory.
## It turns out that in order to get the optimization routine to solve the parameters have to be introduced
## one at a time and also the number of iterations is going to be fairly high.
## 0. Set Up -----------------------------------------------------------------------------------------------------
# The required packages
library(hectorcal)
library(dplyr)
library(tidyr)
library(hector)

# Set up dirs
BASE_DIR   <- getwd()

## 0.5 Define Functions -------------------------------------------------------------------------------------------

# Check all of the diagnostic single ESM calibration output to see if the fits converged or not
#
# Args
#       output_dir : the directory containing the .rds files with the calibration results.
# Returns: nothing if the test passes but will throw an error if a run did not converge
check_convergence <- function(output_dir){

    # Check to see that all of the runs converged.
    lapply(list.files(output_dir, pattern = '.rds', full.names = TRUE), function(input){
        object <- readRDS(input)

        if(!is.null(object) &  object$optim_rslt$convergence == 0){
            0
        }
    }) %>%
        unlist ->
        converge_codes

    if(sum(converge_codes) != 0){stop('issue with runs converging')}


}




## 1. Prep Input Data -----------------------------------------------------------------
# Make the tibble that contains info about the Hector cores we are going to create.
tibble::tibble( file =  c(system.file('input/hector_rcp26.ini', package = 'hector'),
                          system.file('input/hector_rcp85.ini', package = 'hector')),
                experiment = c('esmHistorical', 'esmrcp85')) %>%
    select(ini_file = file, core_name = experiment, experiment) ->
    ini_files_tib

# Format the esm data into of the data to use in the calibraiton, each element in the list should
# contain the experiments (experiment - ensemble combinations) that will be used as the comparison data.
# Sasve a copy of the cmip experiment name to use to generate the weights.
hectorcal::cmip_individual %>%
    filter(grepl(pattern = 'esm', experiment)) %>%
    # We do not want to have an overlap between the historical and the future data time series. Present the
    # cut the historical experiments off at year 2005.
    mutate(keep = if_else(experiment == 'historical' & year >= 2005, FALSE, TRUE)) %>%
    filter(keep) %>%
    select(-keep) %>%
    filter(year <= 2100) %>%
    select(year, model, ensemble, variable, experiment, value)  %>%
    filter(variable %in% c('tas', 'heatflux', 'co2')) %>%
    mutate(keep = if_else(variable %in% c('tas', 'co2'), TRUE, FALSE)) %>%
    mutate(keep = if_else(variable == 'heatflux' & grepl('rcp', experiment), TRUE, keep)) %>%
    filter(keep) %>%
    select(-keep) ->
    data_subsetExperiment

# Now that we have subset the comparison dataset to include the results from the
# correct experiments we must subset it the data frame so that it only contains
# results for the models that have both temperature and heat flux data!
data_subsetExperiment %>%
    filter(grepl('rcp', experiment)) %>%
    group_by(model) %>%
    summarise(n_var = n_distinct(variable)) %>%
    ungroup() %>%
    filter(n_var == 3) %>%
    pull(model) ->
    models_with_temp_heatFlux

data_subsetExperiment %>%
    filter(model %in% models_with_temp_heatFlux) %>%
    ungroup %>%
    split(.$model, drop = TRUE) ->
    esm_data_list

## 2. Calibration 1 ---------------------------------------------------------------------
# Use the heat flux and tempeature data to fit the four cliamte parameters. Results
# from thi calibration exercise are going to be used to initiate the next
# calibration exercise that uses one of the carbon cycle paramters.
OUTPUT1_DIR  <- file.path(BASE_DIR, 'output', 'emiss_temp_heatflux_co2', 'calibration1')
dir.create(OUTPUT1_DIR, recursive = TRUE, showWarnings = FALSE)

fit_params <- c(ECS(), DIFFUSIVITY(), VOLCANIC_SCALE(), AERO_SCALE())

mapply(function(data_name, comp_data){

        comp_data <- filter(comp_data, variable != 'co2')
            print(head(comp_data))
            # Find the inital guess.
            best_guess <- generate_inital_guess(comp_data, fit_params)

            # Calibrate to the ESM comparison data and produce the diagnostic outputs.
            rslt <- singleESM_calibration_diag(inifiles = ini_files_tib$ini_file,
                                               hector_names = ini_files_tib$core_name,
                                               esm_data = comp_data,
                                               normalize = center_scale,
                                               initial_param = best_guess,
                                               n_parallel = 1,
                                               cmip_range = NULL,
                                               maxit = 2000,
                                               showMessages = TRUE)

            saveRDS(rslt, file = file.path(OUTPUT1_DIR, paste0('conc_', data_name, '.rds')))

}, data_name = names(esm_data_list), comp_data = esm_data_list)

# Check to see if all of the runs converged.
check_convergence(OUTPUT1_DIR)



## 3. Calibration 2 ---------------------------------------------------------------------
# Use the heat flux, tempeature, co2 data to fit the four cliamte parameters + the 1 carbon
# paramters. Results from this calibration exercise are going to be used to initiate the next
# calibration exercise that uses one of the carbon cycle paramters.

OUTPUT2_DIR  <- file.path(BASE_DIR, 'output', 'emiss_temp_heatflux_co2', 'calibration2')
dir.create(OUTPUT2_DIR, recursive = TRUE, showWarnings = FALSE)

mapply(function(data_name, comp_data){

    print(head(comp_data))

    # Find the inital guess from the previous calibration experiment.
    object <- readRDS(list.files(OUTPUT1_DIR, data_name, full.names = TRUE))
    best_guess <- c(object$optim_rslt$par, "C0" = 285)

    # Calibrate to the ESM comparison data and produce the diagnostic outputs.
    rslt <- singleESM_calibration_diag(inifiles = ini_files_tib$ini_file,
                                       hector_names = ini_files_tib$core_name,
                                       esm_data = comp_data,
                                       normalize = center_scale,
                                       initial_param = best_guess,
                                       n_parallel = 1,
                                       cmip_range = NULL,
                                       maxit = 2000,
                                       showMessages = TRUE)

    saveRDS(rslt, file = file.path(OUTPUT2_DIR, paste0('conc_', data_name, '.rds')))

}, data_name = names(esm_data_list), comp_data = esm_data_list)

# Check to see if all of the runs converged.
check_convergence(OUTPUT2_DIR)


## 4. Calibration 3  ---------------------------------------------------------------------
# Use the heat flux, tempeature, co2 data to fit the four cliamte parameters + the 2 carbon
# paramters. Results from this calibration exercise are going to be used to initiate the next
# calibration exercise that uses one of the carbon cycle paramters.

OUTPUT3_DIR  <- file.path(BASE_DIR, 'output', 'emiss_temp_heatflux_co2', 'calibration3')
dir.create(OUTPUT3_DIR, recursive = TRUE, showWarnings = FALSE)

mapply(function(data_name, comp_data){

    print(head(comp_data))

    # Find the inital guess from the previous calibration experiment.
    object     <- readRDS(list.files(OUTPUT2_DIR, data_name, full.names = TRUE))
    best_guess <- c(object$optim_rslt$par, "beta" = 0.25)

    # Calibrate to the ESM comparison data and produce the diagnostic outputs.
    rslt <- singleESM_calibration_diag(inifiles = ini_files_tib$ini_file,
                                       hector_names = ini_files_tib$core_name,
                                       esm_data = comp_data,
                                       normalize = center_scale,
                                       initial_param = best_guess,
                                       n_parallel = 1,
                                       cmip_range = NULL,
                                       maxit = 2000,
                                       showMessages = TRUE)

    saveRDS(rslt, file = file.path(OUTPUT3_DIR, paste0('conc_', data_name, '.rds')))

}, data_name = names(esm_data_list), comp_data = esm_data_list)


check_convergence(OUTPUT3_DIR)


## 5. Calibration 4  ---------------------------------------------------------------------
# Use the heat flux, tempeature, co2 data to fit the four cliamte parameters + the 3 carbon
# paramters. Results from this calibration exercise are going to be used to initiate the next
# calibration exercise that uses one of the carbon cycle paramters.

OUTPUT4_DIR  <- file.path(BASE_DIR, 'output', 'emiss_temp_heatflux_co2', 'calibration4')
dir.create(OUTPUT4_DIR, recursive = TRUE, showWarnings = FALSE)

mapply(function(data_name, comp_data){

    print(head(comp_data))

    # Find the inital guess from the previous calibration experiment.
    object     <- readRDS(list.files(OUTPUT3_DIR, data_name, full.names = TRUE))
    best_guess <- c(object$optim_rslt$par, "q10_rh" = 3.5)

    # Calibrate to the ESM comparison data and produce the diagnostic outputs.
    rslt <- singleESM_calibration_diag(inifiles = ini_files_tib$ini_file,
                                       hector_names = ini_files_tib$core_name,
                                       esm_data = comp_data,
                                       normalize = center_scale,
                                       initial_param = best_guess,
                                       n_parallel = 1,
                                       cmip_range = NULL,
                                       maxit = 2000,
                                       showMessages = TRUE)

    saveRDS(rslt, file = file.path(OUTPUT4_DIR, paste0('conc_', data_name, '.rds')))

}, data_name = names(esm_data_list), comp_data = esm_data_list)


check_convergence(OUTPUT4_DIR)





# Okay and think I my final one should be with the calibration of all the paramter vakyes,











# Check to see that all of the runs converged.
lapply(list.files(OUTPUT2_DIR, full.names = TRUE), function(input){
    object <- readRDS(input)

    if(!is.null(object) &  object$optim_rslt$convergence == 0){
        0
    }
}) %>%
    unlist ->
    converge_codes

if(sum(converge_codes) != 0){stop('issue with runs converging')}



## 6. Calibration 5  ---------------------------------------------------------------------
# Use the heat flux, tempeature, co2 data to fit the four cliamte parameters + the 3 carbon
# paramters. Results from this calibration exercise are going to be used to initiate the next
# calibration exercise that uses one of the carbon cycle paramters.

OUTPUT5_DIR  <- file.path(BASE_DIR, 'output', 'emiss_temp_heatflux_co2', 'final')
dir.create(OUTPUT5_DIR, recursive = TRUE, showWarnings = FALSE)

mapply(function(data_name, comp_data){

    print(head(comp_data))

    # Find the inital guess from the previous calibration experiment.
    object     <- readRDS(list.files(OUTPUT4_DIR, data_name, full.names = TRUE))
    best_guess <- object$optim_rslt$par

    # Calibrate to the ESM comparison data and produce the diagnostic outputs.
    rslt <- singleESM_calibration_diag(inifiles = ini_files_tib$ini_file,
                                       hector_names = ini_files_tib$core_name,
                                       esm_data = comp_data,
                                       normalize = center_scale,
                                       initial_param = best_guess,
                                       n_parallel = 1,
                                       cmip_range = NULL,
                                       maxit = 2500,
                                       showMessages = TRUE)

    saveRDS(rslt, file = file.path(OUTPUT5_DIR, paste0('conc_', data_name, '.rds')))

}, data_name = names(esm_data_list)[2], comp_data = esm_data_list[2])


check_convergence(OUTPUT5_DIR)




