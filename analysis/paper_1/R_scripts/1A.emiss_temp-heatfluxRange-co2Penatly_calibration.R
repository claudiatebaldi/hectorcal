## Hector Calibration Paper 1
## Emission driven single ESM calibration using temperature, heatflux, and co2 comparison data.
## This is for the emission driven models that have the incomplete heat flux values, for the
## runs that are missing bot co2 data and heat flux data we are not going to calibrate Hector to
## emualte those models individually because that is not enough information provided from a
## single model.
##
## This one adds the q10 penalty function!

## 0. Set Up -----------------------------------------------------------------------------------------------------
# The required packages
library(hectorcal)
library(dplyr)
library(tidyr)
library(hector)

# Set up dirs
BASE_DIR   <- getwd()
OUTPUT_DIR <- file.path(BASE_DIR, 'output', '1A.calibration_results', 'emiss_temp-heatfluxRange-co2Penalty')
dir.create(OUTPUT_DIR)

## Define the penalty function
penalty <- make_param_penalty_function(list('q10_rh' = function(x){
    fn <- mktruncnorm(a = 0, b = Inf, mu = 1.75, sig = 0.4)
    prod(fn(x), -1)

}))



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
    # These emission driven models have some issue with the
    filter(model %in% c('NorESM1-ME', 'MRI-ESM1', 'MIROC-ESM', 'MPI-ESM-LR')) %>%
    filter(variable == 'heatflux') %>% group_by(model, experiment) %>% summarise(year = max(year))
    split(.$model) ->
    esm_data_list



# Use the heat flux, tempeature, co2 data to fit the four cliamte parameters + the 2 carbon
# paramters. Results from this calibration exercise are going to be used to initiate the next
# calibration exercise that uses one of the carbon cycle paramters.

CAL1 <- file.path(OUTPUT_DIR, 'final2')
dir.create(CAL1, recursive = TRUE, showWarnings = FALSE)

mapply(function(data_name, comp_data){

    print(head(comp_data))

    # Find the inital guess from the previous calibration experiment.
    object     <- readRDS(list.files(file.path(BASE_DIR, 'output', '1A.calibration_results',
                                               'emiss_temp-heatfluxRange-co2', 'final'), data_name, full.names = TRUE))
    best_guess <- c(object$optim_rslt$par)

    # Generate the range
    esm_comparison %>%
        filter(experiment == 'esmrcp85' & variable == 'heatflux') %>%
        mutate(sig = 0.15) %>%
        mutate(lower = mina, upper = maxb) ->
        heatflux_range


    # Calibrate to the ESM comparison data and produce the diagnostic outputs.
    rslt <- singleESM_calibration_diag(inifiles = ini_files_tib$ini_file,
                                       hector_names = ini_files_tib$core_name,
                                       esm_data = comp_data,
                                       normalize = center_scale,
                                       initial_param = best_guess,
                                       n_parallel = 1,
                                       cmip_range = heatflux_range,
                                       maxit = 2000,
                                       showMessages = TRUE)

    saveRDS(rslt, file = file.path(CAL1, paste0('emiss_', data_name, '.rds')))

}, data_name = names(esm_data_list), comp_data = esm_data_list)



## 5. Calibration 4  ---------------------------------------------------------------------
# Use the heat flux, tempeature, co2 data to fit the four cliamte parameters + the 3 carbon
# paramters. Results from this calibration exercise are going to be used to initiate the next
# calibration exercise that uses one of the carbon cycle paramters.

OUTPUT4_DIR  <- file.path(OUTPUT_DIR, 'calibration4')
dir.create(OUTPUT4_DIR, recursive = TRUE, showWarnings = FALSE)

mapply(function(data_name, comp_data){

    print(head(comp_data))

    # Find the inital guess from the previous calibration experiment.
    object     <- readRDS(list.files(OUTPUT3_DIR, data_name, full.names = TRUE))
    best_guess <- c(object$optim_rslt$par, "q10_rh" = 3.5)


    # Generate the range
    esm_comparison %>%
        filter(experiment == 'esmrcp85' & variable == 'heatflux') %>%
        mutate(sig = 0.15) %>%
        mutate(lower = mina, upper = maxb) ->
        heatflux_range

    # Calibrate to the ESM comparison data and produce the diagnostic outputs.
    rslt <- singleESM_calibration_diag(inifiles = ini_files_tib$ini_file,
                                       hector_names = ini_files_tib$core_name,
                                       esm_data = comp_data,
                                       normalize = center_scale,
                                       initial_param = best_guess,
                                       penalty = penalty,
                                       n_parallel = 1,
                                       cmip_range = heatflux_range,
                                       maxit = 2000,
                                       showMessages = TRUE)

    saveRDS(rslt, file = file.path(OUTPUT4_DIR, paste0('conc_', data_name, '.rds')))

}, data_name = names(esm_data_list), comp_data = esm_data_list)







