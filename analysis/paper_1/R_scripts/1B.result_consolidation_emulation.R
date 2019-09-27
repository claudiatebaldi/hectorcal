## This script concolidates all of the calibration results and uses to the fitted Hector
## parameter values to use Hector as an ESM emulator for the different CMIP models.

# This should be set to the hectorcal/analysis/paper_1 as the working directory
BASE_DIR <- getwd()
OUTPUT_DIR <- file.path(BASE_DIR, 'output')
OUTPUT_DIR_1A <- file.path(OUTPUT_DIR, '1A.calibration_results')

library(dplyr)
library(tidyr)
library(hector)
library(hectorcal)

climate_params <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE())
carbon_params  <- c(PREINDUSTRIAL_CO2(), BETA(), Q10_RH())

# 1. Consolidate the fitted Hector calibraiton results -----------------------------------------------

# Format all of the results for the different calibraiton results stored in the
# 1A.calibration_results directory. This will also import results from the
# intermediate emission driven calibrations, those results should be removed beofore
# saving the consolidated results as a csv output.
list.files(OUTPUT_DIR_1A, pattern = '.rds', full.names = TRUE, recursive = TRUE) %>%
    lapply(function(input){

        object <- readRDS(input)

        if(is.null(object) || object$optim_rslt$convergence != 0){

            message(paste0('problem with: ', input))

        } else {

            # Parse out the part of the file name that will tell us about the calbration
            # protocol that was used. This name will be used to add calibration method
            # and model name to the data frame.
            unique_name <- gsub(x = input, pattern = OUTPUT_DIR_1A, replacement = '')

            # Format the model, method, fitted parmeter values, and optmized value into
            # a data frame.
            data.frame(model = gsub(pattern = '.rds', replacement = '', x = basename(unique_name)),
                       method = gsub(pattern = '/', replacement = '', dirname(unique_name))) %>%
                cbind(t(object$optim_rslt$par)) %>%
                mutate(optmized_value = object$optim_rslt$value,
                       convergence = object$optim_rslt$convergence)

        }

    }) %>%
    bind_rows() ->
    calibration_fits_all

# Remove the intermeidate calibration restults, all of these results will have
# calibrationX where X is some number and will not have fits for all of the
# carbon cycle models.
to_exclude       <- calibration_fits_all$method[grepl(pattern = 'emiss', x = calibration_fits_all$model) & grepl(pattern = 'calibration\\d', x = calibration_fits_all$method)]

# There is also output from several models that we want to exclude because the comparison
# data did not actually meet the protocol requirements.
models_to_exclude <- c('conc_CNRM-CM5-2', 'conc_MRI-ESM1', 'emiss_MRI-ESM')


# Subset and format the final dataset of calibration results.
calibration_fits_all %>%
    filter(!method %in% to_exclude) %>%
    filter(!model %in% models_to_exclude) %>%
    mutate(method = gsub(pattern = 'final', replacement = '', method)) ->
    calibration_fits

# Save the results
write.csv(calibration_fits, file = file.path(OUTPUT_DIR, 'calibration_results.csv'), row.names = FALSE)

# 2. Concentration Driven Hector Runs -------------------------------------------------------------

# Generate all of the Hector concentration driven cores.
mapply(newcore,
       inifile = c(system.file('input/hector_rcp26_constrained.ini', package = 'hector'),
                   system.file('input/hector_rcp45_constrained.ini', package = 'hector'),
                   system.file('input/hector_rcp60_constrained.ini', package = 'hector'),
                   system.file('input/hector_rcp85_constrained.ini', package = 'hector')),
       name = c('rcp26', 'rcp45', 'rcp60', 'rcp85')) ->
    corelist


# Run Hector for the rcps and historical experiments. Save the temperature and heatflux
# data in a single data frame.
calibration_fits %>%
    filter(grepl(pattern = 'conc', x = model)) %>%
    # Remove the carbon cycle parameters because for these runs they should be set to the default values.
    select(-carbon_params) %>%
    split(interaction(.$model, .$method), drop = TRUE) %>%
    lapply(function(input){

        lapply(corelist, function(core){

            # Try catch will prevent paramterizations of Hector that don't work
            # for all of the rcp pathways from throwing errors
            tryCatch({

                # Parameterlize the Hector core
                params <- select(input, climate_params)
                parameterize_core(core = core, params = params)
                reset(core)

                # Hector until 2100
                run(core = core, runtodate = 2100)

                # Format results
                data.frame(model = input[['model']],
                           method = input[['method']], stringsAsFactors = FALSE) %>%
                    cbind(fetchvars(core = core, dates = 1850:2100, vars = c(GLOBAL_TEMP(), HEAT_FLUX()))) %>%
                    rename(experiment = scenario) %>%
                    # Change the variable names to reflect the ESM variable nomencalture
                    mutate(variable = if_else(variable == 'Tgav', 'tas', variable)) %>%
                    # Add the historical experiments
                    mutate(experiment = if_else(year <= 2005, 'historical', experiment)) %>%
                    distinct()

            }, error = function(e){

                data.frame(model = input[['model']],
                           method = input[['method']],
                           experiment = core$name,
                           value = NA,
                           stringsAsFactors = FALSE)

            })


        }) %>%
            bind_rows()

    }) %>%
    bind_rows() ->
    Hector_conc_results



# 3. Emission Driven Hector Runs -------------------------------------------------------------

# Make the emission driven core, since there are only two CMIP experiments that did the emission
# driven runs all of these runs can be done with a single Hector core!
core <- newcore(inifile = system.file('input/hector_rcp85.ini', package = 'hector'), name = 'esmrcp85')

calibration_fits %>%
    filter(grepl(pattern = 'emiss', x = model)) %>%
    split(interaction(.$model, .$method), drop = TRUE) %>%
    lapply(function(input){

        tryCatch({
            # Set up the Hector core with the climate and carbon cycle paramters
            params <- select(input, c(carbon_params, climate_params))
            parameterize_core(params = params, core = core)
            reset(core)

            run(core)

            data.frame(model = input[['model']],
                       method = input[['method']], stringsAsFactors = FALSE) %>%
                cbind(fetchvars(core = core, dates = 1850:2100, vars = c(GLOBAL_TEMP(), HEAT_FLUX(), ATMOSPHERIC_CO2()))) %>%
                rename(experiment = scenario) %>%
                # Change the variable names to reflect the ESM variable nomencalture
                mutate(variable = if_else(variable == 'Tgav', 'tas', variable)) %>%
                # Add the historical experiments
                mutate(experiment = if_else(year <= 2005, 'esmHistorical', experiment)) %>%
                distinct()

        }, error = function(e){

            data.frame(model = input[['model']],
                       method = input[['method']],
                       experiment = core$name,
                       value = NA,
                       stringsAsFactors = FALSE)
        })

    }) %>%
    bind_rows() ->
    Hector_emiss_results


# 4. Format and Save Results -----------------------------------------------------------------
Hector_conc_results %>%
 bind_rows(Hector_emiss_results) %>%
    # Numerical instablity means that Hector ouptput for some of the historical
    # runs may be slightly different from one another by 1e-5, to prevent
    # these duplicates from causing errors in subsequent steps average
    # them together.
    group_by(model, method, experiment, year, variable, units) %>%
    summarise(value = mean(value)) %>%
    ungroup ->
    emulated_hector_output

write.csv(emulated_hector_output, file = file.path(OUTPUT_DIR, 'emulated_Hector_output.csv'), row.names = FALSE)




