## Hector Calibration Paper 1
## Concentration driven single ESM calibration at for S, alpha, and volscl at fixed kappa. The point
## of this exercise it to illustrate the symmetry about S and kappa.
## The working directory should be set up to the hectorcal/anaysis/paper_1 direcotry
## 0. Set Up ------------------------------------------------------------------------
# Load the required packages.
library(hector)
library(dplyr)
library(tidyr)
library(hectorcal)


# Set up the script directories, note the base working directory should be set up to
# hectorcal/analysis/paper_1.
BASE_DIR   <- here::here('analysis', 'paper_1')
OUTPUT_DIR <- file.path(BASE_DIR,'output', '2A.symmetry_poster_results')
dir.create(OUTPUT_DIR)

# The models we want to look at, I selected these models to be consistent with the other plots from the
# Hector calibration paper 1.
models_to_fit <- c('CCSM4', 'GFDL-CM3', 'MRI-CGCM3') #c('CMCC-CESM', 'GFDL-CM3', 'MRI-CGCM3')

# The kappa values to fix and the other non kappa values to solve for.
kappa_values <- seq(from = 0.05, to = 8, length.out = 15)
fit_params   <- c(ECS(), AERO_SCALE(), VOLCANIC_SCALE())


## 1. Prep Input Data -----------------------------------------------------------------
# Make the tibble that contains info about the Hector cores we are going to create.
tibble::tibble( file =  c(system.file('input/hector_rcp26_constrained.ini', package = 'hector'),
                          system.file('input/hector_rcp26_constrained.ini', package = 'hector'),
                          system.file('input/hector_rcp45_constrained.ini', package = 'hector'),
                          system.file('input/hector_rcp60_constrained.ini', package = 'hector'),
                          system.file('input/hector_rcp85_constrained.ini', package = 'hector')),
                experiment = c('historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85')) %>%
    select(ini_file = file, core_name = experiment, experiment) ->
    ini_files_tib


# Make a list of the cores to use in make minimize function.
core_list <- mapply(newcore, inifile = ini_files_tib$ini_file, name = ini_files_tib$core_name )



## 2. Temperature-only Fits -----------------------------------------------------------------
# Because the objective of this script is to solve for 3 climate parameters at different values
# of kappa we are not able to use singleESM_calibration or singleESM_calibration_diag functions.
# Instead we are going to have use the make_minimize_function in the optim function.

# Create the temperature only output directory.
TEMP_OUTPUT_DIR <- file.path(OUTPUT_DIR, 'conc_clim_temp')
dir.create(TEMP_OUTPUT_DIR, showWarnings = FALSE)


# Subset the comparison data so that it only contains temperature data for the models
# we are interested in plotting.
cmip_individual %>%
    filter(variable == 'tas' & model %in% models_to_fit) %>%
    split(.$model, drop = TRUE) ->
    temp_comparison_data

# For each of the models tune 3 Hector climate parameters so that Hector reproduces the
# temperature dynamics only using temperature comparison data at fixed values of kappa.
mapply(function(data_name, comp_data){

    # The initial guess for the paramter fits informed by the hectorcal large Hector ensemble.
    best_guess <- generate_inital_guess(comp_data, fit_params)

    # For each kappa value tune the remaining three climate paramters.
    lapply(kappa_values, function(kappa){

        file_name <- paste0('conc_', data_name, '_', kappa, '.rds')
        if(!file.exists(file.path(TEMP_OUTPUT_DIR, file_name))){

            # Set the kappa value.
            names(kappa) <- DIFFUSIVITY()
            lapply(core_list, parameterize_core, params = kappa)
            lapply(core_list, reset)

            # Make the function to minimize
            fn <- make_minimize_function(hector_cores = core_list,
                                         esm_data = comp_data,
                                         normalize = center_scale,
                                         param = best_guess,
                                         cmip_range = NULL,
                                         n = 2,
                                         showMessages = TRUE,
                                         intermediateOutput = FALSE)

            # Use optim to minimize the MSE between Hector and ESM output data
            rslt <- stats::optim(par = best_guess, fn = fn, control = list('maxit' = 800))
            rslt$'diff' <- kappa

            # Save output
            saveRDS(rslt, file = file.path(TEMP_OUTPUT_DIR, file_name))

        }

    })


    if(!file.exists(file.path(TEMP_OUTPUT_DIR, paste0('conc_', data_name, '.rds')))){


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
                                           maxit = 800,
                                           showMessages = TRUE)
        saveRDS(rslt, file = file.path(OUTPUT_DIR, paste0('conc_', data_name, '.rds')))
    }

},
data_name = names(temp_comparison_data),
comp_data = temp_comparison_data)


# Import and format the temperature-only calibration results.
list.files(TEMP_OUTPUT_DIR, pattern = '.rds', full.names = TRUE) %>%
    lapply(function(input){

        object <- readRDS(input)

        kappa_model <- unlist(strsplit(gsub(x = basename(input), pattern = 'conc_|.rds', replacement = ''), split = '_'))

        if(any("convergence" %in% names(object))){

            df <- data.frame(matrix(object$par, nrow = 1, dimnames = list(NULL, names(object$par))))
            df$model <-kappa_model[1]
            df$kappa <- as.numeric(kappa_model[2])
            df$min   <- object$value

            df

        } else {

            data.frame(model = kappa_model[2],
                       min = NA)

        }
    }) %>%
    bind_rows() %>%
    select(model, kappa, S, min) %>%
    arrange(model, kappa) %>%
    na.omit ->
    S_kappa_fits_temp

# Format and Save Output
S_kappa_fits_temp$comp_data <- 'temp'
write.csv(x = S_kappa_fits_temp, file = file.path(TEMP_OUTPUT_DIR, 'summary_fits.csv'), row.names = FALSE)


## 2. Temperature-Heat Flux Fits -----------------------------------------------------------------
# Create the temperature heatflux output directory.
TEMPHF_OUTPUT_DIR <- file.path(OUTPUT_DIR, 'conc_clim_tempHF')
dir.create(TEMPHF_OUTPUT_DIR)

# Subset the comparison data so that it only contains the data for the models
# we are interested in plotting.
cmip_individual %>%
    filter(model %in% models_to_fit) %>%
    mutate(keep = if_else(variable == 'tas', TRUE, FALSE),
           keep = if_else(experiment != 'historical', TRUE, keep)) %>%
    filter(keep) %>%
    select(-keep) %>%
    split(.$model, drop = TRUE) ->
    temp_heatflux_comparison

# For each of the models tune 3 Hector climate parameters so that Hector reproduces the
# temperature and heatflux dynamics of the comparison data.
mapply(function(data_name, comp_data){

    # The initial guess for the parameter fits informed by the hectorcal large Hector ensemble.
    best_guess <- generate_inital_guess(comp_data, fit_params)

    # For each kappa value tune the remaining three climate parameters.
    lapply(kappa_values, function(kappa){

        # Set the kappa value.
        names(kappa) <- DIFFUSIVITY()
        lapply(core_list, parameterize_core, params = kappa)
        lapply(core_list, reset)

        # If the comparison data is missing heatflux data then use the
        # cmip range approach.
        cmip_range <- NULL
        if(!'heatflux' %in% comp_data$variable){

            esm_comparison %>%
                filter(variable == 'heatflux' & experiment %in% comp_data$experiment,
                       experiment != 'historical') %>%
                select(variable, year, lower = mina, upper = maxb, experiment) %>%
                mutate(sig = 0.15) ->
                cmip_range
        }
        # Make the function to minimize
        fn <- make_minimize_function(hector_cores = core_list,
                                     esm_data = comp_data,
                                     normalize = center_scale,
                                     param = best_guess,
                                     cmip_range = cmip_range,
                                     n = 1,
                                     showMessages = TRUE,
                                     intermediateOutput = FALSE)

            # Use optim to minimize the MSE between Hector and ESM output data
            rslt <- stats::optim(par = best_guess, fn = fn, control = list('maxit' = 800))

            table_fn <- make_minimize_function(hector_cores = core_list,
                                               esm_data = comp_data,
                                               normalize = center_scale,
                                               param = best_guess,
                                               cmip_range = cmip_range,
                                               n = 1,
                                               showMessages = TRUE,
                                               intermediateOutput = TRUE)

            rslt$table <- table_fn(rslt$par)

        # Save output
        rslt$diff <- kappa
        file_name <- paste0('conc_', data_name, '_', kappa, '.rds')
        saveRDS(rslt, file = file.path(TEMPHF_OUTPUT_DIR, file_name))

    })


},
data_name = names(temp_heatflux_comparison)[2:3],
comp_data = temp_heatflux_comparison[2:3])

# Import and format the temperature-heat flux calibration results.
list.files(TEMPHF_OUTPUT_DIR, pattern = '.rds', full.names = TRUE) %>%
    lapply(function(input){

        object <- readRDS(input)
print(input)
        kappa_model <- unlist(strsplit(gsub(x = basename(input), pattern = 'conc_|.rds', replacement = ''), split = '_'))

        if(!is.null(object)){

            df <- data.frame(matrix(object$par, nrow = 1, dimnames = list(NULL, names(object$par))))
            df$model     <- kappa_model[1]
            df$kappa     <- as.numeric(kappa_model[2])
            df$min_total <- object$value

            object$table %>%
                dplyr::filter(variable == GLOBAL_TEMP()) %>%
                group_by(experiment) %>%
                summarise(value = mean(value)) %>%
                ungroup %>%
                pull(value) %>%
                mean ->
                temp_MSE

            df$min <- temp_MSE
            df$value <- object$value
            df

        } else {

            data.frame(model = kappa_model[2],
                       min = NA)

        }
    }) %>%
    bind_rows() %>%
    select(model, kappa, S, min, value) %>%
    arrange(model, kappa) %>%
    na.omit ->
    S_kappa_fits_tempHF

# Format and Save Output
S_kappa_fits_tempHF$comp_data <- 'temp heatflux'
write.csv(x = S_kappa_fits_tempHF, file = file.path(TEMPHF_OUTPUT_DIR, 'summary_fits.csv'), row.names = FALSE)

# The End
