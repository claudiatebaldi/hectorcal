## Hector Calibration Paper 1
## Concentration driven single ESM calibration using temperature and indvidual ESM heat flux comparison data.
## Working directory should be set to the analysis/paper_1 directory.
## 0. Set Up -----------------------------------------------------------------------------------------------------
# The required packages
library(hectorcal)
library(dplyr)
library(tidyr)

# Set up dirs
BASE_DIR   <- getwd()
OUTPUT_DIR  <- file.path(BASE_DIR, 'output', '1A.calibration_results', 'conc_temp-heatflux')
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# A vetor of the hector paramters to optmize.
fit_params <- c(hector::ECS(), hector::DIFFUSIVITY(), hector::AERO_SCALE(), hector::VOLCANIC_SCALE())


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

# Format the esm data into of the data to use in the calibraiton, each element in the list should
# contain historical and future temperature results as well as heatflux results over the future
# scenarios.
hectorcal::cmip_individual %>%
  filter(!grepl(pattern = 'esm', experiment) & variable %in% c('tas', 'heatflux')) %>%
  # We do not want to have an overlap between the historical and the future data time series. Present the
  # cut the historical experiments off at year 2005.
  mutate(keep = if_else(experiment == 'historical' & year >= 2005, FALSE, TRUE)) %>%
  filter(keep) %>%
  select(-keep) %>%
  # Keeps all of the temperature data and the future heatflux data.
  mutate(keep = if_else(variable == 'tas', TRUE, FALSE)) %>%
  mutate(keep = if_else(variable == 'heatflux' & experiment != 'historical', TRUE, keep)) %>%
  filter(keep) %>%
  select(-keep) %>%
  # Remove the inmcm4 values becasue we suspect that there is some error with the heatflux values
  # and will opt to use the CMIP5 range to calibrate instead.
  filter(year <= 2100 & model != "inmcm4") ->
  temp_heatflux_subset


# Now indetify the models that have both temperature and heatflux results. If the model
# is missing heatflux data then it will have to be calibrated using the heatflux range.
temp_heatflux_subset %>%
  group_by(model) %>%
  summarise(n_var = n_distinct(variable)) %>%
  filter(n_var == 2) %>%
  pull(model) ->
  fit_these_models

# Subset the temp and heat flux data set so that it only includes output for the models
# that meet the requirement that it contains temp and heatflux comparison data.
temp_heatflux_subset %>%
  filter(model %in% fit_these_models) %>%
  split(.$model, drop = TRUE) ->
  esm_data_list


## 2. Find Best Fits --------------------------------------------------------------------------------------------
# For all of the data sets in the esm_data_list find the parameters that result in the best fits and generate
# the diagnostic figures.
mapply(function(data_name, comp_data){

  if(!file.exists(file.path(OUTPUT_DIR, paste0('conc_', data_name, '.rds')))){


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


}, data_name = names(esm_data_list), comp_data = esm_data_list)


## 3. Summarise Results --------------------------------------------------------------------------------
list.files(OUTPUT_DIR, '.rds', full.names = TRUE) %>%
  lapply(function(input){

    object <- readRDS(input)

    model_name <- gsub('conc_|.rds', '', basename(input))

    if(object$optim_rslt$convergence  == 0) {

      df           <- as.data.frame(matrix(object$optim_rslt$par, nrow = 1, dimnames = list(NULL, names(object$optim_rslt$par))))
      df$min_value <- object$optim_rslt$value
      df$model     <- model_name
      df$converge  <- 0

    } else {

      df <- data.frame(model = model_name,
                       converge = 1)

    }

    df

  }) %>%
  bind_rows() %>%
  select(model, converge, S, diff, alpha, volscl, min_value) ->
  summary_rslts


if(sum(summary_rslts$converge) == 0){

  write.csv(summary_rslts, file = file.path(OUTPUT_DIR, 'fit_summary_table.csv'), row.names = FALSE)


} else {stop('Some of the runs did not converge')}

