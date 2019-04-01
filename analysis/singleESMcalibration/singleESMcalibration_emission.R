## Find the Hector climate system paramters so that Hector emulates individual CMIP5 models.
## Right now this code is set up to emulate the emission driven experiments only to be
## run on pic.
##
## See the setup section for user sepficic changes.

# 0. Set Up -----------------------------------------------------------------------------------

# The directory location of the project on pic -- this will have to be changed by other users.
OUTPUT_DIR <- file.path(getwd(), 'analysis', 'singleESMcalibration', 'rslts')
dir.create(OUTPUT_DIR, recursive = T, showWarnings = F)

# Load hector cal package
devtools::load_all()
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Prep Input Data -----------------------------------------------------------------
# Because not all of models started and ended runs at the same year the center and scale values we
# used to normalize the hector and comparison data from the pc analysis may be missing years for the
# models that ran extra years. Also since we would like to do a joint calibration for the experiment
# and ensembles where we treat each experiment / ensemble combination as comparison data we will
# also need to make sure that the scale and center information contains an entry for each possible
# experiment ensemble combo.

# First extrapolate the center and scale values so that are sufficent entires for all of the
# years of output data.
tibble::tibble(index = names(hectorcal::pc_emiss$scale),
               scale = hectorcal::pc_emiss$scale,
               center = hectorcal::pc_emiss$center) %>%
    tidyr::separate(index, into = c('experiment', 'variable', 'year')) %>%
    dplyr::mutate(year = as.integer(year)) %>%
    dplyr::right_join(hectorcal::cmip_individual %>%
                          dplyr::select(year, variable, experiment) %>%
                          dplyr::filter(grepl(pattern = 'esm', experiment)) %>%
                          dplyr::distinct(),
                      by = c("experiment", "variable", "year")) %>%
    dplyr::filter(grepl(pattern = 'esm', experiment)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(experiment, variable, year) %>%
    split(., interaction(.$experiment, .$variable)) %>%
    lapply(function(input = X){

        input %>%
            dplyr::arrange(experiment, variable, year) %>%
            dplyr::mutate(scale = zoo::na.approx(scale, rule = 2),
                          center =  zoo::na.approx(center, rule = 2))

    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(index = paste0(experiment, '.', variable, '.', year)) ->
    extrapolated_center_scale

# Now duplicate the center and scale values for each experiment / ensemble combination.
#
# Start by determing the combinations of the experimetns and the ensembles. Add a column that will contain
# the new experiment name (a combination of the experiment and ensemble).
hectorcal::cmip_individual %>%
    filter(grepl(pattern = 'esm', experiment)) %>%
    select(experiment, ensemble) %>%
    distinct %>%
    mutate(new_experiment = paste0(experiment, '_', ensemble)) ->
    experiment_ensembles

# Use the list of the experiment/ensemble information to expand the center and the scale values
# used to normalize the data.
extrapolated_center_scale %>%
    full_join(experiment_ensembles,  by = "experiment") %>%
    mutate(new_index = paste0(new_experiment, '.', variable, '.', year)) ->
    exp_en_center_scale

# Format the normalize list.
center        <- exp_en_center_scale$center
names(center) <- exp_en_center_scale$new_index
scale         <- exp_en_center_scale$scale
names(scale)  <- exp_en_center_scale$new_index
normalize     <- list('center' = center, 'scale' = scale)

# Make a mapping file of the hector ini file name and the experiment name
tibble::tibble( file =  c(system.file('input/hector_rcp26.ini', package = 'hector'),
                          system.file('input/hector_rcp85.ini', package = 'hector')),
                experiment = c('esmHistorical', 'esmrcp85')) %>%
    left_join(experiment_ensembles, by = "experiment")  %>%
    select(ini_file = file, core_name = new_experiment, experiment, ensemble) ->
    ini_files_tib

# Default parameters.
param <- c(3.5, 1, 2.7, 0.5, 2.3, 3.5, 285)
names(param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY(), hector::VOLCANIC_SCALE(),
                  hector::BETA(), hector::Q10_RH(), hector::PREINDUSTRIAL_CO2())

# Format the esm data into of the data to use in the calibraiton, each element in the list should
# contain the experiments (experiment - ensemble combinations) that will be used as the comparison data.
hectorcal::cmip_individual %>%
    filter(grepl(pattern = 'esm', experiment)) %>%
    # Save a copy of the experiment information, this is important for the weighted calibration.
    mutate(experiment_copy = experiment,
           experiment = paste0(experiment, '_', ensemble)) %>%
    split(.$model) ->
    esm_data_list

# A. Unweighted Calibration ---------------------------------------------------------------------

system.time(lapply(names(esm_data_list), function(X){

    # Pull out the comparison data.
    comp_data <- esm_data_list[[X]]

    # Ensure that there is comparison data.
    if(nrow(comp_data) > 1){

        # Calibrate to the ESM comparison data and produce the diagnostic outputs.
        rslt <- singleESM_calibration_diag(inifiles = ini_files_tib$ini_file,
                                           hector_names = ini_files_tib$core_name,
                                           esm_data = comp_data,
                                           normalize = normalize,
                                           initial_param = param,
                                           core_weights = NULL,
                                           n_parallel = 6)

        saveRDS(rslt, file = file.path(OUTPUT_DIR, paste0('esm_', X, '.rds')))

    }}))


