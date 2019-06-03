#### Format results of Hector PCA analysis into center and scale values taht can be used in the
#### single ESM calibration functions to noramlize the outputs. The results will be saved as
#### package data.

library(dplyr)
library(tidyr)
library(hectorcal)

# Extend the center and scale values for the concentration driven runs.
# Because not all of models started and ended runs at the same year the center and scale values we
# need to extend the center and scale value from the PCA and treate each experiment and ensemble
# seperaty in the single ESM calibration.

# For the concentration driven runs.
# First extrapolate the center and scale values so that are sufficent entires for all of the
# years of output data.
tibble::tibble(index = names(hectorcal::pc_conc_hflux$scale),
               scale = hectorcal::pc_conc_hflux$scale,
               center = hectorcal::pc_conc_hflux$center) %>%
    tidyr::separate(index, into = c('experiment', 'variable', 'year')) %>%
    dplyr::mutate(year = as.integer(year)) %>%
    dplyr::right_join(hectorcal::cmip_individual %>%
                          dplyr::select(year, variable, experiment) %>%
                          dplyr::filter(!grepl(pattern = 'esm', experiment)) %>%
                          dplyr::distinct(),
                      by = c("experiment", "variable", "year")) %>%
    dplyr::filter(!grepl(pattern = 'esm', experiment) & variable %in% c('tas', 'heatflux')) %>%
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
    filter(!grepl(pattern = 'esm', experiment)) %>%
    select(experiment, ensemble) %>%
    distinct %>%
    mutate(new_experiment = paste0(experiment, '_', ensemble)) ->
    experiment_ensembles

# Use the list of the experiment/ensemble information to expand the center and the scale values
# used to normalize the data.
extrapolated_center_scale ->
 #   full_join(experiment_ensembles,  by = "experiment") %>%
 #   mutate(new_index = paste0(new_experiment, '.', variable, '.', year)) ->
    exp_en_center_scale



# Extend the center and scale values for the emission driven runs.
# Because not all of models started and ended runs at the same year the center and scale values we
# need to extend the center and scale value from the PCA and treate each experiment and ensemble
# seperaty in the single ESM calibration.

# For the concentration driven runs.
# First extrapolate the center and scale values so that are sufficent entires for all of the
# years of output data.
tibble::tibble(index = names(hectorcal::pc_emiss_hflux$scale),
               scale = hectorcal::pc_emiss_hflux$scale,
               center = hectorcal::pc_emiss_hflux$center) %>%
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
    extrapolated_center_scale_emission

# Use the list of the experiment/ensemble information to expand the center and the scale values
# used to normalize the data.
extrapolated_center_scale_emission ->
 #   full_join(experiment_ensembles,  by = "experiment") %>%
 #   mutate(new_index = paste0(new_experiment, '.', variable, '.', year)) ->
    exp_en_center_scale_emissions

# Format results into a list.
center        <- c(exp_en_center_scale$center, exp_en_center_scale_emissions$center)
names(center) <- c(exp_en_center_scale$index, exp_en_center_scale_emissions$index)
scale         <- c(exp_en_center_scale$scale, exp_en_center_scale_emissions$scale)
names(scale)  <- c(exp_en_center_scale$index, exp_en_center_scale_emissions$index)
normalize     <- list('center' = center, 'scale' = scale)

## Save centera and scale values that will be used in the single esm calibration experiment.
center_scale <- normalize
devtools::use_data(center_scale, overwrite=TRUE, compress='xz')

