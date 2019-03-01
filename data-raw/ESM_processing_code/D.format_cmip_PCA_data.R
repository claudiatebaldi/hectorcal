# Convert the cmip strucutre to the PCA ensemble strucutre

# TODO Move out of the raw-data dir? Process the esm comparison full range data?
# What about the individual cmip data?

# 0. Set Up -----------------------------------------------------------------------
library(dplyr)
library(hectorcal)

# 1. Process the Emission driven runs ---------------------------------------------
esm_comparison %>%
    filter(grepl('esm', experiment)) ->
    emiss_runs

# Prevent duplicate years from being saved in the data.
emiss_runs %>%
    filter(experiment == 'esmHistorical') %>%
    # Rename the experiment to match the Hector experiment
    mutate(experiment = 'esmrcp85') ->
    emiss_historical

emiss_runs %>%
    filter(experiment != 'esmHistorical') ->
    future_emiss

future_emiss %>%
    bind_rows(filter(emiss_historical, !year %in% future_emiss$year)) %>%
    arrange(variable, year) %>%
    select(year, variable, experiment, value = cmean) ->
    `PCA_ESMmean-emiss`

# Save as external package data
devtools::use_data(`PCA_ESMmean-emiss`, overwrite = TRUE)

# 2. Process the concentration driven runs ---------------------------------------
# For the concentration runs we are only interested in temperature data and will
# need to combine the historical experiment with the future rcp experiment.
esm_comparison %>%
    filter(!grepl('esm', experiment) & variable == 'tas') ->
    concentration_runs

# Subset the historical experiment and replace with the rcp experiment names.
concentration_runs %>%
    filter(experiment == 'historical') %>%
    mutate(join = 1) %>%
    left_join(tibble(new_experiment = c('rcp26', 'rcp45', 'rcp60', 'rcp85'),
                     join = 1),
              by = 'join') %>%
    select(year, variable, experiment = new_experiment, value = cmean) ->
    historical_concentration_runs

# Combine the future rcp with the newly names hisortical runs and format for output.
concentration_runs %>%
    filter(experiment != 'historical') %>%
    select(year, variable, experiment, value = cmean) ->
    future_concen

# Make sure there aren't duplicate years being saved in the package data.
future_concen %>%
    bind_rows(historical_concentration_runs %>%
                  filter(!year %in% future_concen$year)) %>%
    arrange(variable, experiment, year) ->
    `PCA_ESMmean-concen`

# Save as external package data
devtools::use_data(`PCA_ESMmean-concen`, overwrite = TRUE)


