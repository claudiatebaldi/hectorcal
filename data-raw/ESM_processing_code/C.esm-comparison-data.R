# Select and format the ESM comparison data that will be used in the Bayesian calibration of
# Hector, we are calibrating to the full range and the multi model mean for serveral experiments.
# Because of our experimental design we decided to exclude a number of models from our analysis.
# See comments for reasons why particular models were removed.

# Load the required R pacakges
library('readr')
library('dplyr')
library('tidyr')


# Import the csv of the CMIP5 annual global average values for tas and co2 and the
# global average for heat flux data. Some of these models will need to be removed
# and some units will need to be converted.
alldata <- bind_rows(read_csv('CMIP5_annual_global_average.csv'), read_csv('CMIP5_heat_flux_final.csv'))

# Convert CO2 units when applicable.
fixdata <- mutate(alldata, value = if_else(variable=='co2' & value < 1, value*1e6, value))

## Any models for which the values here are dramatically different from the values reported
## for that model in the table BK provided are deemed "bogus" and dropped.
bogusco2 <- c('BNU-ESM', 'FIO-ESM')

## Some of the models didn't run all the way to 2100, for some reason, or don't have
## a historical run
bogustime  <- c('GFDL-ESM2M', "bcc-csm1-1", "IPSL-CM5A-LR", 'MIROC4h', 'NorESM1-M',
               'GFDL-CM2p1')

## Hadley models have squirrelly formatting.  Drop them for now
bogushadley <- c('HadGEM2-ES')

## Make sure that all of the models have BOTH emission driven historical and rcp85 results
## for a variable.
fixdata %>%
    filter(!model %in% c(bogusco2, bogustime, bogushadley),
           grepl('esm', experiment)) %>%
    select(model, experiment, variable) %>%
    distinct %>%
    mutate(data = TRUE) %>%
    spread(experiment, data) %>%
    filter_all(any_vars(is.na(.))) %>%
    pull(model) ->
    incomplete_esm_runs

## Remove the bogus heat flux inmcm4 emission data.
fixdata %>%
    mutate(remove = 0,
           remove = if_else(variable == 'heatflux' & model == 'inmcm4', 1, remove)) %>%
    filter(remove == 0) %>%
    select(-remove) ->
fixdata


## NB: some of the models seems to have bogus data in 1850 (over the whole
## ensemble, min(co2) = 7.5 ppm, max(co2) = 991,012 ppm).  Everything seems to
## be in order in 1851, so we just drop 1850.
gooddata1 <- filter(fixdata, !(model %in% c(bogusco2, bogustime, bogushadley)),
                    year > 1850,
                    experiment %in% c('rcp26', 'rcp45', 'rcp60', 'rcp85', 'esmrcp85', 'historical', 'esmHistorical'))
## Remove the models that have the incomplete emission driven runs for the emission driven experiments only.
gooddata <- filter(gooddata1, !(model %in% incomplete_esm_runs & grepl('esm', experiment)))


## calculate 1850s baseline temperature for historical and esmHistorical runs.
baseline1850 <- function(d)
{
    mostdata <-
        filter(d, year < 1860, experiment %in% c('historical', 'esmHistorical'), variable=='tas') %>%
        group_by(model, ensemble, experiment) %>%
        summarise(temp1850=mean(value)) %>%
        tidyr::spread(experiment, temp1850) %>%
        ungroup() %>%
        group_by(model) %>%
        summarise(esm1850=mean(esmHistorical, na.rm=TRUE), conc1850=mean(historical, na.rm=TRUE))

    ## GFDL starts its runs in 1861, presumably because they think James Buchanan was a chump.
    gfdldata <- filter(d, year < 1870, experiment %in% c('historical', 'esmHistorical'), variable=='tas',
                       model %in% c('GFDL-CM3', 'GFDL-ESM2G')) %>%
        group_by(model, ensemble, experiment) %>%
        summarise(temp1850=mean(value)) %>%
        tidyr::spread(experiment, temp1850) %>%
        ungroup() %>%
        group_by(model) %>%
        summarise(esm1850=mean(esmHistorical, na.rm=TRUE), conc1850=mean(historical, na.rm=TRUE))

    bind_rows(mostdata, gfdldata)
}
gooddata <- left_join(gooddata, baseline1850(gooddata), by='model')
greatdata <- mutate(gooddata,
                    value = if_else((experiment=='esmrcp85' | experiment=='esmHistorical') & variable=='tas', value-esm1850,  # esm runs temperature
                                    if_else(variable=='tas', value-conc1850,      # concentration runs temperature
                                            value)))   %>%                                # not temperature
    filter(year <= 2100)

## Make sure that all of the data starts at the same time for each experiment
## Start by finiding the start year for each experiment.
greatdata %>%
    group_by(experiment) %>%
    summarise(start_year = min(year)) ->
    experiment_start

## Subset the great data so that the final cmip5 data only includes years which
## all the model reported datat for.
greatdata %>%
    left_join(experiment_start, by = 'experiment') %>%
    filter(year >= start_year) %>%
    select(-start_year) %>%
    distinct ->
    final_cmip5_data

# Make sure that the cmip
final_cmip5_data %>%
    mutate(value = 1) %>%
    distinct %>%
    select(year, model, ensemble, variable, experiment, value) %>%
    spread(variable, value) %>%
    dplyr::filter(any(!is.na(heatflux), !is.na(tas))) %>%
    select(model, ensemble, experiment) %>%
    distinct %>%
    mutate(keep = 1) ->
    complete_ensembles

## Reduce ESM data for use in Bayesian calibration.  For each year, variable,
## and each experiment (i.e. type of run), produce the min and max values (for
## consensus calibration), 10 and 90 percentiles (an alternate definition of
## consensus calibration), of all the model runs that performed that experiment.
## (The multi model mean and median will be calculated below).
esm_comparison_range <-
    final_cmip5_data %>%
    # We are not intersted in data beyond 2100. Also inmcm4 has heat flux values
    # that fall outside the range of the other models. We think that that data
    # is fine to use in the "individual cmip5 best fit exercise" we do not to
    # include it in the esm range comparison data set.
    filter(year <= 2100 & !model %in% c('bcc-csm1-1-m')) %>%
    ## TODO for now remove the ESM models that have the gap in co2 data and are
    ## missing the last 10 years of the future runs.
    mutate(keep = if_else(grepl(pattern = 'esm', experiment) & !model %in% c('inmcm4', 'NorESM1-ME', 'MRI-ESM1'), TRUE, FALSE)) %>%
    filter(keep) %>%
    select(-keep) %>%
    left_join(complete_ensembles, by = c("model", "ensemble", "experiment")) %>%
    group_by(year, variable, experiment) %>%
    summarise(mina=min(value), maxb=max(value), a10=quantile(value, 0.1, na.rm = TRUE),
              b90=quantile(value, 0.9, na.rm = TRUE)) %>%
    ungroup

## In order to prevent unequal weights of models in the multimodel mean find the model
## ensemble average for each of the model to use to calculate the multimodel mean
## and median.
final_cmip5_data %>%
    left_join(complete_ensembles, by = c("model", "ensemble", "experiment")) %>%
    filter(year <= 2100 & !model %in% c('bcc-csm1-1-m')) %>%
    ## TODO for now remove the ESM models that have the gap in co2 data and are
    ## missing the last 10 years of the future runs.
    mutate(keep = if_else(grepl(pattern = 'esm', experiment) & !model %in% c('inmcm4', 'NorESM1-ME', 'MRI-ESM1'), TRUE, FALSE)) %>%
    filter(keep) %>%
    select(-keep) %>%
    ## Calculate the model ensemble mean
    group_by(year, model, variable, experiment, unit) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    ## Calculate the multi model mean and median
    group_by(year, variable, experiment) %>%
    dplyr::summarise(cmean = mean(value), cmedian = median(value)) %>%
    dplyr::ungroup() ->
    esm_comaprison_mean

## Combine the esm comparion range and the mulit model mean information into one tibble.
esm_comparison <- left_join(esm_comparison_range,  esm_comaprison_mean, by = c('year', 'variable', 'experiment'))
devtools::use_data(esm_comparison, overwrite=TRUE, compress='xz')

## Save the full table of model data for use in producing parameters for
## emulating individual models.
cmip_individual <- rename(greatdata, esmbaseline=esm1850, concbaseline=conc1850) %>%
    select(-unit)
devtools::use_data(cmip_individual, overwrite=TRUE, compress='xz')

