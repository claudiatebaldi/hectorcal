# Select and format the ESM comparison data that will be used in the Bayesian calibration of
# Hector, we are calibrating to the full range and the multi model mean for serveral experiments.
# Because of our experimental design we decided to exclude a number of models from our analysis.
# See comments for reasons why particular models were removed.

# kalyn okay so i think that what is going on here is that there is some sort of issue with the
# cleaning up of the co2 data because we lost the temp data when we were doing it so we need to get
# it back also I( dont )
# okay so it looks like the start of the co2 runsare padded with 0s so that need to be reflected in
# the data processing but then


# Load the required R pacakges
library('readr')
library('dplyr')
library('tidyr')


# Import the csv of the CMIP5 annual global average values for tas and co2 and the
# global average for heat flux data. Some of these models will need to be removed
# and some units will need to be converted. The concentration driven data.
concdata <- bind_rows(read_csv('CMIP5_annual_global_average.csv'), read_csv('CMIP5_heat_flux_final.csv')) %>%
    filter(!grepl(pattern = 'esm', x = experiment))

## Remove the bogus heat flux inmcm4 emission data.
concdata %>%
    mutate(remove = 0,
           remove = if_else(variable == 'heatflux' & model == 'inmcm4', 1, remove)) %>%
    filter(remove == 0) %>%
    select(-remove) ->
    concdata


## Some of the models didn't run all the way to 2100, for some reason, or don't have
## a historical run
bogustime  <- c('GFDL-ESM2M', "bcc-csm1-1", "IPSL-CM5A-LR", 'MIROC4h', 'NorESM1-M',
                'GFDL-CM2p1')

## Hadley models have squirrelly formatting.  Drop them for now
bogushadley <- c('HadGEM2-ES', 'HadGEM2-CC')

concdata %>%
    filter(!model %in% c(bogustime, bogushadley)) %>%
    filter(variable != 'co2') %>%
    filter(experiment %in% c("historical", "rcp45", "rcp85", "rcp26", "rcp60")) ->
    final_concdata

## Import the emissions data and make sure that the runs meet the requirements.
emsdata <- bind_rows(read.csv('CMIP5_emission_heat_flux_final.csv'),
                     read.csv('CMIP5_emission_yr_global_tas-co2.csv')) %>%
    filter(is.na(problem_file)) %>%
    select(-problem_file) %>%
    filter(!is.na(experiment))

# Convert CO2 units when applicable.
fixdata <- mutate(emsdata, value = if_else(variable=='co2' & value < 1, value*1e6, value))

## Any models for which the values here are dramatically different from the values reported
## for that model in the table BK provided are deemed "bogus" and dropped. These models
## have co2 concentrations that are too high.
bogusco2 <- c('BNU-ESM', 'FIO-ESM')
fixdata  <- filter(fixdata, !model %in% bogusco2)

## Make sure that the co2 data meets the data requiremtns for each variable / experiment.
fixdata %>%
    filter(variable == 'co2' & experiment %in% c('esmHistorical', 'esmrcp85')) %>%
    group_by(ensemble, model) %>%
    summarise(n_exp = n_distinct(experiment)) %>%
    filter(n_exp == 2) %>%
    select(-n_exp) ->
    sufficent_co2

## Make sure that the tas data meets the requirements.
fixdata %>%
    filter(variable == 'tas' & experiment %in% c('esmHistorical', 'esmrcp85')) %>%
    group_by(ensemble, model) %>%
    summarise(n_exp = n_distinct(experiment)) %>%
    filter(n_exp == 2) %>%
    select(-n_exp) ->
    sufficent_tas

fixdata %>%
    filter(variable == 'heatflux' & experiment %in% c('esmrcp85')) %>%
    select(model, ensemble) %>%
    distinct() ->
    sufficent_hf

fixdata %>%
    inner_join(sufficent_co2) %>%
    inner_join(sufficent_tas) %>%
    inner_join(sufficent_hf) ->
    sufficent_emission_data

## Some of the co2 data is quircky, fix these model specific quirks
## Some of the co2 models padded the start of the run with 0s, these values
## must be removed.
sufficent_emission_data %>%
    filter(!(variable == 'co2' & value <= 100)) %>%
    filter(!(variable == 'tas' & value <= 100)) ->
    emissions_without_0padding

# In some instances there are duplicate values for the models that split up
# data for a single experiment into mulitple netcdf files. Here we are going
# to discard those values and replace them inerploated values is an only if
# these duplicated files are for co2 data.
emissions_without_0padding %>%
    group_by(year, variable, ensemble, model, experiment) %>%
    summarise(n = n()) %>%
    filter(n != 1) %>%
    ungroup %>%
    select(year, ensemble, model, experiment, variable) %>%
    mutate(replace = TRUE) ->
    co2_to_interpolate

list <- split(co2_to_interpolate, co2_to_interpolate$model, drop = FALSE)
for(i in seq_along(list)) {

    info <- list[[i]]

    assertthat::assert_that(unique(co2_to_interpolate$variable) == 'co2')
    assertthat::assert_that(nrow(info) == 1)

    emissions_without_0padding %>%
        filter(model == info[['model']] & variable == 'co2' & experiment == info[['experiment']]) %>%
        left_join(info, by = c("year", "variable", "ensemble", "model", "experiment")) %>%
        mutate(replace = if_else(is.na(replace), FALSE, TRUE),
               value = if_else(replace, NA_real_, value)) %>%
        select(-replace) %>%
        distinct ->
        with_NAs

    with_NAs %>%
        arrange(year) %>%
        mutate(value = zoo::na.approx(value)) ->
        good_co2

    emissions_without_0padding %>%
        filter(!c(model == info[['model']] & variable == 'co2' & experiment == info[['experiment']])) %>%
        bind_rows(good_co2) ->
        emissions_without_0padding


}

## MRI-ESM is missing a chunk of time so interplaote the data.
# Is missing a chunk of year so interpolate the data.
years   <- unique(emissions_without_0padding$year)
MRI_co2 <- filter(emissions_without_0padding, variable == 'co2' & model == "MRI-ESM1")
MRI_co2 %>%
    bind_rows(tibble::tibble(year = setdiff(years, MRI_co2$year),
                             value = NA,
                             model = 'MRI-ESM1',
                             variable = 'co2',
                             ensemble = 'r1i1p1',
                             experiment = 'esmrcp85')) %>%
    arrange(year) ->
    MRI_missing

MRI_missing %>%
    arrange(year) %>%
    mutate(value = zoo::na.approx(value)) ->
    good_MRI_co2

## MIROC-ESM is missing a chunk of time so interplaote the data.
# Is missing a single year of data and we need to replace is.
MIROC_co2 <- filter(emissions_without_0padding, model == 'MIROC-ESM' & variable == 'co2')
MIROC_co2 %>%
    bind_rows(tibble::tibble(year = setdiff(years, MIROC_co2$year),
                             value = NA,
                             model = 'MIROC-ESM',
                             variable = 'co2',
                             ensemble = 'r1i1p1',
                             experiment = 'esmHistorical')) %>%
    arrange(year) ->
    MIROC_missing

MIROC_missing %>%
    arrange(year) %>%
    mutate(value = zoo::na.approx(value)) ->
    good_MIROC_co2

emissions_without_0padding %>%
    # Replace the MRI-ESM1 data
    filter(!(variable == 'co2' & model == "MRI-ESM1")) %>%
    bind_rows(good_MRI_co2) %>%
    # Replace the MIROC co2 data
    filter(!(model == 'MIROC-ESM' & variable == 'co2')) %>%
    bind_rows(good_MIROC_co2) %>%
    filter(!(variable == 'heatflux' & experiment == 'esmHistorical'))  ->
    final_emission_data

gooddata <- bind_rows(final_concdata, final_emission_data)


## calculate 1850s baseline temperature for historical and esmHistorical runs.
baseline1850 <- function(d)
{
    mostdata <-
        filter(d, year < 1860, experiment %in% c('historical', 'esmHistorical'), variable=='tas') %>%
        group_by(model, ensemble, experiment) %>%
        summarise(temp1850=mean(value, na.rm = TRUE)) %>%
        tidyr::spread(experiment, temp1850) %>%
        ungroup() %>%
        group_by(model) %>%
        summarise(esm1850=mean(esmHistorical, na.rm=TRUE), conc1850=mean(historical, na.rm=TRUE))

    ## GFDL starts its runs in 1861
    gfdl_models <- unique(d$model[grepl(pattern = 'gfdl', x = tolower(d$model))])
    gfdldata <- filter(d, year < 1870, experiment %in% c('historical', 'esmHistorical'), variable=='tas',
                       model %in% gfdl_models) %>%
        group_by(model, ensemble, experiment) %>%
        summarise(temp1850=mean(value)) %>%
        tidyr::spread(experiment, temp1850) %>%
        ungroup() %>%
        group_by(model) %>%
        summarise(esm1850=mean(esmHistorical, na.rm=TRUE), conc1850=mean(historical, na.rm=TRUE))

    bind_rows(mostdata, gfdldata)
}


gooddata1 <- left_join(gooddata, baseline1850(gooddata), by='model')
greatdata <- mutate(gooddata1,
                    value = if_else((experiment=='esmrcp85' | experiment=='esmHistorical') & variable=='tas', value-esm1850,  # esm runs temperature
                                    if_else(variable=='tas', value-conc1850,      # concentration runs temperature
                                            value)))   %>%                                # not temperature
    filter(year <= 2100)


## Make sure that all of the data starts at the same time for each experiment
## Start by finiding the start year for each experiment.
greatdata %>%
    # Not all of the models starte in the same year.
    group_by(experiment, model) %>%
    summarise(start_year = min(year)) %>%
    ungroup %>%
    group_by(experiment) %>%
    summarise(start_year = max(start_year)) %>%
    ungroup %>%
    mutate(start_year = if_else(grepl('historical', tolower(experiment)), start_year, 2005)) ->
    experiment_start

## Subset the great data so that the final cmip5 data only includes years which
## all the model reported datat for.
greatdata %>%
    full_join(experiment_start, by = 'experiment') %>%
    filter(year >= start_year) %>%
    select(-start_year) %>%
    distinct ->
    final_cmip5_data

## Check the quailty of the data that will be used to calculate the ensemble range information.
assertthat::assert_that({

    ## Final data quality checks before calculating the cmip range.
    ## Make sure that there are not duplicates of data.
    final_cmip5_data %>%
        group_by(year, model, variable, experiment, ensemble) %>%
        summarise(n = n()) %>%
        ungroup %>%
        filter(n != 1) %>%
        select(model, experiment, variable) %>%
        distinct() ->
        rep_years

    rep_years_check <- nrow(rep_years) == 0

    ## Make sure that the data is continuous for all of the years.
    split(final_cmip5_data,
          interaction(final_cmip5_data$experiment,
                      final_cmip5_data$ensemble,
                      final_cmip5_data$variable,
                      final_cmip5_data$model), drop = TRUE) %>%

        lapply(function(input){

            input %>%
                arrange(year) %>%
                pull(year) %>%
                diff %>%
                unique ->
                gaps

            tibble::tibble(gap = gaps) %>%
                mutate(model = unique(input[['model']]),
                       variable = unique(input[['variable']]),
                       ensemble = unique(input[['ensemble']]),
                       experiment = unique(input[['experiment']]))

        }) %>%
        bind_rows() ->
        gaps

        gap_years_check <- nrow(filter(gaps, gap != 1)) == 0

        all(rep_years_check, gap_years_check)
    })


## Reduce ESM data for use in Bayesian calibration.  For each year, variable,
## and each experiment (i.e. type of run), produce the min and max values (for
## consensus calibration), 10 and 90 percentiles (an alternate definition of
## consensus calibration), of all the model runs that performed that experiment.
## (The multi model mean and median will be calculated below).
esm_comparison_range <-
    final_cmip5_data %>%
    group_by(year, variable, experiment) %>%
    summarise(mina=min(value, na.rm = TRUE), maxb=max(value, na.rm = TRUE), a10=quantile(value, 0.1, na.rm = TRUE),
              b90=quantile(value, 0.9, na.rm = TRUE)) %>%
    ungroup

## In order to prevent unequal weights of models in the multimodel mean find the model
## ensemble average for each of the model to use to calculate the multimodel mean
## and median.
final_cmip5_data %>%
    filter(year <= 2100) %>%
    ## Calculate the model ensemble mean
    group_by(year, model, variable, experiment, unit) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    ## Calculate the multi model mean and median
    group_by(year, variable, experiment) %>%
    dplyr::summarise(cmean = mean(value, na.rm = TRUE), cmedian = median(value, na.rm = TRUE)) %>%
    dplyr::ungroup() ->
    esm_comaprison_mean

## Combine the esm comparion range and the mulit model mean information into one tibble.
esm_comparison <- left_join(esm_comparison_range,  esm_comaprison_mean, by = c('year', 'variable', 'experiment'))

## Save the full table of model data for use in producing parameters for
## emulating individual models.
cmip_individual <- rename(final_cmip5_data, esmbaseline=esm1850, concbaseline=conc1850) %>%
    select(-unit)

## The CESM1-BGC model has a weirdly low value in 2031 for the emissions-driven run.  Replace
## it with a spline interpolation
do_interp <- function(master_table, mod, ens, expt, var, interp_years)
{
    yrstrt <- min(interp_years)-5
    yrend <- max(interp_years)+5
    working_data <- dplyr::filter(master_table, model==mod, ensemble==ens, experiment==expt,
                                  variable==var, year >= yrstrt, year <= yrend, !(year %in% interp_years))
    x <- working_data$year
    y <- working_data$value
    ys <- spline(x, y, xout=interp_years)
    locs <- which(master_table$model==mod & master_table$ensemble==ens &
                      master_table$experiment==expt & master_table$variable==var &
                      master_table$year %in% interp_years)
    ## These locations probably correspond to the interpolation years in order, but we can't guarantee that,
    ## so match up the years from the table to the interpolation years
    tbl_yrs <- dplyr::pull(master_table[locs,], 'year')
    locperm <- match(interp_years, tbl_yrs)
    for(i in seq_along(interp_years)) {
        il <- locs[locperm[i]]
        message('Year = ', interp_years[i], ' loc = ', il, ' Replacement= ', ys$y[i],  ' Data=\n',
                paste(format(master_table[il,]), collapse='\n'))
        master_table[il, 'value'] <- ys$y[i]
        message('Replaced data:\n', paste(format(master_table[il, ]), collapse='\n'))
    }
    master_table
}
cmip_individual <- do_interp(cmip_individual, 'CESM1-BGC', 'r1i1p1', 'esmrcp85', 'co2', 2031)

## Check the quailty of the data that is going to be cmip individual comparison data.
assertthat::assert_that({

    ## Final data quality checks before calculating the cmip range.
    ## Make sure that there are not duplicates of data.
    cmip_individual %>%
        group_by(year, model, variable, experiment, ensemble) %>%
        summarise(n = n()) %>%
        ungroup %>%
        filter(n != 1) %>%
        select(model, experiment, variable) %>%
        distinct() ->
        rep_years

    rep_years_check <- nrow(rep_years) == 0

    ## Make sure that the data is continuous for all of the years.
    split(cmip_individual,
          interaction(cmip_individual$experiment,
                      cmip_individual$ensemble,
                      cmip_individual$variable,
                      cmip_individual$model), drop = TRUE) %>%

        lapply(function(input){

            input %>%
                arrange(year) %>%
                pull(year) %>%
                diff %>%
                unique ->
                gaps

            tibble::tibble(gap = gaps) %>%
                mutate(model = unique(input[['model']]),
                       variable = unique(input[['variable']]),
                       ensemble = unique(input[['ensemble']]),
                       experiment = unique(input[['experiment']]))

        }) %>%
        bind_rows() ->
        gaps

    gap_years_check <- nrow(filter(gaps, gap != 1)) == 0

    all(rep_years_check, gap_years_check)
})


## Write the results into the package data.
usethis::use_data(esm_comparison, overwrite=TRUE, compress='xz')
usethis::use_data(cmip_individual, overwrite=TRUE, compress='xz')
