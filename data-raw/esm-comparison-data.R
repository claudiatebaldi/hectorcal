library('readr')
library('dplyr')

alldata <- read_csv('CMIP5_annual_global_average.csv')
## Some models need their units corrected.
fixdata <- mutate(alldata, value = if_else(variable=='co2' & value < 1, value*1e6, value))

## Any models for which the values here are dramatically different from the values reported
## for that model in the table BK provided are deemed "bogus" and dropped.
bogusco2 <- c('BNU-ESM', 'FIO-ESM')
## Some of the models didn't run all the way to 2100, for some reason, or don't have
## a historical run
bogustime <- c('GFDL-ESM2M', "bcc-csm1-1",   "bcc-csm1-1-m", "IPSL-CM5A-LR", 'MIROC4h', 'NorESM1-M',
               'GFDL-CM2p1')
## Hadley models have squirrelly formatting.  Drop them for now
bogushadley <- c('HadCM3', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES')

## NB: some of the models seems to have bogus data in 1850 (over the whole
## ensemble, min(co2) = 7.5 ppm, max(co2) = 991,012 ppm).  Everything seems to
## be in order in 1851, so we just drop 1850.
gooddata <- filter(fixdata, !(model %in% c(bogusco2, bogustime, bogushadley)),
                   year > 1850,
                   experiment %in% c('rcp26', 'rcp45', 'rcp60', 'rcp85', 'esmrcp85', 'historical', 'esmHistorical'))


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
                            value)))                                   # not temperature


## Reduce ESM data for use in Bayesian calibration.  For each year, variable,
## and each experiment (i.e. type of run), produce the min and max values (for
## consensus calibration), 10 and 90 percentiles (an alternate definition of
## consensus calibration), mean (for traditional calibration), and median
## (because, why not?) of all the model runs that performed that experiment.
esm_comparison <-
    group_by(greatdata, year, variable, experiment) %>%
    summarise(mina=min(value), maxb=max(value), a10=quantile(value, 0.1), b90=quantile(value, 0.9),
              cmean=mean(value), cmedian=median(value)) %>%
    ungroup


devtools::use_data(esm_comparison, overwrite=TRUE, compress='xz')

## Save the full table of model data for use in producing parameters for
## emulating individual models.
cmip_individual <- rename(greatdata, esmbaseline=esm1850, concbaseline=conc1850) %>%
    select(-unit)
devtools::use_data(cmip_individual, overwrite=TRUE, compress='xz')
