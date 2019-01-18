library('readr')
library('dplyr')
library('ggplot2')

alldata <- read_csv('CMIP5_annual_global_average.csv')
fixdata <- mutate(alldata, value = if_else(variable=='co2' & value < 1, value*1e6, value))

## Any models for which the values here are dramatically different from the values reported
## for that model in the table BK provided are deemed "bogus" and dropped.
bogusco2 <- c('BNU-ESM', 'FIO-ESM')
## Some of the models didn't run all the way to 2100, for some reason, or don't have
## a historical run
bogustime <- c('GFDL-ESM2M', "bcc-csm1-1",   "bcc-csm1-1-m", "IPSL-CM5A-LR", 'MIROC4h', 'NorESM1-M',
               'GFDL-CM2p1')
## Hadley models have squirrelly data.  Drop them for now
bogushadley <- c('HadCM3', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES')

gooddata <- filter(fixdata, !(model %in% c(bogusco2, bogustime, bogushadley)), 
                   experiment %in% c('rcp26', 'rcp45', 'rcp60', 'rcp85', 'esmrcp85', 'historical', 'esmHistorical'))


## calculate 1850 baseline temperature for historical and esmHistorical runs.
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


#### More functions that will come in handy later.

## detrend and standardize a sequence 
detstd <- function(x) {
    indx <- 1:length(x)
    lin <- lm(x~indx)
    xdet <- x - predict(lin)

    xctr <- xdet - mean(xdet)
    xctr / sd(xctr)   # standardized values
}

## This thing we're calling "erf" is missing the sqrt(2) factor inside the pnorm
## call. Our calls *should* be erf((b-a)/(sqrt(2)*sig)), so dropping the sqrt(2)
## here allows us to drop it to from the args too, saving us a lot of sqrt(2)
## factors in our code that will only end up canceling out anyway.  It does,
## however, mean that this "erf" is not technically the real erf.
erf <- function(x) {2*pnorm(x)-1}
mesa <- function(x, a, b, sig) {(erf((b-x)/sig) - erf((a-x)/sig))/(2*(b-a))}

## Consensus plot, illustrating the style we will be using.  Colors are from the solarized palette
ggplot(data=filter(esm_comparison, experiment=='esmrcp85', year<2100), aes(x=year)) + 
    geom_ribbon(aes(ymin=mina, ymax=maxb), color='#268bd2', fill='#268bd2', alpha=0.5) + 
    geom_line(aes(y=cmean), color='#dc322f', size=1.1) + 
    facet_wrap(facets=~variable, scales='free', strip.position='left',
               labeller = as_labeller(c(tas = "Temperature anomaly (C)", 
                                        co2 = "CO2 concentration (ppm)"))) + 
    theme_minimal(base_size=44) + ylab('') + theme(strip.placement='outside')

a <- 3.363861
b <- 7.34919
x <- seq(a-1, b+1, length.out=1000)
y <- mesa(x, a, b, 1)
pltdata <- data.frame(x=x, y=y)
ggplot(data=pltdata, aes(x=x, y=y)) + 
    geom_line(size=1.1, color='#268bd2') + 
    theme_minimal(base_size=44) + xlab('T') + ylab('p(T)')
