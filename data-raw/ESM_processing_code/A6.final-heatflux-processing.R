## When we made the intial plots of the CMIP5 ocean heat flux that we cauclated we noticed that there
## were some outlier points (values with heat flux up to a magnitude of 400 W/m2). It is unclear what
## is causing this to happen but it looks like some of the netcdf files when processed with the cdo
## code are return multiple values for the same month (this is happening inconsistently across models /
## ensembles / variables / years). If we discard the entires with mulitple monthly values and then
## calcualte the heat flux with these values we can use these values to replace the unrealistic ones.
##
## What is causing this problem is not clear and we are not going to spend time trying to resolve it
## right now but this method replaces the unrealistic values with reasonable and presumably accurate
## values.
##
## This script should be run with the working directory set to hectorcal/data-raw it requires that
## B2 and B3 have already been run.

# 0. Set Up ----------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)


# 1. Import Data ------------------------------------------------------------------------------------
# Import the heat flux data, because of some issues with the data processing there are some wonky heat
# flux values (heat flux values that have a magnitude greater than 5).
heat_flux <- read.csv("CMIP5_heat_flux_raw.csv", stringsAsFactors = FALSE)

# Import the monthly global fluxes, this data will be used to replace the outliers in the annual heat flux
# time series.
monthly_fluxes <- read_csv("CMIP5_global_flux_means.csv_og.zip")  %>%
    # Parse out the year information from the year_month column.
    rename(year_month = year) %>%
    mutate(year = as.integer(substr(year_month, 1, 4))) %>%
    select(-`X1`)

# 2. Calcualte Annual Heat Flux Monthly Global Means --------------------------------------------------

monthly_fluxes %>%
    # We determined that the months with extra values have entires that are set equal to 0.000, remove
    # these observations before calculating the annual average.
    filter(abs(value) > 1e-3) %>%
    group_by(year, variable, ensemble, model, experiment) %>%
    summarise(value = mean(value, na.rm = TRUE),
              n = n()) %>%
    filter(n >= 11) %>%
    ungroup %>%
    spread(variable, value) %>%
    # Some of the variables start at different years, discard the incomplete
    # years here.
    arrange(model, ensemble, experiment, year) %>%
    na.omit %>%
    # Calcualte the new heat flux value.
    mutate(new_value = rsds-rsus+rlds-rlus-hfss-hfls) ->
    new_heat_flux_fullTibble

# Make a dataframe of the new heat flux value that only contains the information that will be used in
# the join with the original heat flux data frame.
new_heat_flux <- select(new_heat_flux_fullTibble, year, ensemble, model, experiment, new_value)

# 3. Replace Outlier Heat Flux Values -----------------------------------------------------------------

# Determine if the heat flux value is an outlier or not. Here we define outliers as values with a magnitude
# greater than 4.5 for the non rcp85 scenarios and greater than 8 for the rcp85 sceanrio.
heat_flux %>%
    mutate(replace = if_else(abs(value) > 4.5 & experiment != 'rcp85', TRUE, FALSE)) %>%
    mutate(replace = if_else(abs(value) > 8 & experiment == 'rcp85', TRUE, replace)) %>%
    distinct ->
    heatflux_to_replace

heatflux_to_replace %>%
    filter(replace) %>%
    inner_join(new_heat_flux, by = c("year", "ensemble", "model", "experiment")) %>%
    mutate(value = if_else(replace, new_value, value)) %>%
    select(names(heat_flux)) ->
    heatflux_bad_values_replaced

heatflux_to_replace %>%
    filter(!replace) %>%
    bind_rows(heatflux_bad_values_replaced) %>%
    select(names(heat_flux)) %>%
    mutate(unit = 'W/m2') ->
    final_heat_flux


# 4. Save Output ------------------------------------------------------------------------------------
write.csv(x = final_heat_flux, file = 'CMIP5_heat_flux_final.csv', row.names = FALSE)



