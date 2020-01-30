## Process annual and the monhtly ocean emission driven heat flux data for the
# 0. Set Up ----------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

# 1. Import Data ------------------------------------------------------------------------------------
heat_flux <- read.csv("CMIP5_emission_yr_ocean-heat-flux.csv", stringsAsFactors = FALSE)

# Import the monthly global fluxes, this data will be used to replace the outliers in the annual heat flux
# time series.
monthly_fluxes <- read.csv("CMIP5_emission_global_monthly_flux.csv")  %>%
    # Parse out the year information from the year_month column.
    rename(year_month = year) %>%
    mutate(year = as.integer(substr(year_month, 1, 4)))


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


# Upon final inspection it looks like there is an issue with the CanESM2 r3i1p1 heat flux output,
# remove that ensemble member from the data set.
final_heat_flux %>%
    filter(!(model == 'CanESM2' & ensemble == 'r3i1p1')) ->
    final_heat_flux

write.csv(x = final_heat_flux, file = 'CMIP5_emission_heat_flux_final.csv', row.names = FALSE)
