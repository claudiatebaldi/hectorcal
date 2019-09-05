## This file adds extra emissions driven results we downloaded and processed latter to the
## already processed CMIP5 data.

# 0. Set Up ----------------------------------------------------------------
OLD_DATA <- file.path(getwd(), 'data-raw')
NEW_DATA <- file.path(getwd(), 'data-raw', 'extra_emiss_data')

library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Import Data --------------------------------------------------------
# Start by importing the original data
global_data_file <- list.files(OLD_DATA, 'CMIP5_annual_global_average.csv', full.names = TRUE)
heatflux_file    <- list.files(OLD_DATA, 'CMIP5_heat_flux_raw.csv', full.names = TRUE)
global_data_og   <- read.csv(global_data_file, stringsAsFactors = FALSE)
heat_flux_og     <- read.csv(heatflux_file, stringsAsFactors = FALSE)


# Now import the extra emissions data that was processed.
global_data_extra <- read.csv(list.files(NEW_DATA, 'CMIP5_annual_global_average.csv',  full.names = TRUE),
                              stringsAsFactors = FALSE) %>%
    mutate(unit = as.character(unit))
heat_flux_extra   <- read.csv(list.files(NEW_DATA, 'CMIP5_heat_flux.csv', full.names = TRUE),
                              stringsAsFactors = FALSE)


# 2. Combine Data Sets --------------------------------------------------------
global_data_og %>%
    bind_rows(global_data_extra) %>%
    distinct ->
    global_data

heat_flux_og %>%
    bind_rows(heat_flux_extra) %>%
    distinct() ->
    global_heat_flux

# 3. Save the Output --------------------------------------------------------
write.csv(global_data, file = global_data_file, row.names = FALSE)
write.csv(global_heat_flux, file = heatflux_file, row.names = FALSE)
