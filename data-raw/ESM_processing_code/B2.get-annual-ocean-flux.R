

## This file should be launched from the hectorcal/data-raw directory.

# 0. Set Up ----------------------------------------------------------------------------------
# Load required libraries
library(dplyr)
library(tidyr)
library(tibble)
library(assertthat)
library(purrr)
library(rlist)
library(ncdf4)

# Define the directories.
CMIP5_META          <- "/pic/projects/GCAM/CMIP5-CHartin"             # Corinne's CMIP5 directory should have the cell area and the land fraction meta data files.
TEMP_FILE_BASE_NAME <- "/pic/scratch/dorh012"                         # Define a place to store the interminate netcdfs created during the cdo processing.
CDO_EXE             <- "/share/apps/netcdf/4.3.2/gcc/4.4.7/bin/cdo"   # Where the cdo lives

# 1. Define Functions -----------------------------------------------------------------------
# Make a netcdf of ocean weights, where cells that are land based have a value of 0 so that when we take the global average those values do not count.
# This function will save the new netcdf as the out_nc.
# Args
#       cellArea_nc: the netcdf containing the cell area
#       landFrac_nc: the netcdf of the land fraction data
#       out_nc: the full name (including the path) of the ocean area weights
# Returns the path to ocean cell area weights netcdf.
ocean_weights <- function(cellArea_nc, landFrac_nc, out_nc){

    # Although CMIP calls the land fraction data files land fraction they are actually percent land on a scale of 0 to 100.
    # ocean area = cell area * (1 - 0.01 * landPercent)
    system2(CDO_EXE, args = c("-mul", cellArea_nc, "-addc,1","-mulc,-0.01", landFrac_nc, out_nc), stdout = TRUE, stderr = TRUE)

    # Return the path to the weights
    out_nc

    }


# Calculate the ocean heat flux.
# Args
#       input: a dataframe containing the netcdfs for a model / experiment and ensembel memeber. It must also contain
#              path to the ocean area weights.
# Returns a dataframe of the heatflux for a single ensemble / model / experiment
calculate_ocean_heat_flux <- function(input){

    # Make sure that there is data for only one model / experiment / ensemble being used
    # in this function.
    input %>%
        select(model, experiment, ensemble) %>%
        distinct ->
        input_info

    assert_that(nrow(input_info) == 1, msg = 'Only info for a single model / experiment / ensemble can be processed by function')

    # Check to make sure that there is an entry for everysingle required variable.
    required_variables <- c("rsds", "rsus", "rlds", "rlus", "hfss", "hfls")
    assert_that(all(required_variables %in% input$variable), msg = paste('Function needs data for all of the following variables ', required_variables, collapse = ', '))

    # Make sure that the data is from the atmospheric componet.
    assert_that(unique(input$domain) == 'Amon', msg = 'This function can only process monthly atmopheric data')

    input %>%
        mutate(variable = paste0(variable, '_nc')) %>%
        spread(variable, variable_file_path) ->
        wide_cmip_data

    # Determine the weight netcdf.
    weight_nc <- unique(wide_cmip_data[['weight_nc']])

    # Create the netcdfs that will be used as the concacetnated netcdfs.
    cat_rsds_nc <- file.path(TEMP_FILE_BASE_NAME, 'cat_rsds.nc')
    cat_rsus_nc <- file.path(TEMP_FILE_BASE_NAME, 'cat_rsus.nc')
    cat_rlds_nc <- file.path(TEMP_FILE_BASE_NAME, 'cat_rlds.nc')
    cat_rlus_nc <- file.path(TEMP_FILE_BASE_NAME, 'cat_rlus.nc')
    cat_hfss_nc <- file.path(TEMP_FILE_BASE_NAME, 'cat_hfss.nc')
    cat_hfls_nc <- file.path(TEMP_FILE_BASE_NAME, 'cat_hfls.nc')

    # Concatenate the netcdfs for a for the same model / variable / ensemble together,
    # this is important for the models that saved the outputs into different time chunks.
    system2(CDO_EXE, args = c("-a", "copy", wide_cmip_data[['rsds_nc']], cat_rsds_nc), stdout = TRUE, stderr = TRUE)
    system2(CDO_EXE, args = c("-a", "copy", wide_cmip_data[['rsus_nc']], cat_rsus_nc), stdout = TRUE, stderr = TRUE)
    system2(CDO_EXE, args = c("-a", "copy", wide_cmip_data[['rlds_nc']], cat_rlds_nc), stdout = TRUE, stderr = TRUE)
    system2(CDO_EXE, args = c("-a", "copy", wide_cmip_data[['rlus_nc']], cat_rlus_nc), stdout = TRUE, stderr = TRUE)
    system2(CDO_EXE, args = c("-a", "copy", wide_cmip_data[['hfss_nc']], cat_hfss_nc), stdout = TRUE, stderr = TRUE)
    system2(CDO_EXE, args = c("-a", "copy", wide_cmip_data[['hfls_nc']], cat_hfls_nc), stdout = TRUE, stderr = TRUE)

    # Now caculate the monthly downward heat flux.
    # rsds-rsus+rlds-rlus-hfss-hfls

    # Calculate total radiation being absorbed by the oceean and the given off by the oceans,
    # the aggregate totals will be used to calculate flux.
    gain_flux_nc <- file.path(TEMP_FILE_BASE_NAME, 'gain_flux.nc')
    loss_flux_nc <- file.path(TEMP_FILE_BASE_NAME, 'loss_flux.nc')

    system2(CDO_EXE, args = c('add', cat_rsds_nc, cat_rlds_nc, gain_flux_nc), stdout = TRUE, stderr = TRUE)
    system2(CDO_EXE, args = c('add', cat_hfls_nc, '-add', cat_rsus_nc, '-add', cat_rlus_nc, cat_hfss_nc, loss_flux_nc), stdout = TRUE, stderr = TRUE)

    # Calculate annual ocean heat flux (gain - loss).
    total_flux_nc <- file.path(TEMP_FILE_BASE_NAME, 'total_flux.nc')
    system2(CDO_EXE, args = c('fldmean',  paste0("-setgridarea,", weight_nc), '-yearmean', '-sub', gain_flux_nc, loss_flux_nc, total_flux_nc), stdout = TRUE, stderr = TRUE)

    # Extract results from the netcdf and format into a data frame.
    nc  <- nc_open(total_flux_nc)
    tibble(year = as.integer(substr(x = ncvar_get(nc, 'time'), start = 1, stop = 4)),
           value = ncvar_get(nc, 'rsds'),
           variable = 'heatflux',
           ensemble = input_info[['ensemble']],
           model = input_info[['model']],
           experiment = input_info[['experiment']]) ->
        output

    # Clean up the intermediate files
    file.remove(cat_rsds_nc, cat_rsus_nc, cat_rlds_nc, cat_rlus_nc,
                cat_hfss_nc, cat_hfls_nc, total_flux_nc)

    return(output)
}


# 2. Import Data --------------------------------------------------------------------------------
# Import the data to process and break it up into the variable data data and the meta data.
to_process <- read.csv('cmip_flux_data_to_process.csv', stringsAsFactors = FALSE)

# Seperate the cmip data files from the cmip meta data files.
to_process %>%
    filter(!is.na(variable_file_path)) %>%
    select(variable_file_path, domain, model, experiment, ensemble, variable) ->
    cmip_data_files

to_process %>%
    filter(is.na(variable_file_path)) %>%
    select(meta_file_path, model, variable) ->
    meta_data_files

# 3. Calculate the Ocean Area Weights -----------------------------------------------------

# Figure out the cellArea_nc, landFrac_nc, and the OceanAreaWights nc for each model in
# preperation to use the ocean_weights function to actually calculate the ocean area
# weights.
meta_data_files %>%
    mutate(variable = if_else(variable == 'areacella', 'cellArea_nc', variable)) %>%
    mutate(variable = if_else(variable == 'sftlf', 'landFrac_nc', variable)) %>%
    spread(variable, meta_file_path) %>%
    mutate(out_nc = file.path(TEMP_FILE_BASE_NAME, paste0("OceanAreaWeights_", model, '.nc'))) ->
    for_calculating_ocean_area_weights

# Calculate the ocean area cell weights.
oceanWeight_nc <- unlist(mapply(FUN = ocean_weights,
                         cellArea_nc = for_calculating_ocean_area_weights$cellArea_nc,
                         landFrac_nc = for_calculating_ocean_area_weights$landFrac_nc,
                         out_nc = for_calculating_ocean_area_weights$out_nc))

# Add the model information.
oceanWeights_tibble <- tibble(model = for_calculating_ocean_area_weights$model,
                              weight_nc = oceanWeight_nc)

# 4. Calculate Ocean Heat FLux -----------------------------------------------------------------
cmip_data_files %>%
    left_join(oceanWeights_tibble, by = 'model') %>%
    split(., list(.$model, .$experiment, .$ensemble)) %>%
    # Remove the entries of the list that contain 0 rows because the
    # caclulate_ocean_heat_flux function will not work with it.
    list.clean(fun = function(input){nrow(input) == 0}) %>%
    lapply(calculate_ocean_heat_flux) %>%
    bind_rows ->
    heat_flux_results

# Save the output
write.csv(x = heat_flux_results, file = 'CMIP5_heat_flux.csv', row.names = FALSE)

# Clean up
file.remove(oceanWeights_tibble$weight_nc)



