## B2.CMIP5_emission_yr_ocean-heat-flux.R
## Find and process the emission heat flux data from the CMIP5 archive that lives on pic
## to calculate the mean ocean heat flux to be used in the CMIP5 hector calbration.
## This script is set up to be run as a batch file on pic, it will take a while to
## run depending on the number of files that are going to be processed.
# 0. Set Up ------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(assertthat)
library(purrr)
library(rlist)
library(ncdf4)

# Define the directories.
TEMP_FILE_BASE_NAME <- "/pic/scratch/dorh012"                         # Define a place to store the interminate netcdfs created during the cdo processing.
CDO_EXE             <- "/share/apps/netcdf/4.3.2/gcc/4.4.7/bin/cdo"   # Where the cdo lives
BASE                <- '/pic/projects/GCAM/Dorheim/hectorcal/data-raw/ESM_processing_code'
OUTPUT              <- file.path(BASE, 'output'); dir.create(OUTPUT)
# Define the CMIP attributes to process.
VARIABLES       <- c("rsds", "rsus", "rlds", "rlus", "hfss", "hfls")       # variable, if null will search for all variables
EXPERIMENTS     <- c("esmHistorical", 'esmrcp85')                   # experiment names, if null will search for all experiments

# 1. Select the CMIP5 Files to Process -----------------------------------------------
# Find the CMIP5 files climate output data from the CMIP5 archive index. In order to
# calculate ocean heat flux we will need to mulitple heat and radiation output files.
# So we will need to make sure we have sufficent data files for each model / experiment /
# ensemble before we want to process the data.
'/pic/projects/GCAM/CMIP5-KDorheim/cmip5_index.csv' %>%
    read.csv(stringsAsFactors = FALSE) %>%
    filter(variable %in% VARIABLES) %>%
    filter(experiment %in% EXPERIMENTS) %>%
    # Make sure that we are only looking at the monthly data.
    filter(domain %in% 'Amon') ->
    flux_files

# Make sure that a model / experiment / ensemble / time netcdf has as
# netcdf for the six required cmip variables.
flux_files %>%
    select(variable, domain, model, experiment, ensemble, time) %>%
    mutate(exists = 1) %>%
    spread(variable, exists) %>%
    na.omit %>%
    select(domain, model, experiment, ensemble)  ->
    sufficent_data

flux_files %>%
    inner_join(sufficent_data, by = c("domain", "model", "experiment", "ensemble")) %>%
    na.omit ->
    cmip_data_files

# But because we only care about ocean heat flux we need to apply a land mask to
# remove the values over the land.
#
# Select the cell area and the land fraction files from the CMIp5 archive.
'/pic/projects/GCAM/CMIP5-KDorheim/cmip5_index.csv' %>%
    read.csv(stringsAsFactors = FALSE) %>%
    filter(variable %in% c("areacella", "sftlf")) %>%
    rename(meta_file_path = file) %>%
    select(meta_file_path, model, variable, experiment, ensemble) %>%
    distinct ->
    meta_file_tibble

# Determine which models have the complete meta data files needed to create the
# land area mask.
meta_file_tibble %>%
    select(model, variable, ensemble) %>%
    mutate(exists = 1) %>%
    distinct() %>%
    spread(variable, exists)  %>%
    na.omit %>%
    select(model, ensemble) ->
    complete_meta_data

# But some of the models have meta data for mulitple experiments which will cause
# problems latter on, find one set of meta data files for each model.
meta_file_tibble %>%
    inner_join(complete_meta_data) %>%
    na.omit %>%
    split(., .$experiment) ->
    meta_files_experiment

historical_meta <- meta_files_experiment$historical

meta_files_experiment$esmHistorical %>%
    filter(!model %in% historical_meta$model) %>%
    bind_rows(historical_meta) ->
    meta_files

# Models to process, because there is variable and meta data.
models_to_process    <- base::intersect(cmip_data_files$model, meta_files$model)
cmip_data_to_process <- filter(cmip_data_files, model %in% models_to_process)
meta_data_to_process <- filter(meta_files, model %in% models_to_process)

write.csv(cmip_data_to_process, file = file.path(OUTPUT, 'cmip5_flux_to-process.csv'), row.names = FALSE)
write.csv(meta_data_to_process, file = file.path(OUTPUT, 'cmip5_meta-data_to-process.csv'), row.names = FALSE)


# 2. Define Processing Functions -----------------------------------------------------
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


    tryCatch({ assert_that(nrow(input_info) == 1, msg = 'Only info for a single model / experiment / ensemble can be processed by function')
        print(input_info)
        # Check to make sure that there is an entry for everysingle required variable.
        required_variables <- c("rsds", "rsus", "rlds", "rlus", "hfss", "hfls")
        assert_that(all(required_variables %in% input$variable), msg = paste('Function needs data for all of the following variables ', required_variables, collapse = ', '))

        # Make sure that the data is from the atmospheric componet.
        assert_that(unique(input$domain) == 'Amon', msg = 'This function can only process monthly atmopheric data')

        input %>%
            distinct() %>%
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
        system2(CDO_EXE, args = c('fldmean',  paste0("-setgridarea,", weight_nc), '-sub', gain_flux_nc, loss_flux_nc, total_flux_nc), stdout = TRUE, stderr = TRUE)

        # Extract results from the netcdf and format into a data frame.
        nc  <- nc_open(total_flux_nc)
        tibble(year = as.integer(substr(x = ncvar_get(nc, 'time'), start = 1, stop = 4)),
               month = as.integer(substr(x = ncvar_get(nc, 'time'), start = 5, stop = 6)),
               value = ncvar_get(nc, 'rsds'),
               variable = 'heatflux',
               ensemble = input_info[['ensemble']],
               model = input_info[['model']],
               experiment = input_info[['experiment']]) %>%
            group_by(year, variable, ensemble, model, experiment) %>%
            summarise(value = mean(value),
                      month_count = n_distinct(month)) %>%
            ungroup %>%
            filter(month_count == 12) %>%
            select(-month_count)->
            output

        # Clean up the intermediate files
        file.remove(cat_rsds_nc, cat_rsus_nc, cat_rlds_nc, cat_rlus_nc,
                    cat_hfss_nc, cat_hfls_nc, total_flux_nc)

        return(output)
    }, error = function(e){input_info %>%
            mutate(problem = TRUE)})


}


# 3. Process the NETCDFS -----------------------------------------------------
# Figure out the cellArea_nc, landFrac_nc, and the OceanAreaWights nc for each model in
# preperation to use the ocean_weights function to actually calculate the ocean area
# weights.
meta_data_to_process %>%
    select(-experiment) %>%
    distinct %>%
    mutate(variable = if_else(variable == 'areacella', 'cellArea_nc', variable)) %>%
    mutate(variable = if_else(variable == 'sftlf', 'landFrac_nc', variable)) %>%
    spread(variable, meta_file_path) %>%
    mutate(out_nc = file.path(TEMP_FILE_BASE_NAME, paste0("OceanAreaWeights_", model, '.nc'))) %>%
    na.omit ->
    for_calculating_ocean_area_weights

# Calculate the ocean area cell weights.
oceanWeight_nc <- unlist(mapply(FUN = ocean_weights,
                                cellArea_nc = for_calculating_ocean_area_weights$cellArea_nc,
                                landFrac_nc = for_calculating_ocean_area_weights$landFrac_nc,
                                out_nc = for_calculating_ocean_area_weights$out_nc))
oceanWeights_tibble <- tibble(model = for_calculating_ocean_area_weights$model,
                              weight_nc = oceanWeight_nc)

cmip_data_files %>%
    left_join(oceanWeights_tibble, by = 'model') %>%
    filter(!is.na(file) & !is.na(weight_nc)) %>%
    rename(variable_file_path = file)  %>%
    filter(model %in% c('bcc-csm1-1', 'MIROC-ESM', 'GFDL-ESM2G', 'NorESM1-ME')) %>%
    split(., list(.$model, .$experiment, .$ensemble), drop = TRUE) %>%
    lapply(calculate_ocean_heat_flux) %>%
    bind_rows ->
    heat_flux_results

# Save the output
write.csv(x = heat_flux_results, file = file.path(BASE, 'CMIP5_emission_yr_ocean-heat-flux-missing.csv'), row.names = FALSE)

# Clean up
file.remove(oceanWeights_tibble$weight_nc)


