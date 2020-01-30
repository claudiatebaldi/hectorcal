## B1.CMIP5_emission_yr_global_tas-co2.R
## Calculate the annual global average for the emission driven CMIP5 models.
##
## This script is set up to run on pic, see section 0 for set up details.

# 0. Set Up ----------
# Load the libs
library(dplyr)
library(ncdf4)

# Directories
BASE        <- '/pic/projects/GCAM/Dorheim/hectorcal/data-raw/ESM_processing_code' # Directory on pic.
OUTPUT      <- file.path(BASE, 'output'); dir.create(OUTPUT)   # Define a place to store the final output netcdfs created.
TEMP        <- '/people/dorh012/scratch'; dir.create(TEMP, showWarnings = FALSE)
CDO_EXE     <- "/share/apps/netcdf/4.3.2/gcc/4.4.7/bin/cdo"  # Define the path to the CDO executable.

# Import the csv containing the information about the CMIP 5 monthly nc files to process.
'/pic/projects/GCAM/CMIP5-KDorheim/cmip5_index.csv' %>%
    read.csv(stringsAsFactors = FALSE) %>%
    filter(variable %in% c('tas', 'co2') & domain == 'Amon') %>%
    filter(experiment %in% c('esmHistorical', 'esmrcp85')) %>%
    rename(path = file)  ->
    to_process


# 1. Define Functions -------------

# Check that an object has the expected names.
#
# Args:
#   object: a list, tibble, or dataframe. The object to check for expected names,
#   req_names:  vector of the expected object names.
#   object_name: an optional string argument that if set to a string will be incorperated an error message
# Returns: this function throws an error if the object is missing an expected name
check_names <- function(object, req_names, object_name = NULL){

    missing <- !req_names %in% names(object)

    if(any(missing)){

        stop(object_name, ' is missing ', paste(req_names[missing], collapse = ', '))

    }
    NULL
}

# Use CDO commands to calculate the annual global average from gridded CMIP5 data, this uses CDO.
#
# Args:
#   df: the data frame of the cmip5 files to process
#   CDO_EXE: the directory of the cdo
#   temp_dir: the location to save all of the intermediate netcds to
#   output_dir: the locaiton where to save the final netcdf files to
#   showMessages: an optional argument that will show messgaes if set to TRUE
#   cleanUp: an optional argument that default set to TRUE delets all of the intermediate netcds
annual_global_avg <- function(df, CDO_EXE = CDO_EXE, temp_dir = TEMP, output_dir = OUTPUT,
                              showMessages = FALSE, cleanUp = TRUE){

    # Check inputs
    stopifnot(file.exists(CDO_EXE))
    stopifnot(dir.exists(temp_dir))
    stopifnot(dir.exists(output_dir))

    check_names(df, c('model', 'variable', 'ensemble', 'experiment', 'path'))

    cmip_info <- dplyr::distinct(df[ , names(df) %in% c('model', 'variable', 'experiment', 'ensemble')])

    tryCatch({


        if(nrow(cmip_info) > 1) stop('Multiple model, variable, experiment, ensemble information read into function')

        if(nrow(cmip_info) == 1){
            # Create the intermediate output files
            file_basename        <- paste(cmip_info, collapse = '_')
            concat_nc            <- file.path(temp_dir, paste0('concate_', file_basename, '.nc'))
            global_annual_avg_nc <- file.path(output_dir, paste0('global_annual_avg-', file_basename, '.nc'))

            # Remove files
            if(file.exists(concat_nc)) file.remove(concat_nc)
            if(file.exists(global_annual_avg_nc)) file.remove(global_annual_avg_nc)

            # Concatenate the netcdf files together
            if(showMessages) message("Concatenating files and converting to absolute time ", concat_nc)
            cdo_cat_msg <- system2(CDO_EXE, args = c("-a", "cat", df[['path']], concat_nc), stdout = TRUE, stderr = TRUE)

            # Calculate the annual global average
            if(showMessages) message('Calculating the global annual average ', global_annual_avg_nc)
            cdo_pipe_msg <- system2(CDO_EXE, args = c('yearmean', '-fldmean', concat_nc, global_annual_avg_nc), stdout = TRUE, stderr = TRUE)

            if(cleanUp){ unlink(concat_nc) }

            # Return the path to the output file
            global_annual_avg_nc
        }


    }, error = function(e){NULL})




}


# 2. Calculate annual global average ---------------------------------------------------------------------------------------------------

# Calculate the annual global average.
input         <- split(to_process, interaction(to_process$model, to_process$variable, to_process$ensemble, to_process$experiment), drop = TRUE)
global_avg_nc <- lapply(input, annual_global_avg, CDO_EXE = CDO_EXE, temp_dir = TEMP, output_dir = OUTPUT, showMessages = TRUE, cleanUp = TRUE)


# Format the output into a single data frame
files <- unlist(global_avg_nc)

# For all of the final output netcdfs extract the data and concatenate together.
lapply(files, function(file){

    # Add a try catch statement to prevent errors from breaking the code.
    tryCatch({
        if(file.exists(file)){

            # Format the CMIP5 meta data information.
            name             <- gsub(x = basename(file), pattern = '.nc', replacement = '')
            cmip_info        <- unlist(strsplit(gsub(x = name, pattern = 'global_annual_avg-|globalAvg_', replacement = ''), '_'))
            names(cmip_info) <- c('variable', 'model', 'experiment', 'ensemble')

            # Get the data from the netcdf file
            nc   <- ncdf4::nc_open(file)
            data <- ncdf4::ncvar_get(nc, cmip_info[['variable']])

            if(cmip_info[['variable']] == 'co2'){
                # Since CO2 returns results at different levels we want to extract co2 from the
                # first elevation
                data <- data[1, ]
            }

            time <- as.integer(substr(x = ncdf4::ncvar_get(nc, 'time'), start = 1, stop = 4))
            unit <- ncdf4::ncatt_get(nc, cmip_info[['variable']])$units


            # Format as a data frame
            data.frame(year = time,
                       model = cmip_info[['model']],
                       ensemble = cmip_info[['ensemble']],
                       variable = cmip_info[['variable']],
                       experiment = cmip_info[['experiment']],
                       value = data,
                       unit = unit)

        } else {

            message('Could not find ', file)
        }


        # If there was some problem with the netcdf file return it as a problem file.
    }, error = function(e){data.frame(problem_file = file)})


}) %>%
    bind_rows %>%
    bind_rows(data.frame(problem_file = NA)) ->
    output


# 3. Save the output ---------------
write.csv(output, file = file.path(BASE, 'CMIP5_emission_yr_global_tas-co2.csv'), row.names = FALSE)


file.remove(files)



