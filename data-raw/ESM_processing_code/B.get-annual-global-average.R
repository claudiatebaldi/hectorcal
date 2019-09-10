# Calculate the annual global average for the CMIP5 files from the netcdfs and save as a csv.
#
# This script is set up to run on pic, see section 0 for set up details. Also uses the
# CMIP_to_process.csv created by LO1.find_cmip5_files.R. This script should be run
# from the hectocal/dta-raw directory.


# 0. Set Up -------------------------------------------------------------------------------------------------
# Load the libs
library(dplyr)
library(ncdf4)

# Directories
BASE        <- '/pic/projects/GCAM/Dorheim/Dorheim/hectorcal/data-raw' # Directory on pic.
INPUT       <- BASE                                            # Define the place that contains the to process csv file.
OUTPUT      <- file.path(BASE, 'output'); dir.create(OUTPUT)   # Define a place to store the final output netcdfs created.
TEMP        <- '/people/dorh012/scratch'; dir.create(TEMP, showWarnings = FALSE)
CDO_EXE     <- "/share/apps/netcdf/4.3.2/gcc/4.4.7/bin/cdo"  # Define the path to the CDO executable.

# Import the csv containing the information about the CMIP 5 monthly nc files to process.
to_process <- read.csv(file.path(INPUT, 'CMIP_to_process.csv'), stringsAsFactors = FALSE)


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

# Use CDO commands to calculate the annual global average from gridded CMIP5 data, this uses
# CDO to convert the output data top absolute time and calculate the weighted average. But
# because not all of the models use the same annual calander the annual average will be
# calculated in R.
#
# Args:
#   path: the file dir of the monthly CMIP data to process.
#   CDO_EXE: the directory of the cdo
#   temp_dir: the location to save all of the intermediate netcds to
#   output_dir: the locaiton where to save the final netcdf files to
#   showMessages: an optional argument that will show messgaes if set to TRUE
#   cleanUp: an optional argument that default set to TRUE delets all of the intermediate netcds
global_annaul_avg <- function(df, CDO_EXE = CDO_EXE, temp_dir = TEMP, output_dir = OUTPUT,
                              showMessages = FALSE, cleanUp = TRUE){

    # Make sure that the data being read in is only for a single cmip file.
    stopifnot(nrow(df) == 1)

    # Ensure the dirs exsist
    stopifnot(file.exists(CDO_EXE))
    stopifnot(file.exists(df[['path']]))
    stopifnot(dir.exists(temp_dir))
    stopifnot(dir.exists(output_dir))

    # Make sure that the path is monthly data.
    stopifnot(grepl(pattern = 'mon', x = df[['path']]))

    # Create the intermediate output files
    global_monthly_nc        <- file.path(TEMP, paste0(basename(df[['path']]), 'MonthlyGlobal.nc'))

    # Becasue some of the modeling groups report results in relative time, days since 1750
    # convert the time of the netcdfds from the relative time to aboslute time (example 1999, 2000,
    # and so on) and calculate the annual global average.
    if(showMessages) message('Calculating the global average ', global_monthly_nc)
    cdo_pipe_msg <- system2(CDO_EXE, args = c('-a', 'fldmean', df[['path']], global_monthly_nc), stdout = TRUE, stderr = TRUE)

    # Import the monthly data and add data info.
    nc <- nc_open(global_monthly_nc)


    # Depending on the variable the dimentions of the value is
    # going to be differe.
    if(df[['variable']] == 'co2'){
        out <- ncvar_get(nc = nc, varid = df[['variable']])[1, ]

    } else{
        out <- ncvar_get(nc = nc, varid = df[['variable']])

    }


    tibble(value = out,
           time = ncvar_get(nc, 'time'),
           units = ncatt_get(nc, df[['variable']], 'units')[['value']]) %>%
        mutate(model = df[['model']],
               variable = df[['variable']],
               experiment = df[['experiment']],
               ensemble = df[['ensemble']]) ->
        monthly_global_rslts

    # Now format the monthly data into annual averages.
    monthly_global_rslts %>%
        mutate(year = substr(time, 1, 4),
               month = substr(time, 5, 6)) %>%
        group_by(year, units, model,  variable, experiment, ensemble) %>%
        summarise(value = mean(value),
                  month_count = n_distinct(month)) %>%
        ungroup %>%
        filter(month_count == 12) %>%
        select(year, units, model,  variable, experiment, ensemble, value) ->
        annual_global

    if(cleanUp){ unlink(global_monthly_nc) }

    # Return the annual global average
    annual_global


}


# 2. Calculate annual global average ---------------------------------------------------------------------------------------------------

# Calculate the annual global average.
out <- apply(X = to_process, MARGIN = 1, FUN = global_annaul_avg, CDO_EXE = CDO_EXE, temp_dir = TEMP,
             output_dir = OUTPUT, showMessages = TRUE, cleanUp = TRUE)

bind_rows(out) %>%
    write.csv(file = file.path(OUTPUT, 'CMIP5_annual_global_average.csv'), row.names = FALSE)




