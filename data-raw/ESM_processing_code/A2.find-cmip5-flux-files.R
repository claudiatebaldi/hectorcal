## The purpose of this script is to find the location of the CMIP5 downward and upward flux files that
## will be used to caculate the ocean heat flux. This script makes sure that we have all of the
## variable files (for the different types of upwards/downwards radition & heat flux) and the
## cell area and land fraction meta data files that will be used to filter data over land.
##
## This script returns a csv file of the CMIP5 files to prcoess.
## This should be launched on pic from the hectorcal/data-raw directory.

# 0. Set Up -----------------------------------------------------------------------------------------------------
# Load required libraries
library(dplyr)
library(tidyr)
library(tibble)

# Define the directories.
CMIP5_DATA          <- "/pic/projects/GCAM/CMIP5-KDorheim"            # The location of the cmip 5 files to process'
CMIP5_META          <- "/pic/projects/GCAM/CMIP5-KDorheim"             # Corinne's CMIP5 directory should have the cell area and the land fraction meta data files.
BASE                <- "/pic/projects/GCAM/Dorheim/Dorheim/hectorcal/data-raw"

# Select the CMIP5 variables, experiments, ensembles, realms, and time scales to process. Any options defined in
# in each of the following vectors will be processed.

VARIABLES       <- c("rsds", "rsus", "rlds", "rlus", "hfss", "hfls")       # variable, if null will search for all variables
EXPERIMENTS     <- c("esmHistorical", 'esmrcp85', 'historical',
                     'rcp26', 'rcp45', 'rcp60', 'rcp85')                   # experiment names, if null will search for all experiments
ENSEMBLES       <- NULL                                                    # the ensemble name, if null will search for all rensembles
DOMAINS         <- "Amon"                                                  # The modeling domain name
MODELS          <- NULL                                                    # A vector of the models to serach for, if null then search for all

# 1. Define Functions  -----------------------------------------------------------------------
# pattern_gen is a function that converts the cmip search vectors from section
# into strings that will be used to create the regex search pattern. If the
# input vector is set to NULL then the function will return a search all pattern.
pattern_gen <- function(vector){

    if(is.null(vector)){
        "[a-zA-Z0-9-]+"
    } else {
        paste(vector, collapse = "|")
    }

}

# 2. Find the CMIP5 Data Files ----------------------------------------------------------------------
# Generate the cmip5 file name patterns.
varpattern        <- pattern_gen(VARIABLES)
modelpattern      <- pattern_gen(MODELS)
experimentpattern <- pattern_gen(EXPERIMENTS)
domaninpattern    <- pattern_gen(DOMAINS)
ensemblepattern   <- pattern_gen(ENSEMBLES)

# Create the pattern for the isimip file names for the netcdfs to search for.
variable_serach_pattern <- paste("(", varpattern, ")_(", domaninpattern, ")_(", modelpattern,
                                 ")_(", experimentpattern, ")_(", ensemblepattern, ")_([0-9]{6})-([0-9]{6}",   # currently set up for 6 digit dates
                                 ").nc$", sep = "")

# Find all of the files that contains the data of interst.
variable_filelist <- list.files(CMIP5_DATA, pattern = variable_serach_pattern, full.names = TRUE, recursive = TRUE)


# Parse the cmip file information out from the file name.
tibble(variable_file_path = variable_filelist) %>%
    mutate(filename = basename(variable_file_path)) %>%
    separate(filename, into = c("variable", "domain", "model", "experiment",
                                "ensemble", "date"), sep = "_", remove = FALSE) %>%
    mutate(date = gsub('.nc', '', date)) %>%
    select(variable_file_path, variable, domain, model, experiment, ensemble, date) ->
    variable_file_tibble

# Determine which models / experiments / ensemble ares missing data. Each one should have
# data for the following variables("rsds", "rsus", "rlds", "rlus", "hfss", "hfls").
# If any of the variables are missing then cannot caculate the ocean heat flux from the ocean data.
# The missing files will have to be downloaded (is possible).
variable_file_tibble %>%
    select(variable, domain, model, experiment, ensemble, date) %>%
    # Spread the experiment column and fill with the binary this file exits. When we spread the
    # data the model / experiemnts/ ensmeble / variable combo is missing an experment then the new columns
    # named after the variables will contain NA codes.
    mutate(exists = 1) %>%
    spread(variable, exists) ->
    variable_wide

# Which models / experiments / ensembles / varaibles are missing variables? This info may be useful with
# trouble shooting problems.
variable_wide %>%
    filter_all(any_vars(is.na(.))) ->
    missing_variables

# Subset the variable wide so that it only contains information for the models / experiments / ensemble
# members that have the complete variable list and use this to subset the variable_file_tibble to
# make a data frame that contains the data files that we would like to process.
variable_wide %>%
    na.omit %>%
    select(domain, model, experiment, ensemble) %>%
    mutate(keep = 1) %>%
    inner_join(variable_file_tibble, by = c("domain", "model", "experiment", "ensemble")) %>%
    filter(keep ==1) %>%
    # Remove the keep column
    select(-keep) ->
    good_cmip_files

# 3. Find the Meta Data Files -------------------------------------------------------------------------
# Create the pattern for CMIP5 file names related to area, cell area and land fraction to search for.
#
# Search for any experiment and ensemble memeber, what matters most about the area realated ncdfs is
# the model.
META         <- c("areacella", "sftlf")
metapattern  <- pattern_gen(META)

# Since there is no specific experiment or ensmeble pattern associated with the cell area and land fraction
# cells reset the patterns to anything.
experimentpattern <- ensemblepattern <- "[a-zA-Z0-9-]+"

# Use the model patterns to serach for the cell area dn ladn fraction netcdfs.
meta_search_pattern <- paste("(", metapattern, ")_(fx)_(", modelpattern, ")_(", experimentpattern, ")_(", ensemblepattern,").nc$", sep = '')
meta_filelist       <- list.files(CMIP5_META, pattern = meta_search_pattern, full.names = TRUE, recursive = TRUE)


# Parse the cmip file information out from the file name.
tibble(meta_file_path = meta_filelist) %>%
    mutate(filename = basename(meta_filelist)) %>%
    separate(filename, into = c("variable", "domain", "model", "experiment",
                                "ensemble"), sep = "_", remove = FALSE) %>%
    select(meta_file_path, model, variable) %>%
    distinct ->
    meta_file_tibble


# Make sure that there is a cell area and a land fraction netcdf for each model, if not then we will not be
# able to process the model.
meta_file_tibble %>%
    select(model, variable) %>%
    mutate(exists = 1) %>%
    distinct() %>%
    spread(variable, exists) ->
    wide_meta

# Which models have do not have both the cell area and the land fraction data?
wide_meta %>%
    filter_all(any_vars(is.na(.))) ->
    missing_meta

# Subset the wide data frame of meta data so that it only contains information for the models that have
# sufficent meta data. wide so that it only contains information for the models / experiments / ensemble
# members that have the complete variable list and use this to subset the variable_file_tibble to
# make a data frame that contains the data files that we would like to process.
wide_meta %>%
    na.omit %>%
    select(model) %>%
    mutate(keep = 1) %>%
    inner_join(meta_file_tibble, by = c("model")) %>%
    filter(keep ==1) %>%
    # Remove the keep column
    select(-keep) ->
    good_meta_files


# 4. Determine the Files to Process and Save Output ----------------------------------------------------------

# In order to caluclate the ocean heat flux we need to use the meta data to determine the fulx over the
# oceans, so if a model is missing meta data then we cannot porcess it.
models_to_process <- base::intersect(good_cmip_files$model, good_meta_files$model)

# The files that we cannot process becasue it is missing stuff.
good_cmip_files %>%
    bind_rows(good_meta_files) %>%
    filter(!model %in% models_to_process) ->
    cmip_data_not_processed

good_cmip_files %>%
    bind_rows(good_meta_files) %>%
    filter(model %in% models_to_process) ->
    cmip_data_to_process

# Save the cmip data that should be processed in the raw-data intermeidate file.
write.csv(cmip_data_to_process, file = file.path(BASE, 'cmip_flux_data_to_process.csv'), row.names = FALSE)
