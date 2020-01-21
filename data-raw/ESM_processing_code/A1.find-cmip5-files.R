# Locate  the CMIP5 netcdf files on pic and creates a to process file that will be used in L01B to calculate the global
# average tas and pco2.
#
# This script is set up to be run on pic with R/3.4.3, see section 0 for details on set up and section 1 for
# decisions about which CMIP5 files to search for. Should be run from the data-raw directory.


# 0. Set Up --------
# Required libs
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)


# Define the directories.
BASE         <- "/pic/projects/GCAM/Dorheim/Dorheim/hectorcal/data-raw"   # The precip_var project location.
CMIP5_DATA   <- "/pic/projects/GCAM/CMIP5-KDorheim"                # Directory containing the cmip5 files to process.
OUTPUT_DIR   <-  BASE                                             # Define the output directory.


# 1. User Decisions ---------
# Select the CMIP5 variables, experiments, ensembles, realms, and time scales
# to process. These vectors will be used to generate a pattern that will select cmip5
# files of interest from the CMIP5_DATA.  If a vector is set to NULL then the pattern
# will match with any variable, experiment, ensemble, domain pattern.

# cmip5 search vectors
VARIABLES       <- c("tas", "co2")                # cmip5 variables to process
EXPERIMENTS     <- c('rcp26', 'rcp45', 'rcp60',
                     'rcp85', 'historical', 'esmHistorical',
                     'esmrcp85')                # cmip5 experiments to process
ENSEMBLES       <- "[a-zA-Z0-9-]{6}"       # search for an ensemble with p1
DOMAINS         <- "Amon"                  # the modeling domain name
MODELS          <- NULL                    # cmip5 models to process


# 2. Find the CMIP5 files to process ------
# converts the cmip5 search vectors from sectioninto strings that will be used to create the regex search pattern.
#
# Args:
#   vecotr: is a vector containing the CMIP models, ensembles, domains, or experiments to serach for. If set to NULL
# the function witll return a pattern that will serach for everything.
# Return: A regex pattern
pattern_gen <- function(vector){

    if(is.null(vector)){
        "[a-zA-Z0-9-]+"
    } else {
        paste(vector, collapse = "|")
    }

}

# Create the CMIP5 serach patterns
varpattern        <- pattern_gen(VARIABLES)
domaninpattern    <- pattern_gen(DOMAINS)
modelpattern      <- pattern_gen(MODELS)
experimentpattern <- pattern_gen(EXPERIMENTS)
ensemblepattern   <- pattern_gen(ENSEMBLES)


# Create the pattern for CMIP5 file names for the cmip 5 files to search for.
cmip_search_pattern <- paste("(", varpattern, ")_(",
                             domaninpattern, ")_(",
                             modelpattern, ")_(",
                             experimentpattern, ")_(",
                             ensemblepattern,
                             ")_([0-9]{6})-([0-9]{6}", # time, set up to search for any 6 digit date
                             ").nc$", sep = "")

# Search for the files
file_list <- list.files(CMIP5_DATA, pattern = cmip_search_pattern, full.names = TRUE, recursive = TRUE)

if(length(file_list) < 1) stop('Could not find any cmip files files matching ', cmip_search_pattern)


# Parse out the cmip5 meta data information from the file name.
tibble(path = file_list) %>%
    mutate(filename = basename(path)) %>%
    separate(filename, into = c("variable", "domain", "model", "experiment",
                                "ensemble", "date"), sep = "_", remove = FALSE) %>%
    mutate(date = gsub(pattern = '.nc', replacement = '', x = date)) %>%
    select(path, variable, domain, model, experiment, ensemble, date) ->
    file_tibble


# 5. Save output --------------

write.csv(file_tibble, file = file.path(OUTPUT_DIR, "CMIP_to_process.csv"), row.names = F)
