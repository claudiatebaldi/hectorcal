#### Format results of hector ensembles into a useable structure.  Some of these
#### will be saved as package data.

library(assertthat)

indir <- 'PCA_pic_results'
outdir <- 'PCA-data'
ensemble_size <- 1000

### Functions that we will use in the processing

## Read the data from a file
read_data <- function(filename) {
    d <- readRDS(filename)
    row.names(d) <- NULL
    d
}

## Get the table of parameters by runid from the data set.
##      d: data frame
## params: names of parameters included in the input data
harvest_params <- function(d, params) {
    rslt <- unique(d[c('runid', params)])
    row.names(rslt) <- NULL             # These things are really annoying!
    rslt
}

## Strip off the parameters from the data set.
##      d: data frame
## params: names of parameters included in the input data
strip_params <- function(d, params) {
    d[!(names(d) %in% params)]
}


## Find all runs that have data in all scenarios.
## dl: list of data frames, one for each scenario
## return: vector of runids that have valid data for all scenarios
valid_runs <- function(dl) {
    ## list of tables of bogus runs
    bl <-
        lapply(dl, function(d) {
                   dplyr::group_by(d, runid) %>%
                     dplyr::summarise(bogus=any(is.na(value)))
               })

    ## check that we have the same runs in all
    for(b in bl[-1]) {
        assert_that(all.equal(b$runid, bl[[1]]$runid))
    }

    ## matrix of true/false values.  Scenarios in columns, run id in rows
    vm <- sapply(bl, function(b) {b$bogus})

    anybogus <- apply(vm, 1, any)

    ## The runids that do not have any bogus are valid
    bl[[1]]$runid[!anybogus]
}


## Filter a set of runs to the desired size, keeping only valid runs
## dl : list of data frames, one for each scenario
## nkeep : number of runs to keep
## valid : vector of valid runids
## return: list of filtered data frames
filter_runs <- function(dl, nkeep, valid)
{
    keep <- valid[1:nkeep]
    lapply(dl, function(d) {d[d$runid %in% keep,]})
}


## Check that we got the right number of runs and no bad data
## d : data fram containing an ensemble of runs for a scenario
## n : number of runs expected
isgood <- function(d, n)
{
    nrow(d) && !any(is.na(d$value))
}

## Add the experiment column to a dataset.  Also fix variable names
##   (Tgav -> tas, atmos_c -> co2)
## d: data frame with runs for an experiment
## regex: Regular expression for extracting the experiment from the scenario
##        name.
## grp: Group number that will contain the experiment name.  This is using
##      stringr numbering, so 1 is the whole match, 2 is the first group, etc.
## return list of data frames with experiment added
add_experiment <- function(d, regex, grp=2)
{
    scen <- unique(d$scenario)
    assert_that(length(scen)==1)
    d$experiment <- stringr::str_match(scen, regex)[1,grp]
    d$variable <- as.character(d$variable)
    d$variable <- dplyr::if_else(d$variable=='Tgav', 'tas',
                                 dplyr::if_else(d$variable=='atmos_c', 'co2', d$variable))
    d
}


### Start with the concentration-driven runs
concpattern <- '^concen-RCP'
concrunfiles <- list.files(indir, concpattern, full.names=TRUE)
tempparams <- c('S', 'alpha', 'diff')

concdata_full <- lapply(concrunfiles, read_data)
params_conc_l <- lapply(concdata_full, harvest_params, params=tempparams)
for(p in params_conc_l[-1]) {
    assert_that(all.equal(p, params_conc_l[[1]]))
}

data_conc_l <- lapply(concdata_full, strip_params, params=tempparams)
concvalid <- valid_runs(data_conc_l)
assert_that(length(concvalid) >= ensemble_size)
data_conc_l <- filter_runs(data_conc_l, ensemble_size, concvalid) %>%
  lapply(add_experiment, regex='concen-(rcp[0-9]+)')
assert_that(all(sapply(data_conc_l, isgood, n=ensemble_size)))

## Write the output.  Use the scenario as the name
for (d in data_conc_l) {
    fn <- file.path(outdir, paste0(unique(d$scenario), '.rds'))
    saveRDS(d, fn, compress='xz')
}
## I debated whether we should filter the parameter list to just the runids that
## we eventually kept, but ultimately I decided to keep them all, in case we
## wanted them for some reason
params_conc <- params_conc_l[[1]][params_conc_l[[1]]$runid %in% concvalid,]
assert_that(!any(is.na(dplyr::left_join(data_conc_l[[1]], params_conc, by='runid')$S)))
saveRDS(params_conc, file.path(outdir, 'params-concen.rds'), compress='xz')


### emissions-driven runs, carbon cycle parameters constant.
emissconstpattern <- '^emissConstantC-RCP'
emissconstfiles <- list.files(indir, emissconstpattern, full.names=TRUE)
emissconstparams <- tempparams

emissconstdata_full <- lapply(emissconstfiles, read_data)
params_emissconst_l <- lapply(emissconstdata_full, harvest_params, params=emissconstparams)
for(p in params_emissconst_l[-1]) {
    assert_that(all.equal(p, params_emissconst_l[[1]]))
}

data_emissconst_l <- lapply(emissconstdata_full, strip_params, params=emissconstparams)
emissconstvalid <- valid_runs(data_emissconst_l)
assert_that(length(emissconstvalid) >= ensemble_size)
data_emissconst_l <- filter_runs(data_emissconst_l, ensemble_size, emissconstvalid) %>%
  lapply(add_experiment, regex='emissConstantC-(esmrcp[0-9]+)')
assert_that(all(sapply(data_emissconst_l, isgood, n=ensemble_size)))

## Write the output.  Use the scenario as the name
for (d in data_emissconst_l) {
    scen <- unique(d$scenario)
    fn <- file.path(outdir, paste0(scen, '.rds'))
    saveRDS(d, fn, compress='xz')
}
params_emissconst <- params_emissconst_l[[1]][params_emissconst_l[[1]]$runid %in% emissconstvalid,]
assert_that(!any(is.na(dplyr::left_join(data_emissconst_l[[1]], params_emissconst, by='runid')$S)))
saveRDS(params_emissconst, file.path(outdir, 'params-emissConstantC.rds'), compress='xz')


### emissions-driven runs, carbon cycle parameters constant.
emisspattern <- '^emissCC-RCP'
emissfiles <- list.files(indir, emisspattern, full.names=TRUE)
emissparams <- c('S', 'alpha', 'diff', 'C0', 'q10_rh', 'beta')

emissdata_full <- lapply(emissfiles, read_data)
params_emiss_l <- lapply(emissdata_full, harvest_params, params=emissparams)
for(p in params_emiss_l[-1]) {
    assert_that(all.equal(p, params_emiss_l[[1]]))
}

data_emiss_l <- lapply(emissdata_full, strip_params, params=emissparams)
emissvalid <- valid_runs(data_emiss_l)
assert_that(length(emissvalid) >= ensemble_size)
data_emiss_l <- filter_runs(data_emiss_l, ensemble_size, emissvalid) %>%
  lapply(add_experiment, regex='emissCC-(esmrcp[0-9]+)')
assert_that(all(sapply(data_emiss_l, isgood, n=ensemble_size)))

## Write the output.  Use the corrected scenario as the name
for (d in data_emiss_l) {
    scen <- unique(d$scenario)
    fn <- file.path(outdir, paste0(scen, '.rds'))
    saveRDS(d, fn, compress='xz')
}
params_emiss <- params_emiss_l[[1]][params_emiss_l[[1]]$runid %in% emissvalid,]
assert_that(!any(is.na(dplyr::left_join(data_emiss_l[[1]], params_emiss, by='runid')$S)))
assert_that(!any(is.na(dplyr::left_join(data_emiss_l[[1]], params_emiss, by='runid')$beta)))
saveRDS(params_emiss, file.path(outdir, 'params-emissCC.rds'), compress='xz')

