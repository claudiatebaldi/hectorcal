#### Load the previously created package data and install it as package data.

library(usethis)

indir <- 'PCA-data'

### Hector ensembles.  We will group these into named lists to cut down on the
### number of data objects.
rcps <- paste0('rcp', c('26','45','60','85'))
conc_expts <- c('historical', rcps)
conc_files <- file.path(indir, paste0('concen-',conc_expts, '.rds'))
hector_conc_ensemble <- lapply(conc_files, readRDS)
conc_params <- list(readRDS(file.path(indir, 'params-concen.rds')))
hector_conc_ensemble <- c(hector_conc_ensemble, conc_params)
names(hector_conc_ensemble) <- c(conc_expts, 'params')
use_data(hector_conc_ensemble, overwrite=TRUE, compress='xz')

emiss_rcps <- paste0('esm', rcps)
emiss_expts <- c('esmHistorical', emiss_rcps)
emiss_files <- file.path(indir, paste0('emissCC-', emiss_expts, '.rds'))
hector_emiss_ensemble <- lapply(emiss_files, readRDS)
emiss_params <- list(readRDS(file.path(indir, 'params-emissCC.rds')))
hector_emiss_ensemble <- c(hector_emiss_ensemble, emiss_params)
names(hector_emiss_ensemble) <- c(emiss_expts, 'params')
use_data(hector_emiss_ensemble, overwrite=TRUE, compress='xz')

### I don't think we are actually going to use the constant carbon-cycle
### ensembles, and having lots of package data really slows down the package
### load.  The input files are in the input directory under data-raw if we ever
### change our minds.
## emiss_constC_files <- file.path(indir, paste0('emissConstantC-esmrcp', rcps,
##                                               '.rds'))
## hector_emiss_constC_ensemble <- lapply(emiss_constC_files, readRDS)
## emiss_constC_params <- list(readRDS(file.path(indir,
##                                               'params-emissConstantC.rds')))
## hector_emiss_constC_ensemble <- c(hector_emiss_constC_ensemble,
##                                   emiss_constC_params)
## names(hector_emiss_constC_ensemble) <- c(paste0('rcp', c('26','45','60','85')),
##                                          'params')
## ## use_data(hector_emiss_constC_ensemble, overwrite=TRUE, compress='xz')


### Principal components.  The only ones I'm going to add as package data are
### the ones we plan to use in the calibration experiment.  The others will be
### available under data-raw to use in side-experiments, and if we get any
### important results from them, we can add them as package data later.
pc_conc <- readRDS(file.path(indir, 'pc-conc-historical-rcp26-rcp45-rcp60-rcp85.rds'))
use_data(pc_conc, overwrite=TRUE, compress='xz')

pc_emiss <- readRDS(file.path(indir, 'pc-emissCC-esmHistorical-esmrcp85.rds'))
use_data(pc_emiss, overwrite=TRUE, compress='xz')

