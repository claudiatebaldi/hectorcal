#### Compute principal components from ensembles of hector runs

library(assertthat)
library(foreach)
library(hectorcal)
library(dplyr)

### The first thing we do is, we figure out what years we will be using for the
### principal components.  We need a set of years that is present in all of the models.
### We'll restrict the search to the r1i1p1 ensemble, since some of the others cut off
### early.
###
### NB:  I've suppressed this calculation until we get all the new data in place.  For
### the time being, we will test using hand-coded values.
# concvars <- c('tas','heatflux')
# cmip_r1_conc <- filter(cmip_individual, ensemble=='r1i1p1', variable %in% concvars)
# conc_scens <- c('historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85')
# pca_years_conc <-
#     lapply(conc_scens,
#            function(scen) {find_common_years(cmip_r1_conc, scen)})
# names(pca_years_conc) <- conc_scens
#
# esmvars <- c('tas','co2','heatflux')
# cmip_r1_esm <- filter(cmip_individual, ensemble=='r1i1p1', variable %in% esmvars)
# esm_scens <- c('esmHistorical', 'esmrcp85')
# pca_years_esm <-
#     lapply(esm_scens,
#            function(scen) {find_common_years(cmip_r1_esm, scen)})
# names(pca_years_esm) <- esm_scens

histyears_conc <-
    list(tas=1861:2005,
         heatflux=2005)
histyears_emiss <-
    list(tas=1861:2005,
         co2=1861:1990,
         heatflux=2005)

futyears_conc <-
    list(tas=2006:2100,
         heatflux=2006:2095)
futyears_emiss <-
    list(tas=2006:2099,
         co2=2006:2099,
         heatflux=2006:2090)


outdir <- 'PCA-data'

### Get principal components for the concentration-driven runs.
conc_files <- list.files('PCA-data', '^concen-', full.names=TRUE)
conc_scenlist <- lapply(conc_files, readRDS)
conc_expts <- sapply(conc_scenlist, function(d){d$experiment[1]})
## We want to do all combinations of these scenarios.  There are 31 possible
## combinations (not counting the degenerate case where none of them are
## included).  Number them from 1-15 and use the bit patterns to figure out
## which scenarios should be included.
doParallel::registerDoParallel(cores=8)
conc_vars <- 'tas'
emiss_vars <- c('tas','co2')
conc_hflux_vars <- c('tas', 'heatflux')
emiss_hflux_vars <- c('tas', 'co2', 'heatflux')

pca_all_combos <- function(name_stem, vars) {
    pcas <- foreach(i=1:31) %dopar% {
        expts_l <- as.logical(c(bitwAnd(i,1), bitwAnd(i,2), bitwAnd(i,4), bitwAnd(i,8), bitwAnd(i,16)))
        expts <- conc_expts[expts_l]
        scenlist <- conc_scenlist[expts_l]

        ## We don't use heatflux in the historical experiment, so if we have *only* the
        ## historical experiment, we need to remove heat flux from the variable list

        hectorcal::compute_pc(scenlist, expts, vars, futyears_conc, histyears_conc)
    }
    for(pca in pcas) {
        filename <- paste0(name_stem,paste0(pca$meta_data$experiment,
                                            collapse='-'), '.rds')
        saveRDS(pca, file.path(outdir, filename))
    }
    invisible(NULL)
}

## Concentration runs without heat flux
pca_all_combos('pc-conc-', conc_vars)

## Concentration runs with heat flux
pca_all_combos('pc-conc-hflux-', conc_hflux_vars)


## For the emissions run we really only want the rcp85 and historical runs for now
emiss_files <- c('PCA-data/emissCC-esmrcp85.rds', 'PCA-data/emissCC-esmHistorical.rds')
emiss_scenlist <- lapply(emiss_files, readRDS)
emiss_pca <- hectorcal::compute_pc(emiss_scenlist, c('esmrcp85', 'esmHistorical'), emiss_vars ,
                                   futyears_emiss, futyears_emiss)
saveRDS(emiss_pca, file.path(outdir, 'pc-emissCC-esmHistorical-esmrcp85.rds'))

## Also run with heat flux included
emiss_hflux_pca <- hectorcal::compute_pc(emiss_scenlist, c('esmrcp85', 'esmHistorical'),
                                         emiss_hflux_vars, futyears_emiss, histyears_emiss)
saveRDS(emiss_hflux_pca, file.path(outdir,
                                   'pc-emissCC-hflux-esmHistorical-esmrcp85.rds'))
