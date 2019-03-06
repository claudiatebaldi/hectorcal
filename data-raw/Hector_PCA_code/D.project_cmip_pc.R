#### Project the CMIP models onto the Hector principal components


### This function will process each group (emissions-forced
### vs. concentration-forced) of model data into PC projections.
pcify <- function(cmip, pca, msgtag)
{
    ## split the data by model
    cmip_models <- split(cmip, cmip$model)

    ## Filter to just the experiments required by the selected PCA
    cmip_models <- lapply(cmip_models,
                          function(d){dplyr::filter(d, experiment %in%
                                                    pca$meta_data$experiment)})

    ## Some models are missing some of the data needed for the comparison.  Drop
    ## these.
    isok <- sapply(cmip_models, function(x){chkdata(x,pca)})
    badmodels <- sapply(cmip_models[!isok], function(d){d$model[1]})
    message(msgtag, 'Dropping the following models due to incomplete data: ',
            paste(badmodels, collapse=', '))
    cmip_models <- cmip_models[isok]

    ## Some of the cmip inputs have duplicated rows.
    cmip_models <- lapply(cmip_models, unique)

    ## For each model, project onto the principal components
    projections <- lapply(cmip_models, hectorcal::project_climate,
                          principal_components=pca, row_vector=FALSE)

    ## reformat into a data frame
    dfprojections <- lapply(names(projections),
                            function(name) {
                                p <- projections[[name]]
                                index <- seq_along(p)
                                data.frame(model=name, PC=index, value=p, stringsAsFactors = FALSE)
                            })
    dplyr::bind_rows(dfprojections)
}

### Check to see whether a data set has all the data required for this
### projection
chkdata <- function(d, pca)
{
    if(!all(pca$meta_data$experiment %in% d$experiment)) {
        return(FALSE)
    }

    hasvars <-
        dplyr::group_by(d, experiment) %>%
          dplyr::summarise(has_all_vars = all(pca$meta_data$variable %in%
                           variable))
    if(!all(hasvars$has_all_vars)) {
        return(FALSE)
    }

    hasyears <-
        dplyr::mutate(d, ishist=grepl('[Hh]istorical', experiment)) %>%
        dplyr::group_by(experiment, variable) %>%
        dplyr::summarise(has_all_years =
                           (all(ishist) && all(pca$meta_data$histyear %in% year)) ||
                           (all(!ishist) && all(pca$meta_data$year %in% year)))
    if(!all(hasyears$has_all_years)) {
        return(FALSE)
    }
    else {
        return(TRUE)
    }
}

### Join the historical data up to each of the future datasets
histjoin <- function(d)
{
    ishist <- grepl('[Hh]istorical', d$experiment)
    historical <- unique(d[ishist, ])   # some of the rows in the historical
                                        # data were duplicated for some reason
    future <- unique(d[!ishist, ])      # haven't seen any indication that these
                                        # were duplicated, but better safe than
                                        # sorry.

    future_l <- split(future, future$experiment)
    future_l <-
        lapply(future_l,
               function(d) {
                   dplyr::bind_rows(dplyr::mutate(historical,
                                                  experiment=d$experiment[1]),
                                    d)
               })
    dplyr::bind_rows(future_l)
}

### Load the data we will need to use
datadir <- 'PCA-data'
pc_conc_all <- readRDS(file.path(datadir,
                                 'pc-conc-historical-rcp26-rcp45-rcp60-rcp85.rds'))
pc_emiss_rcp85 <- readRDS(file.path(datadir, 'pc-emissCC-esmHistorical-esmrcp85.rds'))

## Ensemble members other than r1i1p1 are a bit spotty.  For many models they
## don't have a full set of scenarios, or don't have a historical run, or are
## otherwise incomplete.  Save ourselves some headaches by restricting our
## analysis to just the r1i1p1 runs.

cmip_data <- dplyr::filter(hectorcal::cmip_individual, ensemble=='r1i1p1')
cmip_data <- split(cmip_data, grepl('^esm', cmip_data$experiment))

cmip_conc_pcproj <- pcify(cmip_data[["FALSE"]], pc_conc_all, 'conc runs: ')
cmip_emiss_pcproj <- pcify(cmip_data[["TRUE"]], pc_emiss_rcp85, 'emiss runs: ')

usethis::use_data(cmip_conc_pcproj, overwrite=TRUE, compress='xz')
usethis::use_data(cmip_emiss_pcproj, overwrite=TRUE, compress='xz')

