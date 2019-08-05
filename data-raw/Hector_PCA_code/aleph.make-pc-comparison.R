## Prepare a hybrid dataset for likelihood function calculations
##
## Combine the principal components based on temperature and (if applicable) CO2
## with year 2100 values for ocean heat flux.  The result will be a table if
## min/mean/max values for the PCs and for the year 2100 heat flux values for
## the rcp85 experiment (only). The format will be the same as esm_comparison.
##
## This data reflects our current plan for doing the envelope calibration for the
## model.  It should be suitable for passing as the compdata argument to build_mcmc_post.
##
## This function assumes that all of the package data for the CMIP ensemble results
## and Hector principal components has been added to the package.  At this stage of
## the development, they are, but if you wanted to rebuild the package data from
## scratch, there would be several stages of running data scripts, followed by
## package reloads.  C'est la guerre.
##
## params:
##   pctbl - Table of ESM principal component projections.  Generally either cmip_conc_pcproj
##           or cmip_emiss_pcproj
##   heatfluxexpt - Name of the experiment to find the heatflux in.  Generally either
##           rcp85 or esmrcp85.  The data will be pulled from the esm_comparison table.
##
library('magrittr')
make_pc_comparison <- function(pctbl, heatfluxexpt)
{
    pc_comparison <- dplyr::group_by(pctbl, PC) %>%
        dplyr::summarise(mina=min(value), maxb=max(value),
                         a10=quantile(value, 0.1), b90=quantile(value, 0.9),
                         cmean=mean(value), cmedian=median(value)) %>%
        dplyr::rename(variable=PC) %>%
        dplyr::mutate(variable = paste0('PC',variable))
    hflux <- dplyr::filter(hectorcal::esm_comparison, variable=='heatflux', year==2100, experiment==heatfluxexpt) %>%
        dplyr::select(variable, mina, maxb, a10, b90, cmean, cmedian)
    dplyr::bind_rows(hflux, pc_comparison)
}

## The result will be stored in the package variables conc_pc_comparison and emiss_pc_comparison.
conc_pc_comparison <- make_pc_comparison(cmip_conc_pcproj, 'rcp85')
usethis::use_data(conc_pc_comparison, compress='xz')
emiss_pc_comparison <- make_pc_comparison(cmip_emiss_pcproj, 'esmrcp85')
usethis::use_data(emiss_pc_comparison, compress='xz')
