#### Experiment 1: Do the results of full-range calibration actually cover the
####               full range of ESM results?  If so, is the apparent falloff in
####               density toward te edge of the range real, and what could be
####               causing it?
#### Part B:  Run the full-range calibration, with very wide priors.  The
####          purpose here is to see whether the edge density falloff is due to
####          the influence of the priors.


library('hectorcal')
library('hector')
library('metrosamp')
library('dplyr')
library('foreach')


## initial guess
pv0 <- c(3, 1, 2.3, 0.3, 2, 285)
## mcmc scale.  Preliminary experiment runs indicated that this gave us a good
## acceptance rate and reasonable effective sample size.
mcmc_scale <- c(0.4, 0.8, 0.4, 0.4, 0.8, 10.0)

plabels <- c('S', 'aero', 'kappa', 'beta', 'q10', 'pre_co2')

comp_esmrcp85 <- filter(esm_comparison, experiment=='esmrcp85')
esmrcp85_ini <- system.file("input/hector_rcp85.ini", package = "hector")
prior_params <- c(ecssig=10, aerosig=10, kappasig=10, betasig=10, q10sig=10,
                  c0sig=60)


nchain <- 8
ccounter <- seq(1,nchain)               # chain counter

lposts <-
    lapply(ccounter, function(ichain) {
               lp <- build_mcmc_post(comp_esmrcp85, esmrcp85_ini,
                                     cal_mean=FALSE,
                                     prior_params=prior_params,
                                     use_lnorm_ecs=FALSE) # nchain identical cores
               function(p) {suppressMessages(lp(p))}
           })

nsamp <- 10000
bsize <- 1
nwarm <- 500

doParallel::registerDoParallel(cores=nchain)

mcrslts_wide_priors <-
    foreach(ichain=ccounter) %dopar% {
        set.seed(867-5309 + ichain)
        warmup <- metrosamp(lposts[[ichain]], pv0, nwarm, 1, mcmc_scale)
        r <- metrosamp(lposts[[ichain]], warmup, nsamp, bsize)
        colnames(r$samples) <- plabels
        r
}

## Convert to coda::mcmc.list to run diagnostics
mcmc_wide_priors <- coda::mcmc.list(
    lapply(mcrslts_wide_priors, function(x) {coda::mcmc(x$samples)}))


saveRDS(mcrslts_wide_priors, 'mcrslts_wide_priors.rds')

doParallel::stopImplicitCluster()
