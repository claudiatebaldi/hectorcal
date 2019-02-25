#### Experiment 1: Do the results of full-range calibration actually cover the
####               full range of ESM results?  If so, is the apparent falloff in
####               density toward te edge of the range real, and what could be
####               causing it?
#### Part B:  Run a calibration to the top 10% of the final temperature range.
####          We will still allow the full CO2 range.


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

pnames <- c('S', 'aero', 'kappa', 'beta', 'q10', 'pre_co2')

comp_esmrcp85_hi10 <-
    filter(esm_comparison, experiment=='esmrcp85', year==2100) %>%
      mutate(a=if_else(variable=='co2', mina, b90), b=maxb)

esmrcp85_ini <- system.file("input/hector_rcp85.ini", package = "hector")

nchain <- 8
ccounter <- seq(1,nchain)               # chain counter

lposts <-
    lapply(ccounter, function(ichain) {
               lp <- build_mcmc_post(comp_esmrcp85_hi10, esmrcp85_ini,
                                     lowcol='a', hicol='b', years=2100,
                                     cal_mean=FALSE) # nchain identical cores
               function(p) {suppressMessages(lp(p))}
           })

nsamp <- 10000
bsize <- 1
nwarm <- 1000

doParallel::registerDoParallel(cores=nchain)

mcrslts_hi10 <-
    foreach(ichain=ccounter) %dopar% {
        set.seed(867-5309 + ichain)
        warmup <- metrosamp(lposts[[ichain]], pv0, nwarm, 1, mcmc_scale)
        r <- metrosamp(lposts[[ichain]], warmup, nsamp, bsize)
        colnames(r$samples) <- pnames
        r
}

## Convert to coda::mcmc.list to run diagnostics
mcmc_hi10 <- coda::mcmc.list(
    lapply(mcrslts_hi10, function(x) {coda::mcmc(x$samples)}))


saveRDS(mcrslts_hi10, 'mcrslts_hi10.rds')

doParallel::stopImplicitCluster()
