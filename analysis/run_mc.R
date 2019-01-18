#### Run a monte carlo test calculation.

library('hectorcal')
library('hector')
library('metrosamp')
library('dplyr')

## set the seed for reproducibility
set.seed(867-5309)


## initial guess
pv0 <- c(3, 1, 2.3, 0.3, 2, 285, 1.0, 10.0)
## mcmc scale for carbon cycle + mean (seems to be effective when scaled by a factor of 0.5)
mcmc_scale <- c(0.2, 0.2, 0.2, 0.1, 0.2, 2.0, 0.1, 2.0)

comp_esmrcp85 <- filter(esm_comparison, experiment=='esmrcp85')
esmrcp85_ini <- system.file("input/hector_rcp85.ini", package = "hector")

lpost <- build_mcmc_post(comp_esmrcp85, esmrcp85_ini)

mcrslt <- metrosamp(lpost, pv0, 100, 1, mcmc_scale*0.5)
