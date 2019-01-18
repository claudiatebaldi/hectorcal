## set the seed for reproducibility
set.seed(867-5309)


## initial guess
pv0 <- c(3, 1, 2.3, 0.3, 2, 285, 1.0, 10.0)
## mcmc scale for carbon cycle + mean (seems to be effective when scaled by a factor of 0.25)
mcmc_scale <- c(0.2, 0.2, 0.2, 0.1, 0.2, 2.0, 0.1, 2.0)

