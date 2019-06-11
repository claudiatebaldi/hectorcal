library(hector)
library(hectorcal)
library(metrosamp)
library(doParallel)

registerDoParallel(cores=4)

## set up the likelihood function
inputdir <- system.file('input', package='hector')
rcps <- c('rcp26', 'rcp45', 'rcp60', 'rcp85')
inputfiles <- file.path(inputdir, sprintf('hector_%s.ini', rcps))
names(inputfiles) <- rcps
lpf <- build_mcmc_post(conc_pc_comparison, inputfiles, pc_conc, cal_mean = FALSE)

p0 <- c(3.0, 2.3, 1.0, 1.0)
names(p0) <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE())

scale <- c(0.75, 1.25, 0.5, 0.5)
names(scale) <- names(p0)

set.seed(867-5309)

ms1 <- metrosamp(lpf, p0, 500, 1, scale)
