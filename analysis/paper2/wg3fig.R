set.seed(867-5309)

library(foreach)
library(doParallel)
library(hector)
library(hectorcal)
library(metrosamp)
library(dplyr)
library(ggplot2)

datadir_conc <- file.path('analysis','mcmc','conc', 'primary')   # XXX Temporary, replace with final runs when available
mcruns_conc <- proc_mc_rslts(datadir_conc, 'hectorcal-conc')
datadir_emiss <- file.path('analysis', 'mcmc', 'emiss', 'primary')
mcruns_emiss <- proc_mc_rslts(datadir_emiss, 'hectorcal-emiss')

figdir <- file.path('analysis', 'figs-paper2')

get_infile <- function(rcp, constrained=TRUE)
{
    indir <- system.file('input', package='hector')
    if(constrained) {
        basefiles <- paste0('hector_', rcp, '_constrained.ini')
    }
    else {
        basefiles <- paste0('hector_', rcp, '.ini')
    }
    file.path(indir, basefiles)
}

rcps_conc <- paste0('rcp', c(26,45,60,85))
inifiles_conc <- get_infile(rcps_conc, TRUE)
rcps_emiss <- paste0('rcp', c(85))
inifiles_emiss <- get_infile(rcps_emiss, FALSE)

hcores_conc <- mapply(newcore, inifile=inifiles_conc, name=rcps_conc)
hcores_emiss <- mapply(newcore, inifile=inifiles_emiss, name=rcps_emiss)
hcores <- c(hcores_conc, hcores_emiss)
scenario_names <- c(rcps_conc, paste0(rcps_emiss,'*'))
ncore_conc <- length(hcores_conc)
ncore <- length(hcores)

setparm <- function(hcore, parms)
{
    pnames <- names(parms)
    for(name in pnames) {
        #cat(paste('pname: ', name, 'value: ', parms[name], '\n'))
        setvar(hcore, NA, name, parms[name], getunits(name))
    }
}

### Run the concentration runs
registerDoParallel(cores=ncore)
nsample <- 1000
samps_conc <- rsample(mcruns_conc$mcobjs[['48']], nsample)
samps_emiss <- rsample(mcruns_emiss$mcobjs[['48']], nsample)
nparm_conc <- 4
nparm_emiss <- 7

tbase <- 0.61     # Based on the caption to figure 6.13, WG-III AR5

tincr <-
    foreach(icore=1:ncore) %dopar% {
        hcore <- hcores[[icore]]
        if(icore>ncore_conc) {
            samps <- samps_emiss[,1:nparm_emiss]
        }
        else {
            samps <- samps_conc[, 1:nparm_conc]
        }
        foreach(isamp=1:nsample, .combine=c) %do% {
            setparm(hcore, samps[isamp,])
            reset(hcore)
            run(hcore, 2100)
            ##tbase <- fetchvars(hcore, 1986:2005, GLOBAL_TEMP())
            t2100 <- fetchvars(hcore, 2100, GLOBAL_TEMP())
            ##t2100$value - mean(tbase$value)
            t2100$value - tbase
        }
    }
names(tincr) <- scenario_names
tistats <- t(sapply(tincr, quantile, probs=c(0.0, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95, 1.0)))
tistats_d <- as_tibble(tistats) %>% mutate(RCP=rownames(tistats))

wg3plt <-
    ggplot(data=tistats_d, aes(x=RCP)) +
    geom_boxplot(aes(ymin=`5%`, lower=`16%`, middle=`50%`, upper=`84%`, ymax=`95%`), stat='identity',
                 fill='SteelBlue', color='MidnightBlue') +
    xlab('') + ylab('Temperature Increase (\u00B0C)') +
    theme_bw(base_size = 12)

ggsave(file.path(figdir,'wg3plt.pdf'), plot=wg3plt, device='pdf', width=4, height=4, units='in')

