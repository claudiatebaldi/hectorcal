set.seed(867-5309)

library(foreach)
library(doParallel)
library(hector)
library(dplyr)

## Load MCMC data
## XXX Using data from pilot runs.  Update to full when available
datadir_conc <- file.path('analysis','mcmc','conc', 'primary')   # XXX Temporary, replace with final runs when available
mcruns_conc <- proc_mc_rslts(datadir_conc, 'hectorcal-conc')

figdir <- file.path('analysis', 'figs-paper2')

## Set up hector cores.
##  inifiles : list of input files
##  names : list of scenario names
##  N : Number of copies of the list to make
##
##  example:
##     setup_cores(c('rcp45.ini','rcp85.ini'), c('rcp45', 'rcp85'), 3)
##  Makes 6 cores, in 3 lists of 2, each list of 2 has an rcp45 and an rcp85 core.
setup_cores <- function(inifile, name, N)
{
    lapply(1:N, function(i) {newcore(inifile=inifile, name=name)})
}

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

rcps <- paste0('rcp', c(85, 26, 45, 60))
protos <- as.character(c(48,32))
pnames <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE())

letters <- gtools::chr(seq(65, 90))
lindx <- 0

data2100 <- list()
length(data2100) <- length(protos) * length(rcps)
names(data2100) <- rep('', length(data2100))

esmdata <- filter(esm_comparison, variable=='tas', year %in% seq(1860, 2100, 1)) %>%
    mutate(variable=GLOBAL_TEMP())

for(rcp in rcps) {
    hcores <- setup_cores(get_infile(rcp), rcp, 8)
    esmrcp <- filter(esmdata, (experiment==rcp & year>2005) | (experiment=='historical' & year < 2006))
    for(proto in protos) {

        lindx <- lindx + 1
        plt <- spaghetti_plot(mcruns_conc$mcobjs[[proto]], 1000, hcores, pnames, GLOBAL_TEMP()) +
            ylab('Global Mean Temperature') +
            theme_bw(base_size=12) +
            geom_ribbon(data=esmrcp, mapping=aes(x=year, ymin=mina, ymax=maxb),
                          fill='grey10', inherit.aes=FALSE, alpha=0.5) +
            labs(title=paste('Protocol', proto), tag=letters[lindx]) +
            theme(strip.background = element_blank(), strip.text = element_blank())
        basefile <- paste0('spaghetti_', rcp, '_', proto, '.pdf')
        filename <- file.path(figdir, basefile)
        ggsave(filename, plot=plt, device='pdf', width=4, height=4, units='in')
        data2100[[lindx]] <- filter(plt$data, year==2100) %>% mutate(protocol=proto, experiment=rcp)
        names(data2100)[lindx] <- paste0(proto,'_',rcp)
    }
}

alldata2100 <- bind_rows(data2100)
esmx <- filter(esm_comparison, variable=='tas') %>% mutate(variable='Tgav') %>% select(variable, year, experiment, mina, maxb)
protosum <- group_by(alldata2100, protocol, variable, year, experiment) %>%
    left_join(esmx, by=c('year', 'variable','experiment')) %>%
    summarise(min=min(value), max=max(value), q10=quantile(value, probs=0.1), q90=quantile(value, probs=0.9), median=median(value),
              `P>5`=sum(value>5)/1000, `P<3`=sum(value<3)/1000, `p<2`=sum(value<2)/1000,
              foutside=sum(value<mina | value > maxb)/1000)

print(protosum)
