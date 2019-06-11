library(hector)
library(hectorcal)
library(metrosamp)
library(doParallel)
library(ggplot2)
library(magrittr)


hfyear <- 2100
hfexpt <- 'rcp85'

#### Definition of runid:
### Bits 0-3: Serial number (used to give each run of the same configuration a
###           unique RNG seed
### Bit 4: Use pcs flag.  If false, we use the outputs directly
### Bit 5: Use heat flux flag.
### Bit 6: Mean calibrtion flag
### Bit 7: Debug flag

production_run_conc <- function(runid, nsamp, filestem='hectorcal',
                                plotfile=TRUE, npc=10)
{
    runid <- as.integer(runid)
    serialnumber <- bitwAnd(runid, 15)
    pcsflag <- as.logical(bitwAnd(runid, 16))
    hfflag <- as.logical(bitwAnd(runid, 32))
    meanflag <- as.logical(bitwAnd(runid, 64))
    debugflag <- as.logical(bitwAnd(runid, 128))


    ## input files for hector initialization
    rcps <- c('rcp26', 'rcp45', 'rcp60', 'rcp85')
    inputfiles <- file.path('input', sprintf('hector_%s.ini', rcps))
    inputfiles <- system.file(inputfiles, package='hector', mustWork = TRUE)
    names(inputfiles) <- rcps

    ## get the comparison data we will be using
    if(pcsflag) {
        ## Drop principal components for which there is no variation in the comparison data.
        nkeep <- npc+1                  # Number of requested PCs, plus the
                                        # heatflux row.
        compdata <- conc_pc_comparison[1:nkeep,]
        pcs <- pc_conc
    }
    else  {
        ## We have all the heat flux for this one, so we need to drop the
        ## years we won't be using.  Also drop the emissions driven runs
        compdata <- dplyr::filter(esm_comparison,
                                  !grepl('^esm', experiment) &
                                      (variable != 'heatflux' | (year==hfyear &
                                                                     experiment==hfexpt)))
        ## We also need to filter out data that doesn't belong, like RCP scenario
        ## data prior to 2006 or co2 data (which doesn't belong in concentration driven runs)
        compdata <- dplyr::filter(compdata,
                                  (experiment=='historical' & year > 1860 & year < 2006) |
                                    (experiment!='historical' & year > 2005 & year <= 2100)) %>%
            dplyr::filter(variable != 'co2')
        pcs <- NULL
    }



    ## if the heat flux flag is cleared, then drop the heatflux value
    if(!hfflag) {
        compdata <- dplyr::filter(compdata, variable != 'heatflux')
    }

    ## Create the log-posterior function
    lpf <- build_mcmc_post(compdata, inputfiles, pcs, cal_mean = meanflag,
                           hflux_year = hfyear, hflux_expt_regex = hfexpt)

    ## initial parameters
    p0 <- c(3.0, 2.3, 1.0, 1.0)
    scale <- c(0.75, 1.25, 0.5, 0.5)    # based on some experimentation with short runs
    pnames <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE())
    if(meanflag) {
        ## Add sigma parameters for mean calibration
        if(pcsflag) {
            p0 <- c(p0, 0.2)
            scale <- c(scale/4, 0.025)
            pnames <- c(pnames, 'sig')
        }
        else {
            p0 <- c(p0, 1.0, 5.0)
            scale <- c(scale, 0.5, 2.0)
            pnames <- c(pnames, 'sigt', 'sigco2')
        }

        if(hfflag) {
            p0 <- c(p0, 0.5)
            scale <- c(scale, 0.05)
            pnames <- c(pnames, 'sighf')
        }
    }
    names(p0) <- pnames
    names(scale) <- pnames


    ## Run monte carlo
    registerDoParallel(cores=4)
    set.seed(867-5309 + serialnumber)


    names(scale) <- names(p0)

    ms <- metrosamp(lpf, p0, nsamp, 1, scale, debug=debugflag)

    ## save output and make some diagnostic plots
    cat('Runid:\n')
    print(runid)
    cat('output file stem:\n')
    print(filestem)
    cat('Acceptance fraction:\n')
    print(ms$accept)
    cat('Effective number of samples:\n')
    print(neff(ms$samples))

    plotfilename <- paste(filestem, runid, 'diagplot.png', sep='-')
    msc <- metrosamp2coda(list(ms))
    if(plotfile) {
        png(plotfilename, 1024, 1024)
    }
    plot(msc)
    if(plotfile) {
        dev.off()
    }
    savefilename <- paste(filestem, runid, 'mcrslt.rds', sep='-')

    saveRDS(ms, savefilename, compress='xz')

    invisible(ms)
}
