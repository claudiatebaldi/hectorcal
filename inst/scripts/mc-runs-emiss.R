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

### To use a restart file, specify the file stem up through the number of samples
### For example, 'testrun-1000'

mc_run_emiss <- function(runid, nsamp, filestem='hectorcal-emiss',
                         plotfile=TRUE, npc=10, restart=NULL)
{
    runid <- as.integer(runid)
    serialnumber <- bitwAnd(runid, 15)
    pcsflag <- as.logical(bitwAnd(runid, 16))
    hfflag <- as.logical(bitwAnd(runid, 32))
    meanflag <- as.logical(bitwAnd(runid, 64))
    debugflag <- as.logical(bitwAnd(runid, 128))


    ## input files for hector initialization
    rcps <- c('esmrcp85')
    inputfiles <- file.path('input', sprintf('hector_%s.ini', canonicalize_expt_names(rcps)))
    inputfiles <- system.file(inputfiles, package='hector', mustWork = TRUE)
    names(inputfiles) <- rcps

    ## get the comparison data we will be using
    if(pcsflag) {
        ## Drop principal components for which there is no variation in the comparison data.
        nkeep <- npc+1                  # Number of requested PCs, plus the
                                        # heatflux row.
        compdata <- emiss_pc_comparison[1:nkeep,]
        pcs <- pc_emiss
    }
    else  {
        ## We have all the heat flux for this one, so we need to drop the
        ## years we won't be using.  Also drop the emissions driven runs
        compdata <- dplyr::filter(esm_comparison,
                                  grepl('^esm', experiment) &
                                      (variable != 'heatflux' | (year==hfyear &
                                                                     experiment==hfexpt)))
        ## We also need to filter out data that doesn't belong, like RCP scenario
        ## data prior to 2006 or co2 data (which doesn't belong in concentration driven runs)
        compdata <- dplyr::filter(compdata,
                                  (experiment=='esmHistorical' & year > 1860 & year < 2006) |
                                    (experiment!='esmHistorical' & year > 2005 & year <= 2100))
        pcs <- NULL
    }



    ## if the heat flux flag is cleared, then drop the heatflux value
    if(!hfflag) {
        compdata <- dplyr::filter(compdata, variable != 'heatflux')
    }

    ## Create the log-posterior function
    lpf <- build_mcmc_post(compdata, inputfiles, pcs, smoothing = 0.5, cal_mean = meanflag,
                           hflux_year = hfyear, hflux_expt_regex = hfexpt)

    ## initial parameters
    if(is.null(restart)) {
        p0 <- c(3.65, 0.98, 1.65, 0.95, 0.04, 1.27, 281.6)
    }
    else {
        restartfile <- paste(restart, runid, 'mcrslt.rds', sep='-')
        p0 <- readRDS(restartfile)
    }

    ## Set an appropriate scale, depending on the protocol.
    if(pcsflag) {
        scale <- c(0.25, 0.25, 0.1, 0.1, 0.1, 0.1, 0.25)    # based on some experimentation with short runs
    }
    else {
        scale <- c(0.25, 0.25, 0.1, 0.1, 0.05, 0.2, 0.05)/2
    }
    pnames <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE(), BETA(), Q10_RH(), PREINDUSTRIAL_CO2())
    if(meanflag) {
        ## Add sigma parameters for mean calibration
        if(pcsflag) {
            if(is.null(restart)) {
                p0 <- c(p0, 0.05)
            }
            sclfac <- c(scale/2, 0.1)
            scale <- diag(nrow=length(sclfac), ncol=length(sclfac))
            scale[1,2] <- scale[2,1] <- 0.9
            scale[5,6] <- scale[6,5] <- 0.9
            scale <- cor2cov(scale, sclfac)
            pnames <- c(pnames, 'sig')
        }
        else {
            if(is.null(restart)) {
                p0 <- c(p0, 0.05, 1.0)      # sigma for temperature and co2
            }
            tempscl <- 1e-4
            co2scl <- 1e-4
            sclfac <- c(scale/20, tempscl, co2scl)
            ## S and kappa distributions are very narrow, so that step needs to be small
            sclfac[1] <- sclfac[1] / 10
            sclfac[2] <- sclfac[2] / 10
            scale <- diag(nrow=length(sclfac), ncol=length(sclfac))
            ## many of the parameters turn out to be highly correlated in this case.
            scale[1,2] <- scale[2,1] <- 0.95
            scale[1,3] <- scale[3,1] <- -0.7
            scale[1,4] <- scale[4,1] <- 0.3
            scale[2,3] <- scale[3,2] <- -0.7
            scale[2,4] <- scale[4,2] <- 0.3
            scale[3,4] <- scale[4,3] <- -0.25
            scale[1,5] <- scale[5,1] <- -0.3
            scale[2,5] <- scale[5,2] <- -0.3
            scale <- cor2cov(scale, sclfac)
            pnames <- c(pnames, 'sigt', 'sigco2')
        }

        if(hfflag) {
            if(is.null(restart)) {
                p0 <- c(p0, 5.0)
            }
            hfscl <- 1e-1
            if(is.matrix(scale)) {
                newscale <- matrix(0, nrow=nrow(scale)+1, ncol=ncol(scale)+1)
                newscale[1:nrow(scale), 1:ncol(scale)] <- scale
                scale <- newscale
                scale[nrow(scale),ncol(scale)] <- hfscl
            }
            else {
                scale <- c(scale, hfscl)
            }
            pnames <- c(pnames, 'sighf')
        }
    }
    else {
        ## Not doing mean calibration.  Still want to add some correlation.
        sclfac <- scale
        scale <- diag(nrow=length(sclfac), ncol=length(sclfac))
        scale[1,2] <- scale[2,1] <- 0.95
        scale[5,6] <- scale[6,5] <- 0.9
        scale <- cor2cov(scale, sclfac)
    }



    if(is.null(restart)) {
        names(p0) <- pnames
    }
    if(is.vector(scale)) {
        names(scale) <- pnames
    }
    else {
        rownames(scale) <- colnames(scale) <- pnames
    }


    ## Run monte carlo
    foreach::registerDoSEQ()
    set.seed(867-5309 + serialnumber)

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

    if(nsamp < 10000) {
        plotfilename <- paste(filestem, runid, 'diagplot.png', sep='-')
        msc <- metrosamp2coda(list(ms))
        if(plotfile) {
            png(plotfilename, 1024, 1024)
        }
        plot(msc)
        if(plotfile) {
            dev.off()
        }
    }
    savefilename <- paste(filestem, runid, 'mcrslt.rds', sep='-')

    saveRDS(ms, savefilename, compress='xz')

    invisible(ms)
}
