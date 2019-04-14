context('log-posterior functions')

library('dplyr')
library('hector')

test_that('log prior function works properly', {
    test_prior_params <- c(
        ecsmu=1.0, ecssig=1.0,
        aeromu=1.0, aerosig=1.0,
        volmu=1.0, volsig=1.0,
        kappamu=1.0, kappasig=1.0,
        betamu=1.0, betasig=1.0,
        q10mu=1.0, q10sig=1.0,
        c0mu=1, c0sig=1.0,
        sigtscale=1.0,
        sigco2scale=1.0)

    testlp <- make_logprior(test_prior_params, FALSE)

    ## Normally distributed
    konst1 <- stats::dnorm(1, 1, 1, log=TRUE)
    for(parm in c(hector::ECS(), hector::AERO_SCALE(), hector::VOLCANIC_SCALE(),
                  hector::DIFFUSIVITY(), hector::PREINDUSTRIAL_CO2())) {
        p <- 1
        names(p) <- parm
        expect_equal(testlp(p), konst1, info=paste('parm=',parm))
    }
    ## Truncated normal
    konst2 <- mktruncnorm(0, Inf, 1, 1)(1)
    for(parm in c(hector::BETA(), hector::Q10_RH())) {
        p <- 1
        names(p) <- parm
        expect_equal(testlp(p), konst2, info=paste('parm=',parm))
    }
    ## Half-Cauchy.  We didn't actually renormalize this one for the
    ## truncation.  Also, it's always centered at zero.
    konst3 <- stats::dcauchy(1, 0, 1, log=TRUE)
    for(parm in c('sigt', 'sigco2')) {
        p <- 1
        names(p) <- parm
        expect_equal(testlp(p), konst3, info=paste('parm=',parm))
    }

    ## Now test that if we pick one of each we get the sum
    p <- rep(1,3)
    names(p) <- c(hector::AERO_SCALE(), hector::BETA(), 'sigt')
    expect_equal(testlp(p), konst1+konst2+konst3)

    ## Check that when we turn on the log-normal distribution that works as
    ## expected.  Note ecsmu is the center for ECS, dlnorm wants the log of
    ## that.
    testlp2 <- make_logprior(test_prior_params, TRUE)
    p <- 1
    names(p) <- hector::ECS()
    konst4 <- stats::dlnorm(1, 0, 1, log=TRUE)
    expect_equal(testlp2(p), konst4)
    p <- rep(1,4)
    names(p) <- c(hector::ECS(), hector::VOLCANIC_SCALE(), hector::Q10_RH(),
                  'sigco2')
    expect_equal(testlp2(p), konst1+konst2+konst3+konst4)
})


test_that('log-likelihood with output comparisons works', {
    ## First we need to make some comparison data.  Run a Hector calculation
    ## with known parameters to get some baseline data.
    nhectorparm <- 7
    parms <- c(2.5, 2.5, 1.0, 1.0, 0.5, 2.0, 280.0, 1.0, 1.0)
    names(parms) <- c(hector::ECS(), hector::DIFFUSIVITY(),
                      hector::AERO_SCALE(), hector::VOLCANIC_SCALE(),
                      hector::BETA(), hector::Q10_RH(),
                      hector::PREINDUSTRIAL_CO2(),
                      'sigt','sigco2')
    ini45 <- system.file('input/hector_rcp45.ini', package='hector')
    ini85 <- system.file('input/hector_rcp85.ini', package='hector')
    hcore45 <- newcore(ini45, name='rcp45')
    parameterize_core(parms[1:nhectorparm], hcore45)
    histyears <- 2000:2005
    futyears <- 2006:2010
    run(hcore45, 2010)
    ## grab data from 1990-2000 as "historical" and 2001-2010 as "future".
    ## Might as well use temperature and CO2 as the comparison variables.
    dh <- fetchvars(hcore45, histyears, c(GLOBAL_TEMP(), ATMOSPHERIC_CO2()),
                    'historical')
    d45 <- fetchvars(hcore45, futyears, c(GLOBAL_TEMP(), ATMOSPHERIC_CO2()))
    shutdown(hcore45)

    ## Get some rcp85 too so that we can test having multiple future scenarios
    hcore85 <- newcore(ini85, name='rcp85')
    parameterize_core(parms[1:nhectorparm], hcore85)
    run(hcore85, 2010)
    d85 <- fetchvars(hcore85, futyears, c(GLOBAL_TEMP(), ATMOSPHERIC_CO2()))
    shutdown(hcore85)

    compdata <-
        bind_rows(dh,d45,d85) %>%
          rename(experiment=scenario, cmean=value) %>%
          mutate(cmedian=cmean, mina=cmean-2, maxb=cmean+2, a10=cmean-1,
                 b90=cmean+1) %>%
          mutate(variable=if_else(variable=='Tgav', 'tas', 'co2')) %>%
          select(year, variable, experiment, mina, maxb, a10, b90, cmean,
                 cmedian) %>%
          arrange(year)                 # This last one is to test that
                                        # scrambled data gets sorted out.

    ## Now we've got comparison data.  Sort of.  It's hard to predict how
    ## changing the parameters will change the output, so we'll always run with
    ## the same set of parameters.  Instead, we'll make new likelihood functions
    ## with the comparison data changed in specific, controlled ways, and we'll
    ## verify that the effect on the likelihood output is what we anticipate.

    doParallel::registerDoParallel()

    ## 1. Mean calibration, perfect model match.
    inifiles <- c(rcp85=ini85, rcp45=ini45)
    llfun1 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata, 0.1, 'maxb','mina',NULL)
    out <- expect_silent(llfun1(parms))
    expected <- nrow(compdata)*stats::dnorm(0,0,parms['sigt'], log=TRUE) # assumes sigt==sigco2
    expect_equal(out, expected)

    ## 2. Mean calibration, discrepancy in the rcp85 scenario
    compdata2 <- dplyr::mutate(compdata,
                               cmean=dplyr::if_else(experiment!='rcp85', cmean,
                               dplyr::if_else(variable=='tas', cmean-1, cmean+1)))
    llfun2 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata2, 0.1,'maxb','mina', NULL)
    out <- expect_silent(llfun2(parms))
    k1 <- stats::dnorm(0,0,parms['sigt'], log=TRUE)
    k2 <- stats::dnorm(1,0,parms['sigt'], log=TRUE)
    expected <- k1*(nrow(dh)+nrow(d45)) + k2*nrow(d85)
    out <- expect_silent(llfun2(parms))
    expect_equal(out, expected)


    ## 3. Envelope calibration, all values in range
    llfun3 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata, 0.1, 'maxb','mina',NULL)
    out <- expect_silent(llfun3(parms))
    expected <- nrow(compdata) * log(mesa(2, 0, 4, 0.4))
    expect_equal(out, expected)

    ## 4. Envelope calibration, all values on the edge of the range
    compdata4 <- dplyr::mutate(compdata,
                               mina=mina+2, maxb=maxb+2)
    llfun4 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata4, 0.1,
                                 'maxb','mina', NULL)
    expected <- nrow(compdata) * log(mesa(0, 0, 4, 0.4))
    out <- expect_silent(llfun4(parms))
    expect_equal(out, expected)

    doParallel::stopImplicitCluster()
})


test_that('log-likelihood with PCA comparison works', {
    pcs <- readRDS('pc-conc-historical-rcp45-rcp85.rds')
    histyears <- pcs$meta_data$histyear
    futyears <- pcs$meta_data$year
    npc <- 10
    ini45 <- system.file('input/hector_rcp45_constrained.ini', package='hector')
    ini85 <- system.file('input/hector_rcp85_constrained.ini', package='hector')

    nhectorparm <- 4
    parms <- c(2.5, 2.5, 1.0, 1.0, 1.0)
    names(parms) <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE(),
                      'sig')
    hcore45 <- newcore(ini45, name='rcp45')
    parameterize_core(parms[1:nhectorparm], hcore45)
    run(hcore45,2100)

    dh <- fetchvars(hcore45, histyears, GLOBAL_TEMP(),
                    'historical')
    d45 <- fetchvars(hcore45, futyears, GLOBAL_TEMP())
    shutdown(hcore45)

    hcore85 <- newcore(ini85, name='rcp85')
    parameterize_core(parms[1:nhectorparm], hcore85)
    run(hcore85, 2100)
    d85 <- fetchvars(hcore85, futyears, GLOBAL_TEMP())
    shutdown(hcore85)

    scendata <- bind_rows(dh,d45,d85) %>% rename(experiment=scenario) %>%
      mutate(variable=hvar2esmvar(variable))
    pcproj <- project_climate(scendata, pcs, FALSE)
    pcproj <- pcproj[1:npc]

    compdata <- data.frame(PC=seq_along(pcproj), cmean=pcproj, mina=pcproj-2,
                           maxb=pcproj+2)

    ## See notes in the previous test; we have 4 cases to run

    doParallel::registerDoParallel()
    inifiles <- c(rcp85=ini85, rcp45=ini45)
    ## 1. Mean calibration, perfect model match
    llfun1 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata, 0.1, 'maxb',
                                 'mina', pcs)
    out <- expect_silent(llfun1(parms))
    expected <- nrow(compdata)*stats::dnorm(0,0,parms['sig'], log=TRUE) # assumes sigt==sigco2
    expect_equal(out, expected)

    ## 2. Mean calibration, imperfect model match
    compdata2 <- dplyr::mutate(compdata, cmean = dplyr::if_else(PC%%2==0, cmean+1,
                                         cmean-1))
    llfun2 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata2, 0.1,
                                 'maxb','mina',pcs)
    expected <- nrow(compdata)*stats::dnorm(1, 0, parms['sig'], log=TRUE)
    out <- expect_silent(llfun2(parms))
    expect_equal(out, expected)

    ## 3. Envelope callibration, centered
    llfun3 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata, 0.1, 'maxb',
                                 'mina', pcs)
    expected <- nrow(compdata) * log(mesa(2, 0, 4, 0.4))
    out <- expect_silent(llfun3(parms))
    expect_equal(out, expected)

    ## 4. Envelope calibration, edge
    compdata4 <- dplyr::mutate(compdata, mina=dplyr::if_else(PC%%2==0, mina+2,
                                         mina-2),
                               maxb=dplyr::if_else(PC%%2==0, maxb+2, maxb-2))
    llfun4 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata4, 0.1, 'maxb',
                                 'mina', pcs)
    expected <- nrow(compdata) * log(mesa(0,0,4,0.4))
    out <- expect_silent(llfun4(parms))
    expect_equal(out, expected)

    doParallel::stopImplicitCluster()
})
