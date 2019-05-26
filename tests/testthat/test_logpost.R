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
        sigco2scale=1.0,
        sigscale=1.0)

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
    for(parm in c('sigt', 'sigco2', 'sig')) {
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

    compdata <- readRDS('out_compdata.rds')

    ## Now we've got comparison data.  Sort of.  It's hard to predict how
    ## changing the parameters will change the output, so we'll always run with
    ## the same set of parameters.  Instead, we'll make new likelihood functions
    ## with the comparison data changed in specific, controlled ways, and we'll
    ## verify that the effect on the likelihood output is what we anticipate.

    doParallel::registerDoParallel()

    ## 1. Mean calibration, perfect model match.
    inifiles <- c(rcp85=ini85, rcp45=ini45)
    llfun1 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata, 0.1,
                                 'maxb','mina',NULL, 2100, 'rcp85', 0.2)
    out <- expect_silent(llfun1(parms))
    expected <- nrow(compdata)*stats::dnorm(0,0,parms['sigt'], log=TRUE) # assumes sigt==sigco2
    expect_equal(out, expected)

    ## 2. Mean calibration, discrepancy in the rcp85 scenario
    compdata2 <- dplyr::mutate(compdata,
                               cmean=dplyr::if_else(experiment!='rcp85', cmean,
                               dplyr::if_else(variable=='tas', cmean-1, cmean+1)))
    llfun2 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata2,
                                 0.1,'maxb','mina', NULL, 2100, 'rcp85', 0.2)
    out <- expect_silent(llfun2(parms))
    k1 <- stats::dnorm(0,0,parms['sigt'], log=TRUE)
    k2 <- stats::dnorm(1,0,parms['sigt'], log=TRUE)
    exptct <- dplyr::count(compdata2, experiment)
    n <- exptct$n
    names(n) <- exptct$experiment
    expected <- as.vector(k1*(n['historical']+n['rcp85']) + k2*n['rcp85'])
    out <- expect_silent(llfun2(parms))
    expect_equal(out, expected)


    ## 3. Envelope calibration, all values in range
    llfun3 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata, 0.1,
                                 'maxb','mina',NULL, 2100, 'rcp85', 0.2)
    out <- expect_silent(llfun3(parms))
    expected <- nrow(compdata) * log(mesa(2, 0, 4, 0.4))
    expect_equal(out, expected)

    ## 4. Envelope calibration, all values on the edge of the range
    compdata4 <- dplyr::mutate(compdata,
                               mina=mina+2, maxb=maxb+2)
    llfun4 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata4, 0.1,
                                 'maxb','mina', NULL, 2100, 'rcp85', 0.2)
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

    parms <- c(2.5, 2.5, 1.0, 1.0, 1.0)
    names(parms) <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE(),
                      'sig')

    compdata <- readRDS('pc_compdata.rds')
    pcidx <- as.integer(substring(compdata$variable, 3))

    ## See notes in the previous test; we have 4 cases to run

    doParallel::registerDoParallel()
    inifiles <- c(rcp85=ini85, rcp45=ini45)
    ## 1. Mean calibration, perfect model match
    llfun1 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata, 0.1, 'maxb',
                                 'mina', pcs, 2100, 'rcp85', 0.2)
    out <- expect_silent(llfun1(parms))
    expected <- nrow(compdata)*stats::dnorm(0,0,parms['sig'], log=TRUE) # assumes sigt==sigco2
    expect_equal(out, expected)

    ## 2. Mean calibration, imperfect model match
    compdata2 <- dplyr::mutate(compdata, cmean = dplyr::if_else(pcidx%%2==0, cmean+1,
                                         cmean-1))
    llfun2 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata2, 0.1,
                                 'maxb','mina',pcs, 2100, 'rcp85', 0.2)
    expected <- nrow(compdata)*stats::dnorm(1, 0, parms['sig'], log=TRUE)
    out <- expect_silent(llfun2(parms))
    expect_equal(out, expected)

    ## 3. Envelope callibration, centered
    llfun3 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata, 0.1, 'maxb',
                                 'mina', pcs, 2100, 'rcp85', 0.2)
    expected <- nrow(compdata) * log(mesa(2, 0, 4, 0.4))
    out <- expect_silent(llfun3(parms))
    expect_equal(out, expected)

    ## 4. Envelope calibration, edge
    compdata4 <- dplyr::mutate(compdata, mina=dplyr::if_else(pcidx%%2==0, mina+2,
                                         mina-2),
                               maxb=dplyr::if_else(pcidx%%2==0, maxb+2, maxb-2))
    llfun4 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata4, 0.1, 'maxb',
                                 'mina', pcs, 2100, 'rcp85', 0.2)
    expected <- nrow(compdata) * log(mesa(0,0,4,0.4))
    out <- expect_silent(llfun4(parms))
    expect_equal(out, expected)

    doParallel::stopImplicitCluster()
})


test_that('Posterior functions are assembled correctly from priors and posteriors', {
    doParallel::registerDoParallel()
    ## test with output comparison data
    compdata <- readRDS('out_compdata.rds')
    ini45 <- system.file('input/hector_rcp45.ini', package='hector')
    ini85 <- system.file('input/hector_rcp85.ini', package='hector')
    inifiles <- c(rcp45=ini45, rcp85=ini85)

    ## test parameters for build_mcmc_post
    test_prior_params <- c(
        ecsmu=2.0, ecssig=1.0,
        kappamu=2.5,
        betamu=0.5, betasig=1.0,
        sigtscale=1.5)

    ## make_logprior needs a full set of parameters (including ones that will be
    ## left at their default values)
    full_prior_params <- c(
        ecsmu=2.0, ecssig=1.0,
        aeromu=1.0, aerosig=1.4,
        volmu=1.0, volsig=1.4,
        kappamu=2.5, kappasig=2.0,
        betamu=0.5, betasig=1.0,
        q10mu=2.0, q10sig=2.0,
        c0mu=285, c0sig=14.0,
        sigtscale=1.5,
        sigco2scale=10.0,
        sigscale=1.0)

    ## output, mean cal
    lpostfunc1 <- build_mcmc_post(compdata, inifiles,
                                  prior_params=test_prior_params)
    lp1 <- make_logprior(full_prior_params, TRUE)
    ll1 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata, 0.1, 'maxb',
                              'mina', NULL, 2100, 'rcp85', 0.2)
    lpostfunc1a <- function(p) {lp1(p) + ll1(p)}

    parms <- c(2.5, 2.5, 1.0, 1.0, 0.5, 2.0, 280.0, 1.0, 1.0)
    names(parms) <- c(hector::ECS(), hector::DIFFUSIVITY(),
                      hector::AERO_SCALE(), hector::VOLCANIC_SCALE(),
                      hector::BETA(), hector::Q10_RH(),
                      hector::PREINDUSTRIAL_CO2(),
                      'sigt','sigco2')
    expect_equal(lpostfunc1(parms), lpostfunc1a(parms))

    ## output, envelope cal
    lpostfunc2 <- build_mcmc_post(compdata, inifiles,
                                  cal_mean=FALSE,
                                  use_lnorm_ecs=FALSE,
                                  prior_params=test_prior_params)
    lp2 <- make_logprior(full_prior_params, FALSE)
    ll2 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata, 0.1, 'maxb',
                              'mina', NULL, 2100, 'rcp85', 0.2)
    lpostfunc2a <- function(p) {lp2(p) + ll2(p)}
    expect_equal(lpostfunc2(parms), lpostfunc2a(parms))

    ## two more tests with PCA calibration
    compdata <- readRDS('pc_compdata.rds')
    pcs <- readRDS('pc-conc-historical-rcp45-rcp85.rds')
    parms <- c(2.5, 2.5, 1.0, 1.0, 1.0)
    names(parms) <- c(hector::ECS(), hector::DIFFUSIVITY(),
                      hector::AERO_SCALE(), hector::VOLCANIC_SCALE(),
                      'sig')

    ## pca, mean cal
    lpostfunc3 <- build_mcmc_post(compdata, inifiles, pcs,
                                  cal_mean=TRUE,
                                  use_lnorm_ecs=FALSE,
                                  prior_params=test_prior_params)
    lp3 <- make_logprior(full_prior_params, FALSE)
    ll3 <- make_loglikelihood(inifiles, FALSE, TRUE, compdata, 0.1, 'maxb',
                              'mina', pcs, 2100, 'rcp85', 0.2)
    lpostfunc3a <- function(p) {lp3(p) + ll3(p)}
    expect_equal(lpostfunc3(parms), lpostfunc3a(parms))


    ## home stretch: pca, envelope cal
    lpostfunc4 <- build_mcmc_post(compdata, inifiles, pcs,
                                  cal_mean=FALSE,
                                  use_lnorm_ecs=TRUE,
                                  prior_params=test_prior_params)
    lp4 <- make_logprior(full_prior_params, TRUE)
    ll4 <- make_loglikelihood(inifiles, FALSE, FALSE, compdata, 0.1, 'maxb',
                              'mina', pcs, 2100, 'rcp85', 0.2)
    lpostfunc4a <- function(p) {lp4(p) + ll4(p)}
    expect_equal(lpostfunc4(parms), lpostfunc4a(parms))

    doParallel::stopImplicitCluster()
})



