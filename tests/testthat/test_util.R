context('Utility functions')

test_that('Run configurations are encoded correctly.', {
    flags <- expand.grid(pcsflag=c(TRUE,FALSE),
                         hfflag=c(TRUE,FALSE),
                         meanflag=c(TRUE,FALSE),
                         debugflag=c(TRUE,FALSE))
    serials <- seq(0, 15)

    flagids <- expand.grid(pcsx=c(16,0),
                           hfx=c(32,0),
                           meanx=c(64,0),
                           dbgx=c(128,0))

    ans <- serials + flagids$pcsx + flagids$hfx + flagids$meanx + flagids$dbgx

    runids <- encode_runid(serials, flags$pcsflag, flags$hfflag, flags$meanflag,
                           flags$debugflag)
    expect_equal(runids, ans)

    expect_error(encode_runid(seq(1,16), flags$pcsflag, flags$hfflag, flags$meanflag,
                              flags$debugflag))
    expect_error(encode_runid(serials, seq(0,1), flags$hfflag, flags$meanflag,
                              flags$debugflag))

})


test_that('Run configurations are decoded correctly.', {
    runids <- seq(1, 127, 16)
    rslt <- decode_runid(runids)

    anssn <- rep(1, 8)
    anspcs <- rep(c(FALSE,TRUE), 4)
    anshf <- rep(c(FALSE, FALSE, TRUE, TRUE), 2)
    ansmean <- rep(c(FALSE, TRUE), c(4,4))
    ansdbg <- rep(FALSE, 8)

    expect_equal(rslt$serialnumber, anssn)
    expect_equal(rslt$pcsflag, anspcs)
    expect_equal(rslt$hfflag, anshf)
    expect_equal(rslt$meanflag, ansmean)
    expect_equal(rslt$debugflag, ansdbg)

    runids <- 130+runids
    rslt <- decode_runid(runids)
    anssn <- rep(3, 8)
    ansdbg <- rep(TRUE, 8)

    expect_equal(rslt$serialnumber, anssn)
    expect_equal(rslt$pcsflag, anspcs)
    expect_equal(rslt$hfflag, anshf)
    expect_equal(rslt$meanflag, ansmean)
    expect_equal(rslt$debugflag, ansdbg)
})

test_that('canonicalize_expt_names converts names.', {
    expt_names <- c('rcp85', 'esmrcp85', 'esmHistorical', 'RCP85')
    cnames <- canonicalize_expt_names(expt_names)
    expect_equal(cnames, c('rcp85', 'rcp85', 'historical', 'rcp85'))
})

test_that('find_common_years works correctly.', {
    ## Make up some test data frames.  We don't need actual values
    mae1vtas_h <- data.frame(year=1900:1909, model='a', ensemble='r1', variable='tas',
                             experiment='esmHistorical', stringsAsFactors = FALSE)
    mae1vco2_h <- data.frame(year=1870:1879, model='a', ensemble='r1', variable='co2',
                             experiment='esmHistorical', stringsAsFactors = FALSE)
    ## model b is identical
    mbe1vtas_h <- mae1vtas_h; mbe1vtas_h$model <- 'b'
    mbe1vco2_h <- mae1vco2_h; mbe1vco2_h$model <- 'b'

    ## model a ensemble 2 has a subset of the years
    mae2vtas_h <- mae1vtas_h[1:5,] ; mae2vtas_h$year <- 1901:1905; mae2vtas_h$ensemble='r2'
    mae2vco2_h <- mae1vco2_h[1:5,] ; mae2vco2_h$year <- 1871:1875; mae2vco2_h$ensemble='r2'

    ## model b ensemble 2 has a partially overlapping set of years for temp, non-overlapping for co2
    mbe2vtas_h <- mbe1vtas_h; mbe2vtas_h$year <- 1906:1915; mbe2vtas_h$ensemble='r2'
    mbe2vco2_h <- mbe1vco2_h; mbe2vco2_h$year <- 1900:1909; mbe2vco2_h$ensemble='r2'

    ## Also create another experiment.  This should be ignored
    mabogo <- dplyr::bind_rows(
        data.frame(year=1800:2000, model='a', ensemble='r1', variable='tas',
                   experiment='bogoscen', stringsAsFactors = FALSE),
        data.frame(year=2000:2100, model='a', ensemble='r1', variable='co2',
                   experiment='bogoscen', stringsAsFactors = FALSE))

    ## Now for the tests.
    ## 1. completely harmonized dates
    testdf <- dplyr::bind_rows(mae1vtas_h, mae1vco2_h, mbe1vtas_h, mbe1vco2_h, mabogo)
    expect_equal(find_common_years(testdf, 'esmHistorical'),
                 list(co2=1870:1879, tas=1900:1909))
    expect_equal(find_common_years(testdf, 'bogoscen'),
                 list(co2=2000:2100, tas=1800:2000))

    ## 2. Add model a ensemble 2
    testdf2 <- dplyr::bind_rows(testdf, mae2vtas_h, mae2vco2_h)
    expect_equal(find_common_years(testdf2, 'esmHistorical'),
                 list(co2=1871:1875, tas=1901:1905))

    ## 3. Add model b ensemble 2 to the base frame
    testdf3 <- dplyr::bind_rows(testdf, mbe2vtas_h, mbe2vco2_h)
    expect_equal(find_common_years(testdf3, 'esmHistorical'),
                 list(co2=integer(0), tas=1906:1909))

    ## 4. All together now
    testdf4 <- dplyr::bind_rows(testdf3, mae2vtas_h, mae2vco2_h)
    expect_equal(find_common_years(testdf4, 'esmHistorical'),
                 list(co2=integer(0), tas=integer(0)))
})
