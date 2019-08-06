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

test_that('cononicalize_expt_names converts names.', {
    expt_names <- c('rcp85', 'esmrcp85', 'esmHistorical', 'RCP85')
    cnames <- canonicalize_expt_names(expt_names)
    expect_equal(cnames, c('rcp85', 'rcp85', 'historical', 'rcp85'))
})
