library('dplyr')
library('hector')
library('hectorcal')

make_pc_compdata <- function()
{
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
    hectorcal:::parameterize_core(parms[1:nhectorparm], hcore45)
    run(hcore45,2100)

    dh <- fetchvars(hcore45, histyears, GLOBAL_TEMP(),
                    'historical')
    d45 <- fetchvars(hcore45, futyears, GLOBAL_TEMP())
    shutdown(hcore45)

    hcore85 <- newcore(ini85, name='rcp85')
    hectorcal:::parameterize_core(parms[1:nhectorparm], hcore85)
    run(hcore85, 2100)
    d85 <- fetchvars(hcore85, futyears, GLOBAL_TEMP())
    shutdown(hcore85)

    scendata <- bind_rows(dh,d45,d85) %>% rename(experiment=scenario) %>%
      mutate(variable=hvar2esmvar(variable))
    pcproj <- project_climate(scendata, pcs, FALSE)
    pcproj <- pcproj[1:npc]

    data.frame(PC=seq_along(pcproj), cmean=pcproj, mina=pcproj-2,
               maxb=pcproj+2)
}


make_output_compdata <- function()
{
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
    hectorcal:::parameterize_core(parms[1:nhectorparm], hcore45)
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
    hectorcal:::parameterize_core(parms[1:nhectorparm], hcore85)
    run(hcore85, 2010)
    d85 <- fetchvars(hcore85, futyears, c(GLOBAL_TEMP(), ATMOSPHERIC_CO2()))
    shutdown(hcore85)

    bind_rows(dh,d45,d85) %>%
      rename(experiment=scenario, cmean=value) %>%
      mutate(cmedian=cmean, mina=cmean-2, maxb=cmean+2, a10=cmean-1,
             b90=cmean+1) %>%
      mutate(variable=if_else(variable=='Tgav', 'tas', 'co2')) %>%
      select(year, variable, experiment, mina, maxb, a10, b90, cmean,
             cmedian) %>%
      arrange(year) # This last one is to test that scrambled data gets sorted out.
}

pc_compdata <- make_pc_compdata()
saveRDS(pc_compdata, '../../tests/testthat/pc_compdata.rds', compress='xz')

out_compdata <- make_output_compdata()
saveRDS(out_compdata, '../../tests/testthat/out_compdata.rds', compress='xz')
