context('Single ESM calibration functions')

# Hector ini files
ini_f1    <- system.file('input/hector_rcp60.ini', package = 'hector')
ini_f2    <- system.file('input/hector_rcp45.ini', package = 'hector')
input_ini <- c(ini_f1, ini_f2)

# Make a list of Hector cores to use in the tests.
out <- setup_hector_cores(inifile = input_ini, name = c('rcp60', 'rcp45'))

test_that('setup_cores works', {

    # We expect that the length of the hector core list returned is equal to the
    # length of the ini files used as setup_hector_cores intput.
    testthat::expect_equal(length(out), length(input_ini))

    # Expect an error to be thrown when there is not a name for every hector core.
    testthat::expect_error(setup_hector_cores(inifile = input_ini, name = c('rcp60')),
                           'inifile and name must be vectors of equal lengths')


})

test_that('parameterize_core throws expected errors', {

    # If params is an unamed vector.
    testthat::expect_error(parameterize_core(params = 1:3, core = out[[1]]), 'params must be named.')

    # If params contains a bad hector parameter.
    params <- 1:3
    names(params) <- c(hector::ECS(), 'fake', hector::AERO_SCALE())
    testthat::expect_warning(testthat::expect_error(parameterize_core(params = params, core = out[[1]]), 'Bogus parameter names'), 'Units entries for the following variables not found: fake')

})

test_that('make sure that parameterize_core works', {

    # Check to see that the temperature parameters are reset correctly.
    temp_param        <- seq(from = 5, by = 5, length.out = 4)
    names(temp_param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY(), hector::VOLCANIC_SCALE())
    lapply(X = out, FUN = parameterize_core, params = temp_param)

    new_cs   <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::ECS())[['value']]})
    new_aero <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::AERO_SCALE())[['value']]})
    new_dif  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::DIFFUSIVITY())[['value']]})
    new_vol  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::VOLCANIC_SCALE())[['value']]})

    testthat::expect(all(new_cs == temp_param[[hector::ECS()]]))
    testthat::expect(all(new_aero == temp_param[[hector::AERO_SCALE()]]))
    testthat::expect(all(new_dif == temp_param[[hector::DIFFUSIVITY()]]))
    testthat::expect(all(new_vol == temp_param[[hector::VOLCANIC_SCALE()]]))

    # Check to make sure that the climate and the carbon parameters are reest correctly.
    temp_carbon_param        <- seq(from = 25, by = 5, length.out = 7)
    names(temp_carbon_param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY(), hector::VOLCANIC_SCALE(), hector::BETA(),
                                  hector::Q10_RH(), hector::PREINDUSTRIAL_CO2())


    lapply(X = out, FUN = parameterize_core, params = temp_carbon_param)

    new_ecs  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::ECS())[['value']]})
    new_aero <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::AERO_SCALE())[['value']]})
    new_dif  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::DIFFUSIVITY())[['value']]})
    new_beta <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::BETA())[['value']]})
    new_q10  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::Q10_RH())[['value']]})
    new_pre  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::PREINDUSTRIAL_CO2())[['value']]})
    new_vol  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::VOLCANIC_SCALE())[['value']]})

    testthat::expect(all(new_ecs == temp_carbon_param[[hector::ECS()]]))
    testthat::expect(all(new_aero == temp_carbon_param[[hector::AERO_SCALE()]]))
    testthat::expect(all(new_dif == temp_carbon_param[[hector::DIFFUSIVITY()]]))
    testthat::expect(all(new_beta == temp_carbon_param[[hector::BETA()]]))
    testthat::expect(all(new_q10 == temp_carbon_param[[hector::Q10_RH()]]))
    testthat::expect(all(new_pre == temp_carbon_param[[hector::PREINDUSTRIAL_CO2()]]))
    testthat::expect(all(new_vol == temp_carbon_param[[hector::VOLCANIC_SCALE()]]))


})

test_that('make_minimize_function throws errors', {

    # Subset the individual cmip data set to use in the minimize functions.
    esm_data <- dplyr::filter(cmip_individual, experiment %in% c('rcp45', 'rcp60') & variable == 'tas')
    pc       <- pc_conc

    # Set up the param vector.
    param <- c(3, 1, 2)
    names(param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY())

    # Should only work with one model.
    testthat::expect_error(make_minimize_function(hector_cores = out, esm_data = cmip_individual, normalize =  pc, param = param),
                           'esm_data can only have input for a single ESM')

    # The Hector cores are set up without all the expected experiments.
    cmip_single_model <- dplyr::filter(cmip_individual, model == 'ACCESS1-0')
    testthat::expect_error(make_minimize_function(hector_cores = out, esm_data = cmip_single_model, normalize  = pc, param = param),
                           'hector_cores are missing cores for the following esm experiments: historical, rcp85')

    # The center and scale values used to normalize the esm and hector data is missing a year.
    esm_data <- dplyr::filter(cmip_individual, model == 'CESM1-CAM5' & experiment %in% c('rcp45', 'rcp60') &
                                  ensemble == 'r1i1p1' & variable == 'tas' & year <= 2102) %>%
        dplyr::bind_rows(tibble::tibble(model ='CESM1-CAM5',
                                        experiment = 'rcp45',
                                        ensemble = 'r1i1p1',
                                        variable = 'tas',
                                        year = 2102))
    testthat::expect_error(make_minimize_function(hector_cores = out, esm_data = esm_data, normalize = pc, param = param),
                           'Missing scale and center values for :rcp45.tas.2102')

    # Make sure that if center and scale are missing names that an error message is thrown.
    missing_scale_names <- pc
    names(missing_scale_names$scale) <- NA
    testthat::expect_error(make_minimize_function(hector_cores = out, esm_data = esm_data, normalize = missing_scale_names, param = param),
                           'normalize scale needs names')

    missing_center_names <- pc
    names(missing_center_names$center) <- NA
    testthat::expect_error(make_minimize_function(hector_cores = out, esm_data = esm_data, normalize = missing_center_names, param = param),
                           'normalize center needs names')

})

test_that('make_minimize_function works with clim parameters', {

    # Hector ini files
    ini_f1    <- system.file('input/hector_rcp60_constrained.ini', package = 'hector')
    ini_f2    <- system.file('input/hector_rcp45_constrained.ini', package = 'hector')
    input_ini <- c(ini_f1, ini_f2)
    new_cores <- setup_hector_cores(inifile = input_ini, name = c('rcp60', 'rcp45'))

    # Define Hector parameters.
    param <- c(1, 1, 1)
    names(param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::VOLCANIC_SCALE())

    # Run Hector with the parameters to generate comparison data.
    lapply(new_cores, function(core){

        # Parameterize the core
        x <- parameterize_core(core = core, params = param)

        # Run Hector
        x <- hector::run(core)

        # Extract the Hector results
        hector::fetchvars(core = core, vars = hector::GLOBAL_TEMP(), dates = 2006:2100) %>%
            dplyr::rename('experiment' = scenario)


    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(model = 'selfTest',
                      variable = 'tas') ->
        comp_data

    # Expect a MSE value of 0 when Hector output data is compared with itself.
    fn <- make_minimize_function(hector_cores = new_cores, esm_data = comp_data, normalize = pc_conc, param, n = 1)
    testthat::expect_equal(fn(param), 0)

    # Make sure that the make minimize function works when parallelized.
    fn <- make_minimize_function(new_cores, comp_data, normalize = pc_conc, param, n = 2)
    testthat::expect_equal(fn(param), 0)

    # Make sure that the minimize function returns a different answer when different parameters are used.
    param2 <- c(2, 1, 1)
    names(param2) <- c(hector::ECS(), hector::AERO_SCALE(), hector::VOLCANIC_SCALE())
    testthat::expect_true(fn(param2) != 0)

    # What happens when the comparison data shifts by a values of 1? If the average
    # difference between the comparison data for each experiment is 1 then we expect that
    # sum of the MSE for the two experiments will equal 2.
    #
    # Shift the comparison data by a set numer.
    shift_by <- 1
    comp_data_shifted <- comp_data
    comp_data_shifted$value <- comp_data_shifted$value - shift_by

    # Change the center values to 0 and the scale values to 1 so that the noramlizing
    # process does not will not change the hector output or comparison data values.
    scale  <- pc_conc$scale / pc_conc$scale
    center <- pc_conc$center * 0
    names(scale) <- names(pc_conc$scale)
    names(center) <- names(pc_conc$center)
    norm = list("scale" = scale, "center" = center)

    fn_shifted    <- make_minimize_function(hector_cores = new_cores, esm_data = comp_data_shifted,
                                            normalize = norm, param, n = 1)
    testthat::expect_equal(fn_shifted(param), shift_by)

})

test_that('make_minimize_function works with co2 parameters', {

    # Hector ini files
    ini_f1    <- system.file('input/hector_rcp60.ini', package = 'hector')
    ini_f2    <- system.file('input/hector_rcp85.ini', package = 'hector')
    input_ini <- c(ini_f1, ini_f2)
    new_cores <- setup_hector_cores(inifile = input_ini, name = c('esmHistorical', 'esmrcp85'))

    # Define Hector parameters.
    param <- c(1, 200, 0.5, 0.5)
    names(param) <- c(hector::BETA(), hector::PREINDUSTRIAL_CO2(), hector::Q10_RH(), hector::ECS())

    # Run Hector with the parameters to generate comparison data.
    lapply(new_cores, function(core){

        # Parameterize the core
        x <- parameterize_core(core = core, params = param)

        # Run Hector
        x <- hector::run(core)

        # Extract the Hector results
        hector::fetchvars(core = core, vars = c(hector::GLOBAL_TEMP(), hector::ATMOSPHERIC_CO2()), dates = 1861:2100) %>%
            dplyr::rename('experiment' = scenario)


    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(model = 'selfTest',
                      variable = dplyr::if_else(variable == 'Tgav', 'tas', 'co2')) ->
        hector_output

    # Subset the comparison data so that the comparison data only contains values there is center / scale
    # values to normalize the data.
    hector_output %>%
        dplyr::filter(experiment == 'esmrcp85' & year > 2006) ->
        comp_data1

    hector_output %>%
        dplyr::filter(experiment == 'esmHistorical' & year < 2006) ->
        comp_data2

    comp_data <- dplyr::bind_rows(comp_data1, comp_data2)


    # Expect a MSE value of 0 when Hector output data is compared with itself.
    fn <- make_minimize_function(hector_cores = new_cores, esm_data = comp_data, normalize = pc_emiss, param, n = 1)
    testthat::expect_equal(fn(param), 0)

})

test_that('make_minimize_function works with heatflux', {

    # Hector ini files
    ini_f1    <- system.file('input/hector_rcp60_constrained.ini', package = 'hector')
    ini_f2    <- system.file('input/hector_rcp45_constrained.ini', package = 'hector')
    input_ini <- c(ini_f1, ini_f2)
    new_cores <- setup_hector_cores(inifile = input_ini, name = c('rcp60_r1i1p1', 'rcp45_r1i1p1'))

    # Define Hector parameters.
    param <- c(1, 1, 1)
    names(param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::VOLCANIC_SCALE())

    # Run Hector with the parameters to generate comparison data.
    lapply(new_cores, function(core){

        # Parameterize the core
        x <- parameterize_core(core = core, params = param)

        # Run Hector
        x <- hector::run(core)

        # Extract the Hector results
        hector::fetchvars(core = core, vars = c(hector::GLOBAL_TEMP(), hector::HEAT_FLUX()), dates = 2006:2100) %>%
            dplyr::rename('experiment' = scenario)


    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(model = 'selfTest',
                      variable = dplyr::if_else(variable == hector::GLOBAL_TEMP(), 'tas', 'heatflux')) ->
        comp_data

    # Expect a MSE value of 0 when Hector output data is compared with itself.
    fn <- make_minimize_function(hector_cores = new_cores, esm_data = comp_data, normalize = center_scale, param, n = 1)
    testthat::expect_equal(fn(param), 0)

    # Make sure that the make minimize function works when parallelized.
    fn <- make_minimize_function(new_cores, comp_data, normalize = center_scale, param, n = 2)
    testthat::expect_equal(fn(param), 0)

    # Make sure that the minimize function returns a different answer when different parameters are used.
    param2 <- c(2, 1, 1)
    names(param2) <- c(hector::ECS(), hector::AERO_SCALE(), hector::VOLCANIC_SCALE())
    testthat::expect_true(fn(param2) != 0)

    # What happens when the comparison data shifts by a values of 1? If the average
    # difference between the comparison data for each experiment is 1 then we expect that
    # sum of the MSE for the two experiments will equal 2.
    #
    # Shift the comparison data by a set numer.
    shift_by <- 1
    comp_data %>%
        dplyr::mutate(value = dplyr::if_else(variable == 'heatflux', value - shift_by, value)) ->
        comp_data_shifted

    # Change the center values to 0 and the scale values to 1 so that the noramlizing
    # process does not will not change the hector output or comparison data values.
    scale  <- center_scale$scale / center_scale$scale
    center <- center_scale$center * 0
    names(scale) <- names(center_scale$scale)
    names(center) <- names(center_scale$center)
    norm = list("scale" = scale, "center" = center)

    fn_shifted    <- make_minimize_function(hector_cores = new_cores, esm_data = comp_data_shifted,
                                            normalize = norm, param, n = 4, showMessages = TRUE)
    testthat::expect_equal(fn_shifted(param), shift_by)

})

test_that('make_minimize_function weights correctly', {

    # Hector ini files
    ini_f1    <- system.file('input/hector_rcp60_constrained.ini', package = 'hector')
    ini_f2    <- system.file('input/hector_rcp60_constrained.ini', package = 'hector')
    input_ini <- c(ini_f1, ini_f2)
    new_cores <- setup_hector_cores(inifile = input_ini, name = c('rcp60_r1i1p1', 'rcp60_r1i1p2'))

    # Define Hector parameters.
    param <- c(1, 1, 1)
    names(param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::VOLCANIC_SCALE())

    # Run Hector with the parameters to generate comparison data.
    lapply(new_cores, function(core){

        # Parameterize the core
        x <- parameterize_core(core = core, params = param)

        # Run Hector
        x <- hector::run(core)

        # Extract the Hector results
        hector::fetchvars(core = core, vars = hector::GLOBAL_TEMP(), dates = 2006:2100) %>%
            dplyr::rename('experiment' = scenario) %>%
            dplyr::mutate(variable = 'tas')


    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(model = 'selfTest') ->
        comp_data

    # Shift the comparison data for a singel ensemble by 1.
    shift_by       <- 1
    shift_ensemble <- dplyr::mutate(comp_data, value = dplyr::if_else(experiment == 'rcp60_r1i1p2', value + shift_by, value))

    # Change the center values to 0 and the scale values to 1 so that the noramlizing
    # process does not will not change the hector output or comparison data values.
    scale  <- center_scale$scale / center_scale$scale
    center <- center_scale$center * 0
    names(scale) <- names(center_scale$scale)
    names(center) <- names(center_scale$center)
    norm = list("scale" = scale, "center" = center)

    # Expect a MSE value of 0 when Hector output data is compared with itself.
    fn <- make_minimize_function(hector_cores = new_cores, esm_data = shift_ensemble, normalize = norm, param, n = 3)
    testthat::expect_equal(fn(param), shift_by/2)

})
