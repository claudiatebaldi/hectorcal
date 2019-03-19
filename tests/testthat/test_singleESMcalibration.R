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
    temp_param        <- seq(from = 5, by = 5, length.out = 3)
    names(temp_param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY())
    lapply(X = out, FUN = parameterize_core, params = temp_param)

    new_cs   <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::ECS())[['value']]})
    new_aero <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::AERO_SCALE())[['value']]})
    new_dif  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::DIFFUSIVITY())[['value']]})

    testthat::expect(all(new_cs == temp_param[[hector::ECS()]]))
    testthat::expect(all(new_aero == temp_param[[hector::AERO_SCALE()]]))
    testthat::expect(all(new_dif == temp_param[[hector::DIFFUSIVITY()]]))

    # Check to make sure that the climate and the carbon parameters are reest correctly.
    temp_carbon_param        <- seq(from = 25, by = 5, length.out = 6)
    names(temp_carbon_param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY(), hector::BETA(),
                                     hector::Q10_RH(), hector::PREINDUSTRIAL_CO2())


    lapply(X = out, FUN = parameterize_core, params = temp_carbon_param)

    new_ecs  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::ECS())[['value']]})
    new_aero <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::AERO_SCALE())[['value']]})
    new_dif  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::DIFFUSIVITY())[['value']]})
    new_beta <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::BETA())[['value']]})
    new_q10  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::Q10_RH())[['value']]})
    new_pre  <- lapply(out, function(x){hector::fetchvars(core = x, dates = NA, vars = hector::PREINDUSTRIAL_CO2())[['value']]})

    testthat::expect(all(new_ecs == temp_carbon_param[[hector::ECS()]]))
    testthat::expect(all(new_aero == temp_carbon_param[[hector::AERO_SCALE()]]))
    testthat::expect(all(new_dif == temp_carbon_param[[hector::DIFFUSIVITY()]]))
    testthat::expect(all(new_beta == temp_carbon_param[[hector::BETA()]]))
    testthat::expect(all(new_q10 == temp_carbon_param[[hector::Q10_RH()]]))
    testthat::expect(all(new_pre == temp_carbon_param[[hector::PREINDUSTRIAL_CO2()]]))

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

    # The center and scale values used to normalize the esm and hector data is missing a few years.
    esm_data <- dplyr::filter(cmip_individual, model == 'CESM1-CAM5' & experiment %in% c('rcp45', 'rcp60') &
                                  ensemble == 'r1i1p1' & variable == 'tas' & year <= 2102)
   testthat::expect_error(make_minimize_function(hector_cores = out, esm_data = esm_data, normalize = pc, param = param),
                          'Missing scale and center values for :rcp45.tas.2101, rcp45.tas.2102, rcp60.tas.2101, rcp60.tas.2102')

})

test_that('make_minimize_function works with clim parameters', {

    # Hector ini files
    ini_f1    <- system.file('input/hector_rcp60_constrained.ini', package = 'hector')
    ini_f2    <- system.file('input/hector_rcp45_constrained.ini', package = 'hector')
    input_ini <- c(ini_f1, ini_f2)
    new_cores <- setup_hector_cores(inifile = input_ini, name = c('rcp60', 'rcp45'))
    other_new <- setup_hector_cores(inifile = input_ini, name = c('rcp60', 'rcp45'))


    # Climate model data to calibrated to.
    esm_data <- dplyr::filter(cmip_individual, model == 'CESM1-CAM5' & experiment %in% c('rcp45', 'rcp60') &
                                  ensemble == 'r1i1p1' & variable == 'tas' & year <= 2100)

    # Run Hector with the default values and format the output.
    lapply(new_cores, function(input){

        hector::reset(input)
        hector::run(input)

        hector::fetchvars(input, dates = 2000:2100, vars = hector::GLOBAL_TEMP()) %>%
            dplyr::mutate(index = paste0(scenario, '.tas.', year))

    }) %>%
        dplyr::bind_rows() ->
        hector_output

    # Normalize the Hector output data.
    hector_output %>%
        dplyr::left_join(tibble::tibble(scale = pc_conc$scale,
                                        index = names( pc_conc$scale)), by = 'index') %>%
        dplyr::left_join(tibble::tibble(center =  pc_conc$center,
                                        index = names( pc_conc$center)), by = 'index') %>%
        na.omit %>%
        dplyr::rename(hector_out = value) %>%
        dplyr::mutate(hec_norm_val = (hector_out - center) / scale) %>%
        dplyr::select(index, hector_out, hec_norm_val) ->
        hector_comp_data

    # Normalize the ESM output data.
    esm_data %>%
        dplyr::mutate(index = paste0(experiment, '.', variable, '.', year)) %>%
        dplyr::left_join(tibble::tibble(scale =  pc_conc$scale,
                                        index = names( pc_conc$scale)), by = 'index') %>%
        dplyr::left_join(tibble::tibble(center =  pc_conc$center,
                                        index = names( pc_conc$center)), by = 'index') %>%
        dplyr::mutate(esm_norm_val = (value - center) / scale) %>%
        dplyr::rename(esm_out = value) %>%
        dplyr::select(index, esm_out, esm_norm_val, experiment) ->
        esm_comp_data

    # Calcualte the difference between the normalized Hector and ESM output data.
    hector_comp_data %>%
        dplyr::left_join(esm_comp_data, by = 'index') %>%
        na.omit() %>%
        dplyr::mutate(SE = (hec_norm_val - esm_norm_val)^2) %>%
        split(.$experiment) %>%
        lapply(function(input){

            mean(input$SE)

        }) %>%
        unlist %>%
        sum ->
        expected_MSE

    # Calculate the MSE using the function returned by make_minimize_function.
    # Use the default Hector S, aero, and diff parameters since we know that Hector should solve with them.
    param_default <- c(3, 1, 2.3)
    names(param_default) <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY())
    fn <- make_minimize_function(hector_cores = other_new, esm_data = esm_data,
                                            normalize = pc_conc, param = param_default)


    # Make sure that the MSE calculated by the function matches the value calcualted by hand.
    testthat::expect_equal(fn(param = param_default), expected_MSE)

    # When the CS is unrealistic we expect the function to return Inf.
    testthat::expect_equal(fn(c(-6, 1, 4)), Inf)

})

test_that('make_minimize_function works with clim and carbon parameters', {

    # Hector ini files
    ini_f1    <- system.file('input/hector_rcp85_constrained.ini', package = 'hector')
    ini_f2    <- system.file('input/hector_rcp60_constrained.ini', package = 'hector')
    input_ini <- c(ini_f1, ini_f2)
    new_cores <- setup_hector_cores(inifile = input_ini, name = c('esmrcp85', 'esmHistorical'))

    # Climate model data to calibrated to.
    esm_data <- dplyr::filter(cmip_individual, model == 'CanESM2' &
                                  ensemble == 'r1i1p1' &
                                  year <= 2100 & year >= 1861 &
                                  grepl('esm', experiment))

    # Run Hector with the default values and format the output.
    lapply(new_cores, function(input){

        hector::reset(input)
        hector::run(input)

        hector::fetchvars(input, dates = 1861:2100, vars = c(hector::GLOBAL_TEMP(), hector::ATMOSPHERIC_CO2())) %>%
            dplyr::mutate(index = dplyr::if_else(variable == 'Tgav', paste0(scenario, '.tas.', year), paste0(scenario, '.co2.', year)))

    }) %>%
        dplyr::bind_rows() ->
        hector_output

    # Normalize the Hector output data.
    hector_output %>%
        dplyr::left_join(tibble::tibble(scale = pc_emiss$scale,
                                        index = names( pc_emiss$scale)), by = 'index') %>%
        dplyr::left_join(tibble::tibble(center =  pc_emiss$center,
                                        index = names( pc_emiss$center)), by = 'index') %>%
        na.omit %>%
        dplyr::rename(hector_out = value) %>%
        dplyr::mutate(hec_norm_val = (hector_out - center) / scale) %>%
        dplyr::select(index, hector_out, hec_norm_val) ->
        hector_comp_data

    # Normalize the ESM output data.
    esm_data %>%
        dplyr::mutate(index = paste0(experiment, '.', variable, '.', year)) %>%
        dplyr::left_join(tibble::tibble(scale =  pc_emiss$scale,
                                        index = names( pc_emiss$scale)), by = 'index') %>%
        dplyr::left_join(tibble::tibble(center =  pc_emiss$center,
                                        index = names( pc_emiss$center)), by = 'index') %>%
        dplyr::mutate(esm_norm_val = (value - center) / scale) %>%
        dplyr::rename(esm_out = value) %>%
        dplyr::select(index, esm_out, esm_norm_val, experiment) ->
        esm_comp_data

    # Calcualte the difference between the normalized Hector and ESM output data.
    hector_comp_data %>%
        dplyr::left_join(esm_comp_data, by = 'index') %>%
        na.omit() %>%
        dplyr::mutate(SE = (hec_norm_val - esm_norm_val)^2) %>%
        split(.$experiment) %>%
        lapply(function(input){

            mean(input$SE)

        }) %>%
        unlist %>%
        sum ->
        expected_MSE

    # Calculate the MSE using the function returned by make_minimize_function.
    # Use the default Hector S, aero, and diff parameters since we know that Hector should solve with them.
    param_default <- c(3, 1, 2.3, 0.36, 2, 276.0897)
    names(param_default) <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY(),
                              hector::BETA(),  hector::Q10_RH(), hector::PREINDUSTRIAL_CO2())
    fn <- make_minimize_function(hector_cores = new_cores, esm_data = esm_data,
                                            normalize = pc_emiss, param = param_default)

    # Make sure that the MSE calculated by the function matches the value calcualted by hand.
    testthat::expect_equal(fn(param = param_default), expected_MSE)

    # Make sure that when the input paramters change a bit that the MSE also changes.
    testthat::expect_true(fn(param = c(3.5, 1, 2.3, 0.36, 2, 276.0897)) != expected_MSE)

    # When the CS is unrealistic we expect the function to return Inf.
    testthat::expect_equal(fn(c(-6, 1, 4, 0.36, 2, 276.0897)), Inf)

})


