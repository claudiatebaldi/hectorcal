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

test_that('make sure that make_param_penalty_function throws errors', {

    testthat::expect_error(make_param_penalty_function(func_list = c( v = 1)),
                           regexp = 'func_list must be a list with names')
    testthat::expect_error(make_param_penalty_function(func_list = list(1)),
                           regexp = 'func_list must be a list with names')
    testthat::expect_error(make_param_penalty_function(func_list = list(fake = 1)),
                           regexp = 'func_list names must be Hector parameters')
    testthat::expect_error(make_param_penalty_function(func_list = list(S = 1)),
                           regexp = 'func_list elements must be functions')

    # If the function is made but is trying to penalize a paramter that is not being optimized
    # the function should throw an error.
    fn          <- make_param_penalty_function(func_list = list('S' = function(x) sum(x)))
    optim_param <- c('beta' = 1)
    testthat::expect_error(fn(optim_param), regexp = 'trying to penalize parameters that are not being optimized')

})

test_that('make sure that make_param_penalty_function works', {

    # Run the make_param_penalty_function
    xx    <- list('S' = function(x){ x * 3})
    fn          <- make_param_penalty_function(xx)
    optim_param <- c('S' = 1)
    rslt        <- fn(optim_param)

    # Make sure that the function returns a data frame with the
    # expected output strucutre and values.
    testthat::expect_equal(nrow(rslt), 1)
    testthat::expect_equal(ncol(rslt), 2)
    testthat::expect_equal(rslt$value, xx$S(optim_param[[1]]))

    # Make sure that the function returns the correct values
    # regardless of the number / order of the paramters being optimized.
    optim_param <- c('S' = 1, 'beta' = 9)
    rslt2       <- fn(optim_param)
    testthat::expect_equal(rslt, rslt2)

    optim_param <- c('beta' = 9, 'S' = 1)
    rslt3       <- fn(optim_param)
    testthat::expect_equal(rslt, rslt3)

    # Test capability with mulitple penalty functions.
    xx <- list(S = function(x) x * 3,
               beta = function(x) x * 9)
    fn          <- make_param_penalty_function(xx)
    optim_param <- c('S' = 1, 'beta' = 1)
    rslt        <- fn(optim_param)

    testthat::expect_equal(dim(rslt), c(2, 2))
    testthat::expect_equal(rslt$value, c(xx$S(1), xx$beta(1)))

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

        # Parameterize and run the core
        x <- parameterize_core(core = core, params = param)
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
    fn <- make_minimize_function(hector_cores = new_cores,
                                 esm_data = comp_data,
                                 cmip_range = NULL,
                                 normalize = pc_conc,
                                 param = param,
                                 n = 1)
    testthat::expect_equal(fn(param), 0)


    # Expect the output to change when intermediateOut is set to TRUE.
    fn <- make_minimize_function(hector_cores = new_cores,
                                 esm_data = comp_data,
                                 cmip_range = NULL,
                                 normalize = pc_conc,
                                 param = param,
                                 n = 1,
                                 intermediateOutput = TRUE)
    rslt <- fn(param)
    testthat::expect_equal(dim(rslt), c(2, 4))
    testthat::expect_equal(rslt$value, c(0, 0))

    # Make sure that the make minimize function works when parallelized.
    fn <- make_minimize_function(hector_cores = new_cores, esm_data = comp_data, normalize = pc_conc, param = param, n = 2)
    testthat::expect_equal(fn(param), 0)

    # Make sure that the minimize function returns a different answer when different parameters are used.
    param2        <- c(2, 1, 1)
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

    # Run the make_param_penalty_function
    func_list <- list('S' = function(x) x * 8)
    penalty   <- make_param_penalty_function(func_list)
    fn_shifted_penalized <- make_minimize_function(hector_cores = new_cores,
                                                   esm_data = comp_data_shifted,
                                                   normalize = norm,
                                                   param = param,
                                                   cmip_range = NULL,
                                                   param_penalty = penalty,
                                                   n = 1)
    # Since output is compared with the shifted comparison data for two experiments that are weighted equally figure out
    # the expected peanlized weighted MSE.
    expected_penalized_MSE <- mean(c(1, 1, func_list$S(param['S'])))
    testthat::expect_equal(expected_penalized_MSE, fn_shifted_penalized(param))

})

test_that('make_minimize_function works with co2 parameters', {

    # Hector ini files
    ini_f1    <- system.file('input/hector_rcp60.ini', package = 'hector')
    ini_f2    <- system.file('input/hector_rcp85.ini', package = 'hector')
    input_ini <- c(ini_f1, ini_f2)
    new_cores <- setup_hector_cores(inifile = input_ini, name = c('esmHistorical', 'esmrcp85'))

    # Define Hector parameters.
    param        <- c(1, 200, 0.5, 0.5)
    names(param) <- c(hector::BETA(), hector::PREINDUSTRIAL_CO2(), hector::Q10_RH(), hector::ECS())

    # Run Hector with the parameters to generate comparison data.
    lapply(new_cores, function(core){

        # Parameterize and run Hector
        x <- parameterize_core(core = core, params = param)
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
    fn <- make_minimize_function(hector_cores = new_cores,
                                 esm_data = comp_data,
                                 normalize = pc_emiss,
                                 param = param,
                                 n = 1)
    testthat::expect_equal(fn(param), 0)

})

test_that('make_minimize_function works with heatflux', {

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
                                            normalize = norm, param, n = 4)
    testthat::expect_equal(fn_shifted(param), c(shift_by/2))

})

test_that('make_minimize_function weights correctly', {

    # Hector ini files
    ini_f1    <- system.file('input/hector_rcp60_constrained.ini', package = 'hector')
    ini_f2    <- system.file('input/hector_rcp60_constrained.ini', package = 'hector')
    input_ini <- c(ini_f1, ini_f2)
    new_cores <- setup_hector_cores(inifile = input_ini, name = c('rcp60', 'rcp602'))

    # Define Hector parameters.
    param <- c(1, 1, 1)
    names(param) <- c(hector::ECS(), hector::AERO_SCALE(), hector::VOLCANIC_SCALE())

    # Run Hector with the parameters to generate comparison data.
    lapply(new_cores, function(core){

        # Parameterize and run Hector
        x <- parameterize_core(core = core, params = param)
        x <- hector::run(core)

        # Extract the Hector results
        hector::fetchvars(core = core, vars = hector::GLOBAL_TEMP(), dates = 2006:2100) %>%
            dplyr::rename('experiment' = scenario) %>%
            dplyr::mutate(variable = 'tas')


    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(model = 'selfTest',
                      ensemble = '1') ->
        comp_data

    # Shift the comparison data for a singel ensemble by 1.
    shift_by       <- 1
    dplyr::mutate(comp_data, ensemble = dplyr::if_else(experiment == 'rcp602', '2', ensemble),
                  experiment = 'rcp60') %>%
        dplyr::mutate(value =  dplyr::if_else(ensemble == '2', value + shift_by, value)) ->
        shift_ensemble

    # Change the center values to 0 and the scale values to 1 so that the noramlizing
    # process does not will not change the hector output or comparison data values.
    scale  <- center_scale$scale / center_scale$scale
    center <- center_scale$center * 0
    names(scale)  <- names(center_scale$scale)
    names(center) <- names(center_scale$center)
    norm = list("scale" = scale, "center" = center)

    # Expect a MSE value of 0 when Hector output data is compared with itself.
    fn <- make_minimize_function(hector_cores = new_cores, esm_data = shift_ensemble, normalize = norm, param, n = 3)
    testthat::expect_equal(fn(param), shift_by/2)

})

test_that('make_minimize_function uses the cmip_range correctly', {

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
        x <- hector::run(core)

        # Extract the Hector results
        hector::fetchvars(core = core, vars = c(hector::GLOBAL_TEMP(), hector::HEAT_FLUX()), dates = 2006:2100) %>%
            dplyr::rename('experiment' = scenario)


    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(model = 'selfTest',
                      variable = dplyr::if_else(variable == hector::GLOBAL_TEMP(), 'tas', 'heatflux')) ->
        comp_data

    # Make sure that make minimize throws the expected error messages.
    testthat::expect_error(make_minimize_function(hector_cores = new_cores, esm_data = comp_data, cmip_range = esm_comparison,
                                                  normalize = center_scale, param = param, n = 1),
                           'Missing columns: lower, upper, sig')

    cmip_range <- dplyr::rename(esm_comparison, lower = mina, upper = maxb, sig = b90)
    testthat::expect_error(make_minimize_function(hector_cores = new_cores, esm_data = comp_data, cmip_range = cmip_range,
                                                  normalize = center_scale, param = param, n = 1),
                           'cmip_range cannot contain data for a variable that is also in esm_data')

    # Make sure that it works as expected.
    temp_data  <- dplyr::filter(comp_data, variable == 'tas')
    range_data <- dplyr::filter(cmip_range, variable == 'heatflux' & experiment %in% c('rcp45', 'rcp60'))

    shift_by <- 1
    comp_data %>%
        dplyr::filter(variable == 'heatflux') %>%
        dplyr::mutate(lower = value - shift_by,
                      upper = value + shift_by,
                      sig = 2) %>%
        dplyr::select(-value) ->
        range_data

    fn <- make_minimize_function(hector_cores = new_cores,
                                 esm_data = temp_data,
                                 cmip_range = range_data,
                                 normalize = center_scale,
                                 param = param,
                                 n = 1)

    expected <-  0.825
    testthat::expect_equivalent(object = fn(param), expected = expected, tolerance = 0.002)


    # Shift the comparison data for a single ensemble by 1.
    shift_by       <- 1
    shift_esm_data <- dplyr::mutate(comp_data, value = value - shift_by) %>%
        dplyr::filter(variable == 'tas')

})

test_that('generate_inital_guess works', {

    # Test to make sure that generate_inital_guess runs.
    comparison_data <- head(hector_conc_ensemble$historical)
    S_only <- generate_inital_guess(comparison_data, 'S')
    testthat::expect_equal(length(S_only), 1)
    S_diff <- generate_inital_guess(comparison_data, c('S', 'diff'))
    testthat::expect_equal(length(S_diff), 2)


    # Test the generate_inital_guess works, returns the expected results
    # for the emission driven and the concentration driven ensmble.

    # Select a run from the emission driven ensemble and use it as comparison data. If the
    # generate_inital_guess functions works then it should return parameter values for the
    # same run.
    find_params <- c(hector::ECS(), hector::DIFFUSIVITY(), hector::AERO_SCALE(), hector::PREINDUSTRIAL_CO2(), hector::Q10_RH())

    # Make sure that generate_inital_guess returns the parameter values for the emisison driven runid 2.
    # Runid 2 was selected arbitrarily.
    this_id <- 2
    expected_param_values <- hector_emiss_ensemble$params[which(hector_emiss_ensemble$params$runid == this_id), ]
    expected_param_values <- expected_param_values[names(expected_param_values) %in% find_params]

    dplyr::bind_rows(hector_emiss_ensemble$esmrcp85, hector_emiss_ensemble$esmHistorical) %>%
        filter(runid == this_id) %>%
        select(value, variable, year, experiment) ->
        emiss_comparison_data

    best_guess <- generate_inital_guess(emiss_comparison_data, find_params)
    testthat::expect_true(setequal(best_guess, expected_param_values))


    # Repeate the test with the concentraiton driven experiments.
    find_params <- c(hector::ECS(), hector::DIFFUSIVITY(), hector::AERO_SCALE())

    # Make sure that generate_inital_guess returns the parameter values for the emisison driven runid 2.
    # Runid 2 was selected arbitrarily.
    this_id <- 2
    expected_param_values <- hector_conc_ensemble$params[which(hector_conc_ensemble$params$runid == this_id), ]
    expected_param_values <- expected_param_values[names(expected_param_values) %in% find_params]

    dplyr::bind_rows(hector_conc_ensemble$rcp26, hector_conc_ensemble$rcp85) %>%
        filter(runid == this_id) %>%
        select(value, variable, year, experiment) ->
        conc_comparison_data

    best_guess <- generate_inital_guess(conc_comparison_data, find_params)
    testthat::expect_true(setequal(best_guess, expected_param_values))

})

test_that('generate_inital_guess works', {

    # Test to make sure that generate_inital_guess runs.
    comparison_data <- head(hector_conc_ensemble$historical)
    S_only <- generate_inital_guess(comparison_data, 'S')
    testthat::expect_equal(length(S_only), 1)
    S_diff <- generate_inital_guess(comparison_data, c('S', 'diff'))
    testthat::expect_equal(length(S_diff), 2)


    # Test the generate_inital_guess works, returns the expected results
    # for the emission driven and the concentration driven ensmble.

    # Select a run from the emission driven ensemble and use it as comparison data. If the
    # generate_inital_guess functions works then it should return parameter values for the
    # same run.
    find_params <- c(hector::ECS(), hector::DIFFUSIVITY(), hector::AERO_SCALE(), hector::PREINDUSTRIAL_CO2(), hector::Q10_RH())

    # Make sure that generate_inital_guess returns the parameter values for the emisison driven runid 2.
    # Runid 2 was selected arbitrarily.
    this_id <- 2
    expected_param_values <- hector_emiss_ensemble$params[which(hector_emiss_ensemble$params$runid == this_id), ]
    expected_param_values <- expected_param_values[names(expected_param_values) %in% find_params]

    dplyr::bind_rows(hector_emiss_ensemble$esmrcp85, hector_emiss_ensemble$esmHistorical) %>%
        filter(runid == this_id) %>%
        select(value, variable, year, experiment) ->
        emiss_comparison_data

    best_guess <- generate_inital_guess(emiss_comparison_data, find_params)
    testthat::expect_true(setequal(best_guess, expected_param_values))


    # Repeate the test with the concentraiton driven experiments.
    find_params <- c(hector::ECS(), hector::DIFFUSIVITY(), hector::AERO_SCALE())

    # Make sure that generate_inital_guess returns the parameter values for the emisison driven runid 2.
    # Runid 2 was selected arbitrarily.
    this_id <- 2
    expected_param_values <- hector_conc_ensemble$params[which(hector_conc_ensemble$params$runid == this_id), ]
    expected_param_values <- expected_param_values[names(expected_param_values) %in% find_params]

    dplyr::bind_rows(hector_conc_ensemble$rcp26, hector_conc_ensemble$rcp85) %>%
        filter(runid == this_id) %>%
        select(value, variable, year, experiment) ->
        conc_comparison_data

    best_guess <- generate_inital_guess(conc_comparison_data, find_params)
    testthat::expect_true(setequal(best_guess, expected_param_values))

})

