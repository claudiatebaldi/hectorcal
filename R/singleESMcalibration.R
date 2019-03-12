#### Functions used in calibrating Hector to emulate an individual ESM.

#' Make a list of Hector cores
#'
#' The list of Hector cores will be passed into the \code{make_fn_minimize}
#'
#' @param inifile A vector of ini file paths that will be used to set up the Hector cores
#' @param name A vector of the Hector core scenario names, it should reflect CMIP5 experiment names.
#' @return A list of hector instances
#' @importFrom assertthat assert_that
#' @importFrom hector newcore
#' @export

setup_hector_cores <- function(inifile, name) {

    # Check inputs
    assertthat::assert_that(length(inifile) == length(name),
                            msg = 'inifile and name must be vectors of equal lengths')

    # Make a hector core of every ini file entry
    lapply(1:(length(inifile)), function(i) {
        hector::newcore(inifile[i], name = name[i], suppresslogging=TRUE)
    })

    }


#' Make a function that will parameterize a list of Hector cores with new Hector parameters.
#'
#' The function
#'
#' @param param a vector of named Hector parameters to used to set up the cores.
#' @importFrom hector ECS AERO_SCALE DIFFUSIVITY BETA Q10_RH PREINDUSTRIAL_CO2 reset setvar
#' @importFrom assertthat assert_that
#' @return A function that can be used to reset a sinlge Hector core with new Hector parameters.
make_parameterize_cores <- function(param){

    # Right now this function is only set up to use the three climate parameters
    # or the six carbon cycle and climate parameters. Make sure that the
    # param vector contains the correct number of parameter values.
    assertthat::assert_that(length(param) == 3 | length(param) == 6, msg = 'param must be a vector of 3 or 6 parameters.')
    assertthat::assert_that(!is.null(names(param)), msg = 'the param vector must be named.')

    # Check to make sure the param vector has the correct names and that the paramters are set up in
    # the correct indexed order.
    if(length(param) == 3){

        expected_params <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY())
        message         <- paste0('param vector elements must be in the following order ', paste0(expected_params, collapse = ', '))
        assertthat::assert_that(!FALSE %in% c(expected_params == names(param)), msg = message)

        # Return a function that will reset the climate system parameters.
        function(core, param){

            # Reset each of the Hector cores with the new paramter values.
            hector::setvar(core = core, dates = NA, var = hector::ECS(), values = param[[1]], unit = 'degC')
            hector::setvar(core = core, dates = NA, var = hector::AERO_SCALE(), values = param[[2]], unit = '(unitless)')
            hector::setvar(core = core, dates = NA, var = hector::DIFFUSIVITY(), values = param[[3]], unit = 'cm2/s')
            hector::reset(core = core)

            }


        }  else {

            # Make sure that the param vector is set up correctly.
            expected_params <- c(hector::ECS(), hector::AERO_SCALE(), hector::DIFFUSIVITY(), hector::BETA(), hector::Q10_RH(), hector::PREINDUSTRIAL_CO2())
            message         <- paste0('param vector elements must be in the following order ', paste0(expected_params, collapse = ', '))
            assertthat::assert_that(!FALSE %in% c(expected_params == names(param)), msg = message)

            # Return a function that will reset the climate and carbon cycle parameters.
            function(core, param){

                    hector::setvar(core = core, dates = NA, var = hector::ECS(), values = param[[1]], unit = 'degC')
                    hector::setvar(core = core, dates = NA, var = hector::AERO_SCALE(), values = param[[2]], unit = '(unitless)')
                    hector::setvar(core = core, dates = NA, var = hector::DIFFUSIVITY(), values = param[[3]], unit = 'cm2/s')
                    hector::setvar(core = core, dates = NA, var = hector::BETA(), values = param[[4]], unit = '(unitless)')
                    hector::setvar(core = core, dates = NA, var = hector::Q10_RH(), values = param[[5]], unit = '(unitless)')
                    hector::setvar(core = core, dates = NA, var = hector::PREINDUSTRIAL_CO2(), values = param[[6]], unit = 'ppmv CO2')
                    hector::reset(core = core)

            }

        }

    }



#' Make the function to minimize
#'
#' Create the function that will calculate the mean squared error between Hector and ESM CMIP5 output data. THis
#' function will be minimized using optim.
#'
#' @param hector_cores A list of Hector cores that will be used.
#' @param esm_data A data frame of ESM data for a single model that contains the following columns, year, model, variable, experiment.
#' @param normalize A list of center and the scale values to use to noramlize the Hector and ESM output data.
#' @param param A vector of parameter values.
#' @return A function that will calculate the mean squared error between a Hector run and ESM data.
#' @importFrom foreach %do%
#' @importFrom dplyr %>%
#' @export

make_minimize_function <- function(hector_cores, esm_data, normalize, param){

    # Check inputs
    check_columns(input = normalize, req_cols = c('scale', 'center'))
    check_columns(input = esm_data, req_cols = c('model', 'experiment', 'variable', 'year'))

    # Make sure that ESM data only contains information for a single model.
    assertthat::assert_that(length(unique(esm_data$model)) == 1, msg = 'esm_data can only have input for a single ESM')

    # Make sure that the ESM data only contains values for tas and co2 (other variables may be added in the future).
    assertthat::assert_that(any(unique(esm_data$variable) %in% c('tas', 'co2')), msg = 'esm_data can only contain tas or co2 as variables.')

    # Make sure there is a hector core for each experiment in the ESM data. Then subset the core list so that is only
    # contains the cores there is ESM data for.
    hector_experiments  <- unlist(lapply(hector_cores, function(x){hector::fetchvars(core = x, dates = NA, vars = ECS())[['scenario']]}))
    missing_experiments <- unique(esm_data$experiment)[!unique(esm_data$experiment) %in% hector_experiments]
    message             <- paste0('hector_cores are missing cores for the following esm experiments: ', paste(missing_experiments, collapse = ', '))
    assertthat::assert_that(length(missing_experiments) == 0, msg = message)
    cores_to_use <- hector_cores[which(hector_experiments %in% unique(c(esm_data$experiment)))]

    # Use the order of the experiments in the Hector cores to determine the order of experiments in the esm data.
    experiment_order <- hector_experiments[which(hector_experiments %in% unique(c(esm_data$experiment)))]
    esm_data$experiment <- factor(x = esm_data$experiment, levels = experiment_order, ordered = TRUE)

    # Center and scale the ESM data then split up the normalized esm data.
    esm_data %>%
        dplyr::mutate(id = paste0(experiment, '.', variable, '.', year)) %>%
        dplyr::left_join(tibble::tibble(id = names(normalize$scale), scale = normalize$scale), by = 'id') %>%
        dplyr::left_join(tibble::tibble(id = names(normalize$center), center = normalize$center), by = 'id') %>%
        dplyr::mutate(esm_norm = (value - center)/scale) ->
        normalized_esm

    # Check to make sure there is a center and scale value for each of the entries in the esm data frame.
    normalized_esm %>%
        dplyr::filter_at(.vars = c('center', 'scale'), dplyr::any_vars(is.na(.))) %>%
        dplyr::pull(id) ->
        missing_entries

    # Split up the esm data by experiment.
    normalized_esm %>%
        dplyr::arrange(experiment, variable, year) %>%
        dplyr::select(year, model, variable, experiment, esm_norm, scale, center) %>%
        split(.$experiment) ->
        esm_experiment_list

    # Make a message about the misisng scale and center values.
    message <- paste0('Missing scale and center values for :', paste(missing_entries, collapse = ', '))
    assertthat::assert_that(length(missing_entries) == 0, msg = message)

    # Make the function that will reset the Hector cores with new parameters.
    # This function makes sure that the parameters have been set up in the correct order.
    reset_hector_param_fn <- make_parameterize_cores(param)

    # The function that will be used as fn input to optim.
    function(param){

        # Run the Hector cores and calculate the mean squared error between the Hector output and ESM data.
        out <- foreach::foreach(exp = 1:length(cores_to_use), .combine = 'c') %do% {

            # Pull out the years and variable names from the esm data for the experiment.
            yrs <- unique(esm_experiment_list[[exp]]$year)
            var <- ifelse(unique(esm_experiment_list[[exp]]$variable) == 'co2', hector::ATMOSPHERIC_CO2(), hector::GLOBAL_TEMP())

            # Run the Hector core and calculate the MSE for the Hector and ESM output data.
            MSE <- suppressWarnings(suppressMessages(tryCatch({

                # Reset the core.
                reset_hector_param_fn(core = cores_to_use[[exp]], param = param)

                # Run the Hector core.
                hector::run(core = cores_to_use[[exp]], runtodate = max(yrs))

                # Fetch the output data of intrest for the years to compare with the ESM data and arrange to match order
                # of the information in the ESM data frame.
                hector::fetchvars(core = cores_to_use[[exp]], yrs, var) %>%
                    dplyr::arrange(scenario, variable, year) %>%
                    dplyr::mutate(variable = dplyr::if_else(variable == 'Tgav', 'tas', 'co2')) %>%
                    dplyr::left_join(esm_experiment_list[[exp]],
                                     by = c('scenario' = 'experiment', 'year', 'variable')) %>%
                    na.omit() %>%
                    dplyr::mutate(hector_norm = (value - center) / scale) %>%
                    dplyr::mutate(SE = (hector_norm - esm_norm)^2) %>%
                    dplyr::pull(SE) %>%
                    mean

                }, error = errhandler)))

            if(is.null(MSE)){
                Inf
            } else {
                MSE
            }

        }

        sum(out)
    }
}



