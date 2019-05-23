#### Functions used in calibrating Hector to emulate an individual ESM.

#' Parameterize a Hector core.
#'
#' @param params A named vector of Hector parameters.
#' @param core A Hector core.
#' @return A Hector set up to run with a new set of paramters.
parameterize_core <- function(params, core) {

    pnames <- names(params)
    assertthat::assert_that(!is.null(pnames), msg = 'params must be named.')

    punits <- hector::getunits(pnames)
    assert_that(!any(is.na(punits)), msg='Bogus parameter names')

    mapply(function(x, y, z){hector::setvar(core = core, dates = NA, var = x, values = y, unit = z)},
           x = pnames, y = params, z = punits)

    hector::reset(core = core)

}


#' Switch from Hector variable names to ESM variable names.
#'
#' Hector and ESM use different nomenclature for the variables, this function changes the Hector
#' variable name to the esm variable name for (tas, co2, and heatflux but other variables may be added).
#'
#' @param hector_output A dataframe of Hector output data.
#' @importFrom dplyr %>%
#' @return A dataframe of Hector output data but with the ESM variable names.
translate_variable_name <- function(input){

    # Make the mapping file for the variable names.
    variable_mapping_tib <- tibble::tibble(variable = c(hector::GLOBAL_TEMP(), hector::ATMOSPHERIC_CO2(), hector::HEAT_FLUX()),
                                           esm_variable = c('tas', 'co2', 'heatflux'))

    # Make sure that the Hector output data only contains data for the variables in the mapping tibble.
    if(all(input$variable %in% variable_mapping_tib$variable)){

        # Replace the Hector variable names with the ESM variable names.
        input %>%
            dplyr::left_join(variable_mapping_tib, by = 'variable') %>%
            dplyr::select(-variable) %>%
            dplyr::rename(variable = esm_variable)

    } else if (all(input$variable %in% variable_mapping_tib$esm_variable)) {

        # Replace the Hector variable names with the ESM variable names.
        input %>%
            dplyr::rename(esm_variable = variable) %>%
            dplyr::left_join(variable_mapping_tib, by = 'esm_variable') %>%
            dplyr::select(-esm_variable)

    } else {

        stop('There is some problem with how the variable names are being mapped between Hector and ESM variable names')

    }


}



#' Make the function to minimize
#'
#' Create the function that will calculate the mean squared error between Hector and ESM CMIP5 output data. THis
#' function will be minimized using optim.
#'
#' @param hector_cores A list of Hector cores that will be used.
#' @param esm_data A data frame of ESM data for a single model that contains the following columns, year, model, variable, experiment, weight.
#' @param normalize A list of center and the scale values to use to noramlize the Hector and ESM output data.
#' @param param A vector of parameter values.
#' @param cmip_range Default set to NULL and the make_minim_function will only look at the error between esm_data and Hector data.
#' But if cmip_range is set to a dataframe containing columns (variable, year, sig, lower, and upper) then the minimize function will also
#' minimize the value returned by the -log of the mesa function.
#' @param n The number of cores to parallelize the Hector runs over, unless sepcified will use the number of cores detected by \code{detectCores}
#' @param showMessages Default set to FALSE, will supress Hector error messages.
#' @return A function that will calculate the mean squared error between a Hector run and ESM data.
#' @importFrom foreach %do% %dopar%
#' @importFrom dplyr %>%
#' @export
make_minimize_function <- function(hector_cores, esm_data, normalize, param, cmip_range = NULL, n = NULL, showMessages = FALSE){

    # Check data used to normalize the Hector and ESM output data.
    check_columns(input = normalize, req_cols = c('scale', 'center'))
    check_columns(input = esm_data, req_cols = c('model', 'experiment', 'variable', 'year'))
    assertthat::assert_that(all(!is.na(names(normalize$center))), msg = 'normalize center needs names')
    assertthat::assert_that(all(!is.na(names(normalize$scale))), msg = 'normalize scale needs names')

    # Check the ESM comparison data.
    check_columns(input = comp_data, req_cols = c('year', 'model', 'ensemble', 'variable', 'experiment', 'value', 'cmip_experiment'))
    assertthat::assert_that(length(unique(esm_data$model)) == 1, msg = 'esm_data can only have input for a single ESM')
    assertthat::assert_that(all(esm_data$variable %in% c('tas', 'co2', 'heatflux')), msg = 'esm_data can only contain tas, heatflux, or co2 as variables.')

    # Store a list of the experiments that are going to be used as the hector core names
    comparison_experiments <-


    # If calibrating to the cmip range check that the cmip range df contains the required columns and that the
    # the esm_data and the cmip_range objects do not contains values for the same variables.
    # TODO or should this check only prevent the use of esm and the cmip range for the years / experiment or variables?
    if(!is.null(cmip_range)){

        check_columns(input = cmip_range, req_cols = c('year', 'variable', 'cmip_experiment', 'lower', 'upper', 'sig'))
        assertthat::assert_that(!any(esm_data$variable %in% cmip_range$variable), msg = 'cmip_range cannot contain data for a variable that is also in esm_data')
        assertthat::assert_that(any(unique(cmip_range$variable) %in% c('tas', 'co2', 'heatflux')), msg = 'cmip_range can only contain tas, heatflux, or co2 as variables.')

    }



    # Make sure there is a hector core for each experiment in the ESM data. Then subset the core list so that is only
    # contains the cores there is ESM data for.
    hector_experiments <- unlist(lapply(hector_cores, hector::getname))
    missing_experiments <- unique(esm_data$experiment)[!tolower(unique(esm_data$experiment)) %in% tolower(hector_experiments)]
    message             <- paste0('hector_cores are missing cores for the following esm experiments: ', paste(missing_experiments, collapse = ', '))
    assertthat::assert_that(length(missing_experiments) == 0, msg = message)

    # Select the cores and the core weights to use, assume that they are in the same order.
    cores_to_use   <- hector_cores[which(hector_experiments %in% unique(esm_data$experiment))]

    # Prep Data -----
    # Use the order of the experiments in the Hector cores to determine the order of experiments in the esm data.
    experiment_order    <- hector_experiments[which(hector_experiments %in% unique(esm_data$experiment))]
    esm_data$experiment <- factor(x = esm_data$experiment, levels = experiment_order, ordered = TRUE)

    # Center and scale the ESM data then split up the normalized esm data.
    esm_data %>%
        dplyr::mutate(id = paste0(experiment, '.', variable, '.', year)) %>%
        dplyr::left_join(tibble::tibble(id = names(normalize[['scale']]), scale = normalize[['scale']]), by = 'id') %>%
        dplyr::left_join(tibble::tibble(id = names(normalize[['center']]), center = normalize[['center']]), by = 'id') %>%
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
        dplyr::mutate(cmip_experiment = base::gsub(x = experiment, pattern = '_[a-zA-Z0-9-]{6}', replacement = '')) %>%
        split(.$experiment) ->
        esm_experiment_list

    # Make a message about the misisng scale and center values.
    message <- paste0('Missing scale and center values for :', paste(missing_entries, collapse = ', '))
    assertthat::assert_that(length(missing_entries) == 0, msg = message)


    # Define Internal Functions -----
    # The debug_errhandler is designed to return messages from inside the tryCatch statements if
    # the show messages statement is true.
    #
    # Args:
    #    e - an error message thrown by code within a tryCatch statement.
    # Returns: NULL
    debug_errhandler <- function(e){

        if(showMessages){

            message(conditionMessage(e))

        }
        NULL
    }

    # Set Up Parallel Cores ------
    # Parllel runs set up.
    if(is.null(n)){
        n <- parallel::detectCores()}
    if(n > length(cores_to_use)){
        n <- length(cores_to_use)
    }

    doParallel::registerDoParallel(cores=n)


   # The function that will be used as fn input to optim.
    function(param){

        if(showMessages){message(paste0(signif(param, 3), collapse = ', '))}

        # Run the Hector in parallel for all of the scnearios and calculate the
        # mean squared error for each variable / experiment / ensemble.
        rslt <- foreach::foreach(i = 1:length(cores_to_use), .combine = 'rbind') %dopar% {

            # Pull out the years and variable names from the esm data for the experiment.
            yrs <- unique(esm_experiment_list[[i]]$year)
            var <- unique(dplyr::pull(translate_variable_name(esm_experiment_list[[i]]), variable))

            # Run the Hector core and calculate the MSE for the Hector and ESM output data.
            MSE <- tryCatch({
                # Reset the core with the new parameters.
                parameterize_core(core = cores_to_use[[i]], params = param)

                # Run the Hector core.
                hector::run(core = cores_to_use[[i]], runtodate = max(yrs))

                # Fetch the output data of intrest for the years to compare with the ESM data and arrange to match order
                # of the information in the ESM data frame.
                hector::fetchvars(core = cores_to_use[[i]], yrs, var) %>%
                    dplyr::arrange(scenario, variable, year) %>%
                    dplyr::mutate(scenario = as.character(scenario)) %>%
                    translate_variable_name %>%
                    dplyr::left_join(esm_experiment_list[[i]] %>%
                                         dplyr::mutate(experiment = as.character(experiment)),
                                     by = c('scenario' = 'experiment', 'year', 'variable')) %>%
                    stats::na.omit() %>%
                    dplyr::mutate(hector_norm = (value - center) / scale) %>%
                    dplyr::mutate(SE = (hector_norm - esm_norm)^2) %>%
                    dplyr::group_by(scenario, variable) %>%
                    dplyr::summarise(value = mean(SE)) %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(cmip_experiment = base::gsub(x = scenario, pattern = '_[a-zA-Z0-9-]{6}', replacement = '')) ->
                    rslt_esm_comparison

                core_cmip_exp <- unique(rslt_esm_comparison$cmip_experiment)


                # If also calibrating to the cmip rnage then calculate the -log(mesa) for the cmip range.
                rslt_cmip_range_comparison <- NULL
                if(!is.null(cmip_range) & any(core_cmip_exp %in% cmip_range$cmip_experiment)){

                    range_yrs        <- cmip_range$year
                    range_vars       <- unique(pull(translate_variable_name(cmip_range), variable))

                    tibble::as_tibble(hector::fetchvars(core = cores_to_use[[i]], range_yrs, range_vars)) %>%
                        dplyr::mutate(cmip_experiment = base::gsub(x = scenario, pattern = '_[a-zA-Z0-9-]{6}', replacement = '')) %>%
                        dplyr::left_join(cmip_range,  by = c("cmip_experiment", "year", "variable")) %>%
                        na.omit ->
                        hector_cmip_range

                    hector_cmip_list <- split(hector_cmip_range, list(hector_cmip_range$scenario, hector_cmip_range$variable), drop = TRUE)

                    lapply(hector_cmip_list, function(input){

                        mesa_value <- -log(mesa(x = input$value, a = input$lower, b = input$upper, sig = input$sig))

                        input %>%
                            dplyr::select(scenario, year, variable) %>%
                            dplyr::bind_cols(value = mesa_value) %>%
                            dplyr::group_by(scenario, variable) %>%
                            dplyr::summarise(value = mean(value)) %>%
                            dplyr::ungroup() %>%
                            dplyr::mutate(cmip_experiment = unique(core_cmip_exp))

                    }) ->
                        rslt_cmip_range_comparison
                }

                dplyr::bind_rows(rslt_esm_comparison, rslt_cmip_range_comparison)



            },error=debug_errhandler)

            if(is.null(MSE)){
                tibble::tibble(value = Inf,
                               scenario = NA_character_,
                               variable = NA_character_)
            } else {
                MSE
            }

        }


        # Clean up partallel clusters.
        doParallel::stopImplicitCluster()

        # Calculate the value to minimize for each experiment / variable.
        # Take the average for each experiment / variable where variables and ensemble memebers have equal weight to the
        # value to minimize
        dplyr::bind_rows(rslt) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(experiment = gsub(pattern = '_[a-zA-Z0-9-]{6}', replacement = '', x = scenario)) %>%
            dplyr::group_by(experiment, variable) %>%
            dplyr::summarise(value = mean(value)) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(experiment) %>%
            dplyr::summarise(value = sum(value)) %>%
            dplyr::ungroup() %>%
            dplyr::pull(value) %>%
            mean()
    }
}


#' Emulate a single ESM by calibrating Hector to ESM output data.
#'
#' @param inifiles A vector of the Hector inifile cores to use, there should be a core for each esm CMIP comparison data.
#' @param hector_names A vector of the Hector core names, should reflect experiment names in the esm CMIP comparion data.
#' @param esm_data A data frame of ESM data for a single model that contains the following columns, year, model, variable, experiment.
#' @param normalize A list of center and the scale values to use to noramlize the Hector and ESM output data.
#' @param initial_param A named vector of inital paramters to be optimized over.
#' @param cmip_range Default set to NULL and the make_minim_function will only look at the error between esm_data and Hector data.
#' But if cmip_range is set to a dataframe containing columns (variable, year, sig, lower, and upper) then the minimize function will also
#' minimize the value returned by the -log of the mesa function.
#' @param maxit The max number of itterations for optim, default set to 500.
#' @param n_parallel The max number of cores to parallize the runs over,  unless sepcified will use the number of cores detected by \code{detectCores}.
#' @param showMessages Default set to FALSE, will supress Hector error messages.
#' @return An object returned by \code{optim}
#' @export
singleESM_calibration <- function(inifiles, hector_names, esm_data, normalize, initial_param, cmip_range = NULL, maxit = 500, n_parallel = NULL, showMessages = FALSE){

    # Set up the Hector cores.
    cores <- setup_hector_cores(inifile = inifiles, name = hector_names)

    # Make the function that will calculate the mean squared error between Hector output and the esm comparison data,
    # this function will be minimized by optim.
    fn <- make_minimize_function(hector_cores = cores,
                                 esm_data = esm_data,
                                 normalize = normalize,
                                 param = initial_param,
                                 cmip_range = cmip_range,
                                 n = n_parallel,
                                 showMessages = showMessages)

    # Use optim to minimize the MSE between Hector and ESM output data
    stats::optim(par = initial_param, fn = fn, control = list('maxit' = maxit))

}


#' Emulate a single ESM by calibrating Hector to ESM output data and return diagnostic materials.
#'
#' Calibrate Hector so that is emulates a single ESM. Calcaulte the mean squared error between the calibrated Hector and the
#' ESM comparison data, plot the Hector output vs the ESM comparison data, and a plot of the residuals.
#'
#' @param inifiles A vector of the Hector inifile cores to use, there should be a core for each esm CMIP comparison data.
#' @param hector_names A vector of the Hector core names, should reflect experiment names in the esm CMIP comparion data.
#' @param esm_data A data frame of ESM data for a single model that contains the following columns, year, model, variable, experiment.
#' @param normalize A list of center and the scale values to use to noramlize the Hector and ESM output data.
#' @param initial_param A named vector of inital paramters to be optimized over.
#' @param cmip_range Default set to NULL and the make_minim_function will only look at the error between esm_data and Hector data.
#' But if cmip_range is set to a dataframe containing columns (variable, year, sig, lower, and upper) then the minimize function will also
#' minimize the value returned by the -log of the mesa function.
#' @param maxit The max number of itterations for optim, default set to 500.
#' @param n_parallel The max number of cores to parallize the runs over,  unless sepcified will use the number of cores detected by \code{detectCores}.
#' @param showMessages Default set to FALSE, will supress Hector error messages.
#' @return A list containing the following elements, copmarison_plot a plot comparing Hector and ESM output data, residual_plot a plot comparing the
#' normalized residuals, MSE a data frame of the mean squared error for each experiment and variable optim minimizes the sum of the MSE values, and
#' optim_rslt is the object returned by \code{optim}.
#' @export
singleESM_calibration_diag <- function(inifiles, hector_names, esm_data, normalize, initial_param, cmip_range = NULL, maxit = 500, n_parallel = NULL, showMessages = FALSE){

    # Make an empty list to return the output in.
    output <- list()

    # Parse out information from the ESM comparison data. This will be used to extract data from the Hector cores.
    yrs <- unique(esm_data$year)
    var <- unique(dplyr::pull(translate_variable_name(esm_data), variable))
    esm_model_name <- unique(esm_data$model)                                      # Save the esm model name for latter.


    # Get the best parameter fit for Hector and the ESM.
    calibration_rslts <- singleESM_calibration(inifiles = inifiles, hector_names = hector_names, esm_data = esm_data,
                                               normalize = normalize, initial_param = initial_param,
                                               cmip_range = cmip_range,
                                               maxit = maxit, n = n_parallel, showMessages = showMessages)

    if(calibration_rslts$convergence == 0){

        # Make a new set of Hector cores that will be parameterized with the results returned from singleESM_calibration.
        cores <- setup_hector_cores(inifile = inifiles, name = hector_names)

        # Run all of the Hector cores with the new parameters and extract results.
        lapply(cores, function(x){

            # Reset the Hector cores.
            parameterize_core(core = x, params = calibration_rslts$par)

            # Run Hector.
            hector::run(core = x, runtodate = max(yrs))

            # Extract the results.
            translate_variable_name(hector::fetchvars(core = x, dates = yrs, vars = var))

        }) %>%
            dplyr::bind_rows() %>%
            dplyr::rename(experiment = scenario,
                          hector = value) ->
            hector_output

        # Combine the Hector output and esm comparion data into a single data frame.
        esm_data %>%
            dplyr::select(esm = value, experiment, variable, year) %>%
            dplyr::ungroup() %>%
            dplyr::left_join(hector_output, by = c("experiment", "variable", "year")) %>%
            dplyr::mutate(copy_var = variable,
                          variable = paste0(variable, ' ', units)) %>%
            na.omit() ->
            esm_hector_df

        # Make a caption for the plots out the parameter values.
        param_caption <- paste(paste0(names(calibration_rslts$par), ' = ', signif(calibration_rslts$par, digits = 4)), collapse = ', ')

        # Compare the Hector and ESM output data.
        esm_hector_df %>%
            tidyr::gather(model, value, hector, esm) %>%
            dplyr::mutate(model = dplyr::if_else(model == 'hector', model, esm_model_name)) %>%
            dplyr::mutate(cmip_experiment = gsub(pattern = '_[a-zA-Z0-9-]{6}', replacement = '', x = experiment))->
            to_plot

        # Define the colors and assign factor levels to use in plotting.
        model_names   <- c('hector', unique(dplyr::pull(dplyr::filter(to_plot, model != 'hector'), model)))
        colors        <- c('blue', 'grey')
        names(colors) <- model_names

        to_plot$model <- factor(to_plot$model, levels = rev(model_names), ordered = TRUE)

        ggplot(data = to_plot, aes(year, value, color = model, linetype = cmip_experiment, group = interaction(experiment, model))) +
            geom_line() +
            facet_wrap('variable', scales = 'free_y') +
            labs(title = paste0(esm_model_name, '  vs. Hector Output'),
                 caption = param_caption,
                 y = NULL,
                 x = 'Year') +
            theme_bw() +
            scale_color_manual(values = colors) ->
            output[['comparison_plot']]

        # Calculate the residuals
        esm_hector_df %>%
            dplyr::ungroup() %>%
            dplyr::mutate(index = paste0(experiment, '.', copy_var, '.', year)) %>%
            dplyr::left_join(tibble::tibble(index = names(normalize$scale),
                                            scale = normalize$scale,
                                            center = normalize$center), by = 'index') %>%
            dplyr::mutate(esm_norm = (esm - center) / scale,
                          hec_norm = (hector - center) / scale,
                          diff = hec_norm - esm_norm) %>%
            na.omit() ->
            esm_hector_df

        # Plot the residuals.
        ggplot(data = esm_hector_df) +
            geom_point(aes(year, diff, color = experiment)) +
            facet_wrap('variable', scales = 'free') +
            labs(x = 'Year',
                 y = paste0('Difference Between Normalized Hector & ', esm_model_name),
                 caption = param_caption,
                 title = paste0('Residuals for Hector Calibrated to ', esm_model_name)) +
            theme_bw()->
            output[['residual_plot']]

        if(!is.null(cmip_range)){
            output[['extra_data']] <- cmip_range
            }


        # Make a data frame of the MSE for each experiment and variable.
        esm_hector_df %>%
            dplyr::ungroup() %>%
            dplyr::group_by(experiment, variable) %>%
            dplyr::summarise(value = mean(diff^2)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(cmip_experiment = base::gsub(experiment, pattern = '_[a-zA-Z0-9-]{6}', replacement = '')) %>%
            dplyr::group_by(experiment, variable) %>%
            dplyr::summarise(value = mean(value)) %>%
            dplyr::ungroup() ->
            output[['MSE']]

    } else {

        output[['message']] <- 'did not converge'

    }

    # Add the calibration results to the output list.
    output[['optim_rslt']] <- calibration_rslts


    # Return the output.
    return(output)

}




