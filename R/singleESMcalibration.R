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


#' Make the function to minimize
#'
#' Create the function that will calculate the mean squared error between Hector and ESM CMIP5 output data. THis
#' function will be minimized using optim.
#'
#' @param hector_cores A list of Hector cores that will be used.
#' @param esm_data A data frame of ESM data for a single model that contains the following columns, year, model, variable, experiment, weight.
#' @param normalize A list of center and the scale values to use to noramlize the Hector and ESM output data.
#' @param param A vector of parameter values.
#' @param core_weights An optional vector of experiment weights to get the weighted sum of experiment MSE, the default it set to NULL the no weights are used.
#' @param n The number of cores to parallelize the Hector runs over, unless sepcified will use the number of cores detected by \code{detectCores}
#' @param showMessages Default set to FALSE, will supress Hector error messages.
#' @return A function that will calculate the mean squared error between a Hector run and ESM data.
#' @importFrom foreach %do% %dopar%
#' @importFrom dplyr %>%
#' @export
make_minimize_function <- function(hector_cores, esm_data, normalize, param, core_weights = NULL, n = NULL, showMessages = FALSE){

    # Check Inputs -----
    # Check inputs for names.
    check_columns(input = normalize, req_cols = c('scale', 'center'))
    check_columns(input = esm_data, req_cols = c('model', 'experiment', 'variable', 'year'))

    # Make sure that the normalize center and scale arguments are all named.
    assertthat::assert_that(all(!is.na(names(normalize$center))), msg = 'normalize center needs names')
    assertthat::assert_that(all(!is.na(names(normalize$scale))), msg = 'normalize scale needs names')

    # Make sure that ESM data only contains information for a single model.
    assertthat::assert_that(length(unique(esm_data$model)) == 1, msg = 'esm_data can only have input for a single ESM')

    # Make sure that the ESM data only contains values for tas and co2 (other variables may be added in the future).
    assertthat::assert_that(any(unique(esm_data$variable) %in% c('tas', 'co2')), msg = 'esm_data can only contain tas or co2 as variables.')

    # Make sure that if weights are povided that there is one for every experiment (hector core). If the
    # weights argument is left to default use 1 as the weight value.
    if(is.null(core_weights)){
        core_weights <- rep(x = 1, length(hector_cores))
    }
    assertthat::assert_that(length(core_weights) == length(hector_cores), msg = 'core_weights and hector_cores must have the same length')


    # Make sure there is a hector core for each experiment in the ESM data. Then subset the core list so that is only
    # contains the cores there is ESM data for.
    hector_experiments <- unlist(lapply(hector_cores, hector::getname))
    missing_experiments <- unique(esm_data$experiment)[!tolower(unique(esm_data$experiment)) %in% tolower(hector_experiments)]
    message             <- paste0('hector_cores are missing cores for the following esm experiments: ', paste(missing_experiments, collapse = ', '))
    assertthat::assert_that(length(missing_experiments) == 0, msg = message)

    # Select the cores and the core weights to use, assume that they are in the same order.
    cores_to_use   <- hector_cores[which(hector_experiments %in% unique(esm_data$experiment))]
    weights_to_use <- core_weights[which(hector_experiments %in% unique(esm_data$experiment))]

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

        ## Organize the runs into batches that will be run in parallel
        nbatch  <- as.integer(floor(length(cores_to_use)/n))
        nextra  <- as.integer(length(cores_to_use)%%n)
        nparams <- rep()

        # Calculate the mean squared error between the Hector output and ESM data for the runs in the
        # batches.
        if(nbatch > 0){

            rslt1 <- foreach::foreach(k = 1:nbatch, .combine = 'c') %do% {

                # Parallel batch
                foreach::foreach(i = 1:n, .combine = 'c') %dopar% {


                    exp <- (k-1) * n + i        # Define the experiment index based  on the setup of the parallel batches.

                    # Pull out the years and variable names from the esm data for the experiment.
                    yrs <- unique(esm_experiment_list[[exp]]$year)
                    var <- ifelse(unique(esm_experiment_list[[exp]]$variable) == 'co2', hector::ATMOSPHERIC_CO2(), hector::GLOBAL_TEMP())

                    # Run the Hector core and calculate the MSE for the Hector and ESM output data.
                    MSE <- tryCatch({
                        # Reset the core with the new parameters.
                        parameterize_core(core = cores_to_use[[exp]], params = param)

                        # Run the Hector core.
                        hector::run(core = cores_to_use[[exp]], runtodate = max(yrs))

                        # Fetch the output data of intrest for the years to compare with the ESM data and arrange to match order
                        # of the information in the ESM data frame.
                        hector::fetchvars(core = cores_to_use[[exp]], yrs, var) %>%
                            dplyr::arrange(scenario, variable, year) %>%
                            dplyr::mutate(variable = dplyr::if_else(variable == 'Tgav', 'tas', 'co2'),
                                          scenario = as.character(scenario)) %>%
                            dplyr::left_join(esm_experiment_list[[exp]] %>%
                                                 dplyr::mutate(experiment = as.character(experiment)),
                                             by = c('scenario' = 'experiment', 'year', 'variable')) %>%
                            na.omit() %>%
                            dplyr::mutate(hector_norm = (value - center) / scale) %>%
                            dplyr::mutate(SE = (hector_norm - esm_norm)^2) %>%
                            dplyr::pull(SE) %>%
                            mean

                    },error=debug_errhandler)

                    if(is.null(MSE)){
                        Inf
                    } else {
                        MSE
                    }

                }

            }

        } else {rslt1 <- NA}


        # Caclculate left over MSE.
        if(nextra > 0) {
            rslt2 <- foreach::foreach(i=1:nextra, .combine='c') %dopar% {

                exp <- nbatch * n + i       # Define the experiment list inex based on batch position.

                # Pull out the years and variable names from the esm data for the experiment.
                yrs <- unique(esm_experiment_list[[exp]]$year)
                var <- ifelse(unique(esm_experiment_list[[exp]]$variable) == 'co2', hector::ATMOSPHERIC_CO2(), hector::GLOBAL_TEMP())

                # Run the Hector core and calculate the MSE for the Hector and ESM output data.
                MSE <- tryCatch({
                    # Reset the core with the new parameters.
                    parameterize_core(core = cores_to_use[[exp]], params = param)

                    # Run the Hector core.
                    hector::run(core = cores_to_use[[exp]], runtodate = max(yrs))

                    # Fetch the output data of intrest for the years to compare with the ESM data and arrange to match order
                    # of the information in the ESM data frame.
                    hector::fetchvars(core = cores_to_use[[exp]], yrs, var) %>%
                        dplyr::arrange(scenario, variable, year) %>%
                        dplyr::mutate(variable = dplyr::if_else(variable == 'Tgav', 'tas', 'co2'),
                                      scenario = as.character(scenario)) %>%
                        dplyr::left_join(esm_experiment_list[[exp]] %>%
                                             dplyr::mutate(experiment = as.character(experiment)),
                                         by = c('scenario' = 'experiment', 'year', 'variable')) %>%
                        na.omit() %>%
                        dplyr::mutate(hector_norm = (value - center) / scale) %>%
                        dplyr::mutate(SE = (hector_norm - esm_norm)^2) %>%
                        dplyr::pull(SE) %>%
                        mean

                },error=debug_errhandler)




                if(is.null(MSE)){
                    Inf
                } else {
                    MSE
                }

            }

        } else {
            rslt2 <- NA
        }

        doParallel::stopImplicitCluster()

        # Calculate the weighted sum for the non NA values.
        final_values <- c(rslt1, rslt2)[!is.na(c(rslt1, rslt2))]
        weighted.mean(x = final_values, w = weights_to_use)

    }
}


#' Emulate a single ESM by calibrating Hector to ESM output data.
#'
#' @param inifiles A vector of the Hector inifile cores to use, there should be a core for each esm CMIP comparison data.
#' @param hector_names A vector of the Hector core names, should reflect experiment names in the esm CMIP comparion data.
#' @param esm_data A data frame of ESM data for a single model that contains the following columns, year, model, variable, experiment.
#' @param normalize A list of center and the scale values to use to noramlize the Hector and ESM output data.
#' @param initial_param A named vector of inital paramters to be optimized over.
#' @param core_weights An optional vector of experiment weights to get the weighted sum of experiment MSE, the default it set to NULL the no weights are used.
#' @param maxit The max number of itterations for optim, default set to 500.
#' @param n_parallel The max number of cores to parallize the runs over,  unless sepcified will use the number of cores detected by \code{detectCores}.
#' @param showMessages Default set to FALSE, will supress Hector error messages.
#' @return An object returned by \code{optim}
#' @export
singleESM_calibration <- function(inifiles, hector_names, esm_data, normalize, initial_param, core_weights = NULL, maxit = 500, n_parallel = NULL, showMessages = FALSE){

    # Set up the Hector cores.
    cores <- setup_hector_cores(inifile = inifiles, name = hector_names)

    # Make the function that will calculate the mean squared error between Hector output and the esm comparison data,
    # this function will be minimized by optim.
    fn <- make_minimize_function(hector_cores = cores, esm_data = esm_data, normalize = normalize, core_weights = core_weights, param = initial_param, n = n_parallel, showMessages = showMessages)

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
#' @param core_weights An optional vector of experiment weights to get the weighted sum of experiment MSE, the default it set to NULL the no weights are used.
#' @param maxit The max number of itterations for optim, default set to 500.
#' @param n_parallel The max number of cores to parallize the runs over,  unless sepcified will use the number of cores detected by \code{detectCores}.
#' @param showMessages Default set to FALSE, will supress Hector error messages.
#' @return A list containing the following elements, copmarison_plot a plot comparing Hector and ESM output data, residual_plot a plot comparing the
#' normalized residuals, MSE a data frame of the mean squared error for each experiment and variable optim minimizes the sum of the MSE values, and
#' optim_rslt is the object returned by \code{optim}.
#' @export
singleESM_calibration_diag <- function(inifiles, hector_names, esm_data, normalize, initial_param, core_weights = NULL, maxit = 500, n_parallel = NULL, showMessages = FALSE){

    # Make an empty list to return the output in.
    output <- list()

    # Parse out information from the ESM comparison data. This will be used to extract data from the Hector cores.
    yrs <- unique(esm_data$year)
    var <- unique(esm_data$variable)
    var <- ifelse(var == 'tas', hector::GLOBAL_TEMP(), hector::ATMOSPHERIC_CO2()) # Rename to Hector strings.
    esm_model_name <- unique(esm_data$model)                                      # Save the esm model name for latter.


    # Get the best parameter fit for Hector and the ESM.
    calibration_rslts <- singleESM_calibration(inifiles = inifiles, hector_names = hector_names, esm_data = esm_data,
                                               normalize = normalize, initial_param = initial_param, core_weights = core_weights,
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
        hector::fetchvars(core = x, dates = yrs, vars = var) %>%
            dplyr::mutate(variable = if_else(variable == hector::GLOBAL_TEMP(), 'tas', 'co2'))

         }) %>%
        dplyr::bind_rows() %>%
        dplyr::rename(experiment = scenario,
                      hector = value) ->
        hector_output

    # Combine the Hector output and esm comparion data into a single data frame.
    esm_data %>%
        dplyr::select(esm = value, experiment, variable, year) %>%
        dplyr::left_join(hector_output, by = c("experiment", "variable", "year")) %>%
        dplyr::mutate(variable = paste0(variable, ' ', units)) %>%
        na.omit() ->
        esm_hector_df

    # Make a caption for the plots out the parameter values.
    param_caption <- paste(paste0(names(calibration_rslts$par), ' = ', signif(calibration_rslts$par, digits = 4)), collapse = ', ')

    # Compare the Hector and ESM output data.
    esm_hector_df %>%
        tidyr::gather(model, value, hector, esm) %>%
        dplyr::mutate(model = dplyr::if_else(model == 'hector', model, esm_model_name)) %>%
        ggplot(aes(year, value, color = model, linetype = experiment)) +
        geom_line() +
        facet_wrap('variable', scales = 'free') +
        labs(title = paste0(esm_model_name, '  vs. Hector Output'),
             caption = param_caption,
             y = NULL,
             x = 'Year') +
        theme_bw() ->
        output[['comparison_plot']]

    # Calculate the residuals
    esm_hector_df %>%
        dplyr::mutate(var1 = substr(variable, 1, 3)) %>%
        dplyr::mutate(index = paste0(experiment, '.', var1, '.', year)) %>%
        dplyr::left_join(tibble::tibble(index = names(normalize$scale),
                                        scale = normalize$scale,
                                        center = normalize$center), by = 'index') %>%
        dplyr::mutate(esm_norm = (esm - center) / scale,
                      hec_norm = (hector - center) / scale,
                      diff = hec_norm - esm_norm) %>%
        na.omit() ->
        esm_hector_df

    # plot the residuals.
    ggplot(data = esm_hector_df) +
        geom_point(aes(year, diff, color = experiment)) +
        facet_wrap('variable', scales = 'free') +
        labs(x = 'Year',
             y = paste0('Difference Between Normalized Hector & ', esm_model_name),
             caption = param_caption,
             title = paste0('Residuals for Hector Calibrated to ', esm_model_name)) +
        theme_bw()->
        output[['residual_plot']]

    # Make a data frame of the MSE for each experiment and variable.
    esm_hector_df %>%
        dplyr::group_by(experiment, variable) %>%
        dplyr::summarise(MSE = mean(diff^2)) %>%
        dplyr::ungroup() ->
        output[['MSE']]

    } else {

        output[['message']] <- 'did not converge'

    }

    # Add the calibration results to the output list.
    output[['optim_rslt']] <- calibration_rslts

    output[['weights']]    <- core_weights


    # Return the output.
    return(output)

}




