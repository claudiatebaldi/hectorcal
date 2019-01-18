
# Required libs
library(hector)
library(mcmc)
library(dplyr)


# Format the params as a matrix that can be used as hector inputs
# args: 
#   params: a vector of the paramters being calibrated by the mcmc
#   carbon_cycle: TRUE or FALSE argument to indicate if the carbon cycle paramters are being calibrated or not
#   calibration: "mean" or "consensus" a string indicating what comparison data to compare hector results with
# returns: a dataframe of the parameters
#
# Internal function
make_params_df <- function(params, carbon_cycle, calibration){
  
  # Check the length of the params vecotr, if it is off then throw and error. 
  expected_length <- c(3, 4, 5, 6, 8) 
  if(!length(params) %in% expected_length){stop('Unexpected length of params')}
  
  # Create a mapping matrix of the params, it is important to note that the order of the paramters matters 
  # here, the params must always start with and include the ECS, AERO_SCALE, and DIFFUSIVITY params
  params_matrix <- matrix(nrow = length(params), ncol = 5, byrow = TRUE)
  params_df <- tibble::as_tibble(params_matrix)
  params_df[1, ] <- c(NA, ECS(), params[1], 'degC', 'hector')          # ECS
  params_df[2, ] <- c(NA, AERO_SCALE(), params[2], '(unitless)','hector')   # AERO_SCALE 
  params_df[3, ] <- c(NA, DIFFUSIVITY(), params[3], 'cm2/s', 'hector') # DIFFUSIVITY
  
  if(carbon_cycle){
    # If the carbon_cycle argument is true then the next paramters in the params_matrix will be 
    # the carbon cycle parameters, PREINDUSTRIAL_CO2(), Q10_RH(), and BETA(). 
    params_df[4, ] <- c(NA, PREINDUSTRIAL_CO2(), params[4], 'ppmv CO2', 'hector') # PREINDUSTRIAL_CO2
    params_df[5, ] <- c(NA, Q10_RH(), params[5], '(unitless)', 'hector')          # Q10_RH
    params_df[6, ] <- c(NA, BETA(), params[6], '(unitless)', 'hector')            # BETA
    
  } 
  
  if (carbon_cycle == FALSE & calibration == 'mean'){
    # If the carbon_cycle paramters are not being used and we are calibrating to the comparison 
    # data mean value then params 4 is the sigma_tgav, NOTE this is not a hector paramter and there is 
    # no sigma_co2. 
    params_df[4, ] <- c(NA, 'sigma_tgav', params[4], NA, NA)  # variability of tgav
    
  } 
  
  if (carbon_cycle == TRUE & calibration == 'mean'){
    # If the carbonc_cycle parameters are being used and we are calibrating to the comparison 
    # data mean value then params 6 and 7 are simga_c02 and sigma_tgav, these are NOT Hector 
    # paramters. 
    params_df[7, ] <- c(NA, 'sigma_tgav', params[7], NA, NA)  # variability of tgav
    params_df[8, ] <- c(NA, 'sigma_co2', params[8], NA, NA)   # varibaility of co2
    
  }
  
  names(params_df) <- c('date', 'param', 'value', 'units', 'source')
  
  # If we use a better method to make this then we could drop this step
  mutate(params_df, value = as.numeric(value)) %>%  
    filter(!is.na(value))
  
}


# Calculate the log prior probability density 
#
# Args: 
#   params_df: a data frame of the mcmc parameters created by make_params_df
# Returns: a vector of the log prior probabilty density for each parameter 
#
# ENHANCEMENT IDEA -- I would like this to be something that users can decide what a param probabilty looks like
# also would like to use a different sources for priors
#
# Internal function
get_log_priors <- function(params_df){
  
  # All priors are coming from https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1600-0870.2010.00471.x
  # will want to revist how we define them post PI meeting 
  
  # Make an empty vector that will be returned as output
  num_param <- nrow(params_df)
  out <- rep(NA, num_param)
  
  # The first three parameters from the params martix should be 
  # ECS, AERO_SCALE, and DIFFUSIVITY 
  sampled_ecs <- params_df$value[1]
  out[1]      <- dnorm(sampled_ecs, mean = 5.05, sd = 6, log = TRUE) # range of 0.1 to 10 reported
  
  sampled_aero <- as.numeric(params_df[2, 3])
  out[2]       <- dnorm(sampled_aero, mean = 1, sd = 1.4, log = TRUE) # range of 0 to 2 reported
  
  sampled_diff <- as.numeric(params_df[3, 3])
  out[3]       <- dnorm(sampled_diff, mean = 2.05, sd = 2.7, log = TRUE ) # range of .1 to 4
  
  # The remaining order of the params the in matrix depends on the carbon cycle being turned on or off
  # and the calibration method being used. Now determine the prior based on the parameter's array index. 
  if(PREINDUSTRIAL_CO2() %in% params_df$param){
    
    pCO2_ind <- which(params_df$param == PREINDUSTRIAL_CO2())
    out[pCO2_ind] <- dnorm(params_df$value[pCO2_ind], mean = 285, sd = 14.1, log = TRUE) # range 275 to 295
    
  }
  
  if(Q10_RH() %in% params_df$param){
    
    q10_ind <- which(params_df$param ==  Q10_RH())
    out[q10_ind] <- dnorm(params_df$value[q10_ind], mean = 2.6, sd = 3.4, log = TRUE) # range 0.2 to 5 
    
  }
  
  if(BETA() %in% params_df$param){
    
    beta_ind <- which(params_df$param ==  BETA())
    out[beta_ind] <- log(truncnorm::dtruncnorm(params_df$value[beta_ind],  a=0, b=Inf, mean = 0.5, sd = 0.7))
    
  }
  
  if('sigma_tgav' %in% params_df$param){
    
    sigmaTgav_ind <- which(params_df$param == 'sigma_tgav')
    prior_sigma <- dcauchy(params_df$value[sigmaTgav_ind], location = 0, scale = 1)
    # Since we are instrested in the half cauchy distribtuion reset anything 
    # less than 0 to 0 before determining the log prior
    prior_sigma <- ifelse(prior_sigma <= 0, 0, prior_sigma)
    
    out[sigmaTgav_ind] <- log(prior_sigma)
    
  }
  
  if('sigma_co2' %in% params_df$param){
    if(params_df$value[sigmaCo2_ind]) {
        out[sigmaCo2_ind] <- -Inf
    }
    else {
      
        sigmaCo2_ind <- which(params_df$param == 'sigma_co2')
        prior_sigma <- dcauchy(params_df$value[sigmaCo2_ind], location = 0, scale = 10)
        # Since we are instrested in the half cauchy distribtuion reset anything 
        # less than 0 to 0 before determining the log prior
        prior_sigma <- ifelse(prior_sigma <= 0, 0, prior_sigma)
    
        out[sigmaCo2_ind] <- log(prior_sigma)
    
    }
  }
  
  out 
  
}


# Thes
## This thing we're calling "erf" is missing the sqrt(2) factor inside the pnorm
## call. Our calls *should* be erf((b-a)/(sqrt(2)*sig)), so dropping the sqrt(2)
## here allows us to drop it to from the args too, saving us a lot of sqrt(2)
## factors in our code that will only end up canceling out anyway.  It does,
## however, mean that this "erf" is not technically the real erf.
erf <- function(x) {2*pnorm(x)-1}
consensus_mesa <- function(x, a, b, sig) {(erf((b-x)/sig) - erf((a-x)/sig))/(2*(b-a))}

# Generate the log posteior function 
#
# Args:
#   ini_file: a string dir of the ini file to use in the hector core 
#   cdata: a data frame of data that will be compared with Hector output. The contents of this data frame depends on identify of 
#     calibration. If cdata does not meet calibration requirements then an error will be thrown. 
#   cabron_cycle: T/F to indicate if calibrating to co2, if set to F then carbon parameters should not be sampled by metrop
#   calibration: a string indicating what calibration method to used, mean or consensus (IDEA: change to some sort of function)
#   step: the number of years to skip in when comparing data to decrease/prevent correlation between point
#   showMessages: if set to F then using the function returned by make_log_post_function will supress messages
# returns: the log posterior function, this function requires a vector of parameters. 
make_log_post_function <- function(ini_file, cdata, carbon_cycle, calibration, step, showMessages = TRUE){
  
  # Make a hector core
  core <- hector::newcore(ini_file, suppresslogging = TRUE)
  
  # Check the calibration methods 
  # TODO this will probably change when the calibration method argument switches from a string to a function
  if(!calibration %in% c('mean', 'consensus')){ stop('calibration method selected is not supported')}
  
  
  if(calibration == 'mean'){
    
    req_cols <- c('cmean', 'year', 'experiment')
    missing  <- !req_cols %in% names(cdata)
    if(any(missing)){
      stop('cdata is missing ', paste(req_cols[missing], collapse = ', ', ' which are required for the mean calibration method'))
    }
    
  }
  
  
  if(calibration == 'consensus'){
    
    req_cols <- c('mina', 'maxb', 'experiment')
    missing  <- !req_cols %in% names(cdata)
    if(any(missing)){
      stop('cdata is missing ', paste(req_cols[missing], collapse = ', ', ' which are required for the consensus calibration method'))
    }
    
  }
  
  
  # Check the variables of the cdata 
  if(!'tas' %in% cdata$variable){stop('cdata is missing variable tas')}

  if(carbon_cycle){
    
    if(!'co2' %in% cdata$variable){stop('cdata is missing variable co2 which is required when carbon_cycle = TRUE')}
    
  }
  
  
  
  # Make sure that the comparison data has results from a single experiment 
  # TODO this requiremetn might change 
  if(length(unique(cdata$experiment)) > 1){ stop('cdata can only contain results from a single experiment')}
  
  
  
  function(params){
    
    if(!any(calibration %in% c('mean', 'consensus'))) {
      
      stop('only mean and consensus calibration are currently supported')
      
    }
    
    # First calibrate make the params matrix
    params_df <- make_params_df(params, carbon_cycle, calibration)
    
    if(showMessages){
      
     # message(params_df, appendLF = TRUE)
      #cat(params_df)
      print(params_df)
      
    }
    
    
    # Then claculate the log prior for the sampled parameters, we are going 
    # to multiply the log liklihood by these values.
    lprior <- get_log_priors(params_df)
    

    # If the sum pof the lprior is -Inf the bail and continue to the next itteration
    if(is.infinite(sum(lprior))){

      if(showMessages) message('lprior contains -Inf\n', paste(lprior, collapse = ', '))
      return(-Inf)

      }

    
    # Determine which entires in the params_matrix are hector (not sigma parameters) and 
    # use those entries to calibrate the hector core.
    hector_params <- which(params_df$source == 'hector')
    
    # Set the hector parameters
    for(i in hector_params){
      
      # Reset the hector core, date, paramater, value, and units for each hector parameter in the params matrix.
      setvar(core = core, dates = NA, var = params_df[['param']][i], values = params_df[['value']][i], unit = params_df[['units']][i])
      
    }
    
    # Reset the core so that spin up uses the new param values
    safe_reset <- purrr::safely(reset)
    rslt       <- safe_reset(core)
    if(is.null(rslt$result)){ 
      
      if(showMessages) message('hector core failed to reset\n')
      if(showMessages) message(rslt$error)
      if(showMessages) message(core)
      
      return(-Inf) 
      
      }
    
    # Safely run hecotr, if the parameters cause Hector to fail return -Inf for the posterior
    safe_hector <- purrr::safely(hector::run)
    rslt        <- safe_hector(core)
    if(is.null(rslt$result)){ 
      
      if(showMessages) message('hector core failed to solve')
      
      return(-Inf) 
      }
    
    # Determine the years to compare hector vs output data with 
    start_yr <- min(cdata$year)
    end_yr   <- max(cdata$year)
    cyears   <- seq(from = start_yr, to = end_yr, by = step)
    
    if(showMessages) message(paste(cyears, collapse = ', '))
    
    # Extract hector results 
    hector::fetchvars(core, dates = cyears, vars = c(GLOBAL_TEMP(), ATMOSPHERIC_CO2())) %>% 
      dplyr::arrange(year) -> 
      hector_all
    #hector::shutdown(core)
    
    hector_tgav <- hector_all[hector_all$variable == GLOBAL_TEMP(), ]
    hector_co2  <- hector_all[hector_all$variable ==  ATMOSPHERIC_CO2(), ]
    
    # Format the esm comparison data that will be used to calculate the log likelihood
    cdata <- dplyr::arrange(cdata, year)
    esm_tgav <- cdata[cdata$year %in% cyears & cdata$variable == 'tas', ]
    if(carbon_cycle){ 
      esm_co2 <- cdata[cdata$year %in% cyears & cdata$variable == 'co2', ] 
    }
    
    # Calculate the log likelihood from the hector and esm results base.d on the 
    # calibration method.
    if(calibration == 'mean'){
      
      sigma_tgav <- params_df$value[which(params_df$param == 'sigma_tgav')]
      diff       <- hector_tgav$value - esm_tgav$cmean; diff
      ll_tgav    <- dnorm(x = diff , mean = 0, sd = sigma_tgav, log = TRUE); ll_tgav
      ll_co2     <- 0
      
      if(showMessages) message('tgav diff ', paste(diff, collapse = ', '))
      
      if(carbon_cycle){
        
        simga_co2 <- params_df$value[which(params_df$param == 'sigma_co2')]
        diff      <- hector_co2$value - esm_co2$cmean
        ll_co2    <- dnorm(x = diff, mean = 0, sd = simga_co2, log = TRUE)
        
        if(showMessages) message('co2 diff ', paste(diff, collapse = ', '))
        
      } 
      
      
    }
    
    if(calibration == 'consensus'){
      
      mesa_out <- consensus_mesa(hector_tgav$value, esm_tgav$mina, esm_tgav$maxb, .4)
      ll_tgav <- sum(log(mesa_out)) # WANT TO CHECK OUT THE SIGMA 
      ll_co2 <- 0
      if(showMessages) message('mesa out tgav ', mesa_out)
      
      if(carbon_cycle){
        
        mesa_out <- consensus_mesa(hector_co2$value, esm_co2$mina, esm_co2$maxb, 15)
        ll_co2 <- sum(log(mesa_out)) # WANT TO CHECK OUT THE SIGMA 
        if(showMessages) {  
          message('mesa out co2 ', mesa_out)
          message('\n\n')}
        
      }
      
    }
    
    # Return the sum of log likelihood * log prior
    if(any(is.infinite(c(lprior, ll_tgav, ll_co2)))){
      
      if(showMessages){
        
        message('lprior ', paste(lprior, collapse = ', '))
        message('ll_tgav ', paste(ll_tgav, collapse = ', '))
        message('ll_co2 ', paste(ll_co2, collapse = ', '))
        message('\n\n')
        
        
      }
      
      return(-Inf)
      
    }
    
    if(showMessages){
      
      message(paste('log like tgav', ll_tgav))
      message(paste('log like co2', ll_co2))
      message('-------')
    }

    
    # log prior * log liklihood 
    sum(lprior, ll_tgav, ll_co2)
    
  }
  
  
}





