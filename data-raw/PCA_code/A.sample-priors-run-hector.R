# This script creates an ensemble of Hector runs that will be for the PCA.
# Right now this is set up to run on pic.

# 0. Set Up --------------------------------------------------------------------------------
Sys.time()

pic_dir <- '.'

# Load the required packages
library(hector)
library(hectorcal)
library(doParallel)
library(dplyr)


# Define the prior hyper-parameters, these should be the same as the hyper-parameters defined in
# the likelihood function. These will be used as default inputs into functions defined
# in section 1 to create the prior distirbtuion.
hyper_params <- list(ecsmu = 3,.0,
                     ecssig = 3.0,
                     aeromu = 1.0,
                     aerosig = 1.4,
                     kappamu = 2.3,
                     kappasig = 2.0,
                     betamu = 0.3,
                     betasig = 0.7,
                     q10mu = 2.0,
                     q10sig = 2.0,
                     c0mu = 285,
                     c0sig = 14.0)


# Set the seed to ensure reproducible results
set.seed(867-5309)

# This should be equal to the number of nodes selected in the
nodes <- 6

# The number of Hector runs per rcp scenario
n_runs <- 2000

# Years to extract the Hector output for
years  <- 1850:2100


# 1. Define functions ---------------------------------------------------------------------
# Randomly generate carbon values based on the prior distribution. These are the
# parameters used in the emission driven runs. Will be used in the hector_sample function.
#
# Args
#   ecsmu: the climate sensitivity mean, default is set to the value in the hyper parameter defined in section 0
#   ecssig: the climate sensitivity  standard deviation
#   aeromu: the aerosol scalar mean
#   aerosig: the aerosol standard deviation
# Returns: a list of functions that will sample the climate sensitivity, aero scaler, and ocean diffusivity
# parameter spaces. Each function in the list will sample the prior distribution n times.
make_concen_dist <- function(ecsmu = hyper_params$ecsmu, ecssig = hyper_params$ecssig,
                             aeromu = hyper_params$aeromu, aerosig = hyper_params$ecssig,
                             kappamu = hyper_params$kappamu, kappasig = hyper_params$kappasig){

    list <- list(function(n) {rlnorm(n = n, meanlog = log(ecsmu), sdlog = log(ecssig))},
                 function(n) {rnorm(n = n, mean = aeromu, sd = aerosig)},
                 function(n) {rnorm(n = n, mean = kappamu, sd = kappasig)})
    # Name the functions returned in the list by hector parameter
    names(list) <- c(ECS(), AERO_SCALE(), DIFFUSIVITY())

    list
}

# Randomly generate carbon cycle and climate parameter values based on the prior distribution. These are the
# parameters used in the emission driven runs. Will be used in the hector_sample function.
#
# Args
#   ecsmu: the climate sensitivity mean, default is set to the value in the hyper parameter defined in section 0
#   ecssig: the climate sensitivity  standard deviation
#   aeromu: the aerosol scalar mean
#   aerosig: the aerosol standard deviation
#   c0mu: mean preindustrial co2
#   c0sig: sd preindustrial co2
#   q10mu: mean temperature effect on heterotrophic respiration
#   q10sig: sd effect on heterotrophic respiration
#   betamu: mean co2 fertilization
#   betasig: sd co2 fertilization
# Returns: a list of functions that will sample the climate sensitivity, aero scaler, ocean diffusivity
# beta, q10, and pre industrial CO2. Each function in the list will sample the prior distribution n times.
make_emiss_dist <- function(ecsmu = hyper_params$ecsmu, ecssig = hyper_params$ecssig,
                            aeromu = hyper_params$aeromu, aerosig = hyper_params$ecssig,
                            kappamu = hyper_params$kappamu, kappasig = hyper_params$kappasig,
                            c0mu = hyper_params$c0mu,  c0sig = hyper_params$c0sig,
                            q10mu = hyper_params$q10mu, q10sig = hyper_params$q10sig,
                            betamu = hyper_params$betamu, betasig = hyper_params$betasig){


    concentration_dists <- make_concen_dist(ecsmu = ecsmu, ecssig = ecssig,
                                            aeromu = hyper_params$aeromu, aerosig = hyper_params$ecssig,
                                            kappamu = hyper_params$kappamu, kappasig = hyper_params$kappasig)

    emission_dists <- append(concentration_dists, list(function(n){ rnorm(n = n, mean = c0mu, sd = c0sig) },
                                                       function(n){ hectorcal::rtruncnorm(n = n, a = 0, b = Inf, mu = q10mu, sig = q10sig) },
                                                       function(n){ hectorcal::rtruncnorm(n = n, a = 0, b = Inf, mu = betamu, sig = betasig)}))
    # Name the functions returned in the list by hector parameter
    names(emission_dists) <- c(ECS(), AERO_SCALE(), DIFFUSIVITY(), PREINDUSTRIAL_CO2(), Q10_RH(), BETA())

    emission_dists
}

# Create a list n long with a hector core for every inifile
#
# Args
#   inifile: a vector of the paths to the hector ini files to use set up the cores in the each element of the list
#   names: a vector of the core names, should be equal in length to the inifile argument
#   n: the length of the list of hector cores to return, this should correspond to the number of nodes
#       the hector runs are going to parallelized  over
# Returns: the a list of hector cores to parallelize over
setup_cores <- function(inifile,name, n) {
    lapply(1:n, function(i) {newcore(inifile, name = name, suppresslogging=TRUE)})
}

# run Hector once using a sampled parameter set, this function is meant to be used
#
# Args
#   core: a Hector core to use
#   idx: an index of the parameter values to use, determined by the set up in the hector_sample function
#   param_vals: a list of the parameter values to run
#   keeptimes: a vector of the years to return the Hector output
#   keepvars: the paramter values
#   runid: a vector of run id
run_sample <- function(core, idx, param_vals, keeptimes, keepvars, runid)
{
    for(param in names(param_vals)) {

        # Figure out the paramter units
        units_list        <- c('degC', '(unitless)', 'cm2/s', 'ppmv CO2', '(unitless)', '(unitless)')
        names(units_list) <- c(ECS(), AERO_SCALE(), DIFFUSIVITY(), PREINDUSTRIAL_CO2(), Q10_RH(), BETA())
        unit_indx <- which(names(units_list) == param)

        # Set the Hector paramter values
        setvar(core, NA, param, param_vals[[param]][idx], units_list[[unit_indx]])
    }
    reset(core)
    stat <- tryCatch(
        run(core, max(keeptimes)),
        error = function(e){NULL})

    if(is.null(stat)) {
        ## Hector run failed, probably due to excessively high or low temperature.
        data.frame(runid = runid[idx],
                   value = NA,
                   variable = NA,
                   scenario = NA,
                   year = NA)

    }
    else {
        fetchvars(core, keeptimes, keepvars) %>%
            mutate(runid = runid[idx]) %>%
            select(runid, value, variable, scenario, year)



    }
}



# Run hector multiple times with sampled parameter values
#
# Args
#   n: the number of parameter samples to use as Hector inputs
#   hcores: the number of hector cores to use/parallize the runs over
#   keeptime: the Hector output result years to keep
#   keepvars: a vector of the Hector variable names to return output for, default is set to hector::GLOBAL_TEMP
#   priorDist: a list of the functions that will sample the prior distribution n times
hectorSample <- function(n, hcores, keeptime =1850:2100, keepvars =c(GLOBAL_TEMP()),
                         priorDist = NULL, param_vals = NULL) {

    # priorDist and param_vals cannot both be NULL
    if(is.null(param_vals) & !is.null(priorDist)){

        # Generate matrix of outputs for hector runs sampled from the supplied distribution
        param_vals <- lapply(priorDist, function(f){f(n)})
    } else if (is.null(priorDist) & !is.null(param_vals)) {
        # Do nothing because the param_vals being passed into the function is being use
    } else {
        # However if the conditions of the if else if statement are not met then
        # throw some error
        stop('problem with priorDist and param_vals')
    }

    # Organize the runs into batches that will be run in parallel
    max_cores <- parallel::detectCores()
    ncore     <- length(hcores)
    if(ncore >= max_cores) stop('hcores cannot exceed ', max_cores)
    nbatch <- as.integer(floor(n/ncore))
    nextra <- as.integer(n%%ncore)

    # Register the cores for parallelizing the hector runs
    doParallel::registerDoParallel(cores=ncore)


    # Figure out the run names, this will be added to the Hector
    # results and the parameter values
    runid <- 1:n

    # For each of the complete batches parallelize the Hector runs
    batch_rslts <- foreach(k=1:nbatch, .combine=rbind) %do% {
        ## Run a parallel batch
        foreach(i=1:ncore, .combine=rbind) %dopar% {
            core <- hcores[[i]]
            idx <- (k-1)*ncore + i
            run_sample(core, idx, param_vals, keeptime, keepvars, runid)
        }
    }

    # If there are any runs that did not fit into the complete batches
    # run them now.
    if(nextra > 0) {
        extra_rslts <- foreach(i=1:nextra, .combine=rbind) %dopar% {

            core <- hcores[[i]]
            idx <- nbatch*ncore + i
            run_sample(core, idx, param_vals, keeptime, keepvars, runid)

        }

        # Combine the batch and the extra run results
        rslt <- rbind(batch_rslts, extra_rslts)
    } else {

        # Return only the batch results
        rslt <- batch_rslts
    }

    # Kill the parallel cluster
    doParallel::stopImplicitCluster()

    # Format the results
    row.names(rslt) <- NULL

    # Determine the status of the run, if it failed or not
    param_df     <- data.frame(runid = runid , do.call(cbind, param_vals))

    dplyr::left_join(rslt, param_df, by = 'runid')
}


# 2. Runs ------------------------------------------------------------------------------------------------
# Make the output directory
output_dir <- file.path(pic_dir, 'PCA_pic_results')
dir.create(output_dir)

# Make the concentration and the emission prior distribution functions.
concentration_dists <- make_concen_dist()
emission_dists      <- make_emiss_dist()

# Find the ini files and set up the ini vectors to used in the runs
all_inifiles <- list.files(system.file('input', package = 'hector'), full.names = TRUE)
emission_ini <- all_inifiles[grepl(pattern = 'rcp[0-9]{2}.ini', all_inifiles)]
concen_ini   <- all_inifiles[grepl(pattern = 'rcp[0-9]{2}_constrained', all_inifiles)]


# Emission Driven Runs with the CC parameters varied -----
# RCP 26
hcores    <- setup_cores(emission_ini[grepl('rcp26', emission_ini)], n = nodes, name ='emissCC-esmrcp26')
rcp26_emissionCC <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), priorDist = emission_dists)
saveRDS(object = rcp26_emissionCC, file = file.path(output_dir, 'emissCC-RCP26.rds'))

# Format the parameters back into a list to reuse in the hectorSample function
emission_params <- list('S' = unique(rcp26_emissionCC$S),
                        'alpha' = unique(rcp26_emissionCC$alpha),
                        'diff' = unique(rcp26_emissionCC$diff),
                        'C0' = unique(rcp26_emissionCC$C0),
                        'q10_rh' = unique(rcp26_emissionCC$q10_rh),
                        'beta' = unique(rcp26_emissionCC$beta))

# RCP 45
hcores    <- setup_cores(emission_ini[grepl('rcp45', emission_ini)], n = nodes, name ='emissCC-esmrcp45')
rcp45_emissionCC <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), param_vals = emission_params)
saveRDS(object = rcp45_emissionCC, file = file.path(output_dir, 'emissCC-RCP45.rds'))

# RCP 60
hcores    <- setup_cores(emission_ini[grepl('rcp60', emission_ini)], n = nodes, name ='emissCC-esmrcp60')
rcp60_emissionCC <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), param_vals = emission_params)
saveRDS(object = rcp60_emissionCC, file = file.path(output_dir, 'emissCC-RCP60.rds'))

# RCP 85
hcores    <- setup_cores(emission_ini[grepl('rcp85', emission_ini)], n = nodes, name ='emissCC-esmrcp85')
rcp85_emissionCC <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), param_vals = emission_params)
saveRDS(object = rcp85_emissionCC, file = file.path(output_dir, 'emissCC-RCP85.rds'))



# Concentration Driven Runs -----
# RCP 26
hcores    <- setup_cores(concen_ini[grepl('rcp26', concen_ini)], n = nodes, name ='concen-rcp26')
rcp26_concen <- hectorSample(n_runs, hcores, keeptime = years, keepvars = GLOBAL_TEMP(), priorDist = concentration_dists)
saveRDS(object = rcp26_concen, file = file.path(output_dir, 'concen-RCP26.rds'))

# Format the parameters back into a list to reuse in the hectorSample function
concen_params <- list('S' = unique(rcp26_concen$S),
                      'alpha' = unique(rcp26_concen$alpha),
                      'diff' = unique(rcp26_concen$diff))

# RCP 45
hcores    <- setup_cores(concen_ini[grepl('rcp45', concen_ini)], n = nodes, name ='concen-rcp45')
rcp45_concen <- hectorSample(n_runs, hcores, keeptime = years, keepvars = GLOBAL_TEMP(), param_vals = concen_params)
saveRDS(object = rcp45_concen, file = file.path(output_dir, 'concen-RCP45.rds'))

# RCP 60
hcores    <- setup_cores(concen_ini[grepl('rcp60', concen_ini)], n = nodes, name ='concen-rcp60')
rcp60_concen <- hectorSample(n_runs, hcores, keeptime = years, keepvars = GLOBAL_TEMP(), param_vals = concen_params)
saveRDS(object = rcp60_concen, file = file.path(output_dir, 'concen-RCP60.rds'))

# RCP 85
hcores    <- setup_cores(concen_ini[grepl('rcp85', concen_ini)], n = nodes, name ='concen-rcp85')
rcp85_concen <- hectorSample(n_runs, hcores, keeptime = years, keepvars = GLOBAL_TEMP(), param_vals = concen_params)
saveRDS(object = rcp85_concen, file = file.path(output_dir, 'concen-RCP85.rds'))
# Emission Driven Runs with Constant C cycle -------

# Emission Driven Runs with the CC parameters varied -----

# Format the parameters back into a list to reuse in the hectorSample function
emission_params_Constant_carbon <- list('S' = unique(rcp26_emissionCC$S),
                                        'alpha' = unique(rcp26_emissionCC$alpha),
                                        'diff' = unique(rcp26_emissionCC$diff))


# RCP 26
hcores    <- setup_cores(emission_ini[grepl('rcp26', emission_ini)], n = nodes, name ='emissConstantC-esmrcp26')
time_test <- system.time(rcp26_emissionConstantC <- hectorSample(n_runs, hcores, keeptime = years, keepvars = GLOBAL_TEMP(), param_vals = emission_params_Constant_carbon))
saveRDS(object = rcp26_emissionConstantC, file = file.path(output_dir, 'emissConstantC-RCP26.rds'))


# RCP 45
hcores    <- setup_cores(emission_ini[grepl('rcp45', emission_ini)], n = nodes, name ='emissConstantC-esmrcp45')
rcp45_emissionConstantC <- hectorSample(n_runs, hcores, keeptime = years, keepvars = GLOBAL_TEMP(), param_vals = emission_params_Constant_carbon)
saveRDS(object = rcp45_emissionConstantC, file = file.path(output_dir, 'emissConstantC-RCP45.rds'))

# RCP 60
hcores    <- setup_cores(emission_ini[grepl('rcp60', emission_ini)], n = nodes, name ='emissConstantC-esmrcp60')
rcp60_emissionConstantC <- hectorSample(n_runs, hcores, keeptime = years, keepvars = GLOBAL_TEMP(), param_vals = emission_params_Constant_carbon)
saveRDS(object = rcp60_emissionConstantC, file = file.path(output_dir, 'emissConstantC-RCP60.rds'))

# RCP 85
hcores    <- setup_cores(emission_ini[grepl('rcp85', emission_ini)], n = nodes, name ='emissConstantC-esmrcp85')
rcp85_emissionConstantC <- hectorSample(n_runs, hcores, keeptime = years, keepvars = GLOBAL_TEMP(), param_vals = emission_params_Constant_carbon)
saveRDS(object = rcp85_emissionConstantC, file = file.path(output_dir, 'emissConstantC-RCP85.rds'))


Sys.time()
# End
