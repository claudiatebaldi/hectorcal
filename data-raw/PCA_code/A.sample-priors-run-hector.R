# This script creates an ensemble of Hector runs that will be for the PCA.
# Right now this is set up to run on pic.

# 0. Set Up --------------------------------------------------------------------------------
Sys.time()

pic_dir <- '.'

# Load the required pacakges
library(hector)
library(hectorcal)
library(doParallel)


# Define the prior hyperparameters, these should be the same as the hyperparameters
# defined in the liklihood function. These will be used to inform the prior distribtuion
# sampled in this script.
ecsmu   <- 3.0; ecssig   <- 3.0
aeromu  <- 1.0; aerosig  <- 1.4
kappamu <- 2.3; kappasig <- 2.0
betamu  <- 0.3; betasig  <- 0.7
q10mu   <- 2.0; q10sig   <- 2.0
c0mu    <- 285; c0sig    <- 14.0

# Set the seed to ensure reproduceable results
set.seed(867-5309)

# This should be equal to the number of nodes selected in the
nodes <- 20

# The number of Hector runs
n_runs <- 5000

# Years to extract the Hector output for
years  <- 1860:2100

# 1. Define functions ---------------------------------------------------------------------
# TODO turn these into acutal functions...
# Randomly generate climate system parameter values based on the prior distribution. These are
# the parameters used in the concentration driven runs. Will be used in the hector sample function.
#
# Args
#   n: the number of observations to sample
#   * the hyperparamters for the prior distrubtions are defined in section 0.
# Returns: a list of functions that will sample the climate senstivity, aero scaler, and ocean diffustvity
# parameter spaces. Each function in the list will sample the prior distrubtion n times.
concentration_dists <- list(function(n) {rlnorm(n = n, meanlog = log(ecsmu), sdlog = log(ecssig))},
                            function(n) {rnorm(n = n, mean = aeromu, sd = aerosig)},
                            function(n) {rnorm(n = n, mean = kappamu, sd = kappasig)})
# Name the functions returned in the list by hector parameter
names(concentration_dists) <- c(ECS(), AERO_SCALE(), DIFFUSIVITY())

# Randomly generate carbon cycle parameter values based on the prior distribtuion. These are the
# parameters used in the emission driven runs. Will be used in the hector_sample function.
#
# Args
#   n: the number of observations to sample
#   * the hyperparamters for the prior distrubtions are defined in section 0.
# Returns: a list of functions that will sample the climate senstivity, aero scaler, ocean diffustvity
# beta, q10, and preindustrical CO2. Each function in the list will sample the prior distrubtion n times.
emission_dists <- append(concentration_dists, list(function(n){ rnorm(n = n, mean = c0mu, sd = c0sig) },
                                                   function(n){ hectorcal::rtruncnorm(n = n, a = 0, b = Inf, mu = q10mu, sig = q10sig) },
                                                   function(n){ hectorcal::rtruncnorm(n = n, a = 0, b = Inf, mu = betamu, sig = betasig)}))
# Name the functions returned in the list by hector parameter
names(emission_dists) <- c(ECS(), AERO_SCALE(), DIFFUSIVITY(), PREINDUSTRIAL_CO2(), Q10_RH(), BETA())

# Create list of hector cores
#
# Args
#   inifile: the path to the hector ini file to use in the hector core
#   n: the number of hector cores to return
# Returns: the a list of hector cores to parallize over
setup_cores <- function(inifile, n) {
    lapply(1:n, function(i) {newcore(inifile, suppresslogging=TRUE)})
}

# run Hector once using a sampled paramter set, this function is meant to be used
#
# Args
#   core: a Hector core to use
#   idx: an index of the paramter values to use, determined by the set up in the hector_sample function
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
        rep(NA, length(keeptimes)*length(keepvars))
    }
    else {
        rslt <- fetchvars(core, keeptimes, keepvars)
        cbind(run_id = runid[[idx]], year = rslt$year, variable = rslt$variable, value = rslt$value)
    }
}



# Run hector mulitple times with sampled paramter values
#
# Args
#   n: the number of parameter samples to use as Hector intpus
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
    # results and the paramter values
    runid <- 1:n

    # For each of the compelte batches parallelize the Hector runs
    batch_rslts <- foreach(k=1:nbatch, .combine=rbind) %do% {
        ## Run a parallel batch
            foreach(i=1:ncore, .combine=rbind) %dopar% {
                core <- hcores[[i]]
                idx <- (k-1)*ncore + i
                run_sample(core, idx, param_vals, keeptime, keepvars, runid)

                } # close dopar
        } # close do

    # If there are any runs that did not fit into the complete batches
    # run them now.
    if(nextra > 0) {
        extra_rslts <- foreach(i=1:nextra, .combine=rbind) %dopar% {
            core <- hcores[[i]]
            idx <- nbatch*ncore + i
            run_sample(core, idx, param_vals, keeptime, keepvars, runid)
        } # close dopar

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

    list(rslt=rslt, params=cbind(run_id = runid , do.call(cbind, param_vals)))
}


# 2. Runs ------------------------------------------------------------------------------------------------
# Make the output directory
output_dir <- file.path(pic_dir, 'pic_results')
dir.create(output_dir)

# RCP 85
# Right now lets start off with 5,000 emission driven runs for a single rcp
indir  <- system.file('input',package = 'hector')
hcores <- setup_cores(file.path(indir, 'hector_rcp85.ini'), n = nodes)
rslt1  <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), priorDist = concentration_dists)
saveRDS(object = rslt1, file = file.path(output_dir, 'emission_rcp85.rds'))

# Format the parmaters back into a list to reuse in the hectorSample function
emission_params <- list('S' = rslt1$params[,2],
                        'alpha' = rslt1$params[,3],
                        'diff' = rslt1$params[,4])
                       # 'C0' = rslt1$params[, 5],
                       # 'q10_rh' = rslt1$params[, 6],
                       # 'beta' = rslt1$params[, 7])
# RCP 60
indir  <- system.file('input',package = 'hector')
hcores <- setup_cores(file.path(indir, 'hector_rcp60.ini'), n = nodes)
rslt2  <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), param_vals = emission_params)
saveRDS(object = rslt2, file = file.path(output_dir, 'emission_rcp60.rds'))

# RCP 45
indir  <- system.file('input',package = 'hector')
hcores <- setup_cores(file.path(indir, 'hector_rcp45.ini'), n = nodes)
rslt2  <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), param_vals = emission_params)
saveRDS(object = rslt2, file = file.path(output_dir, 'emission_rcp45.rds'))

# RCP 26
indir  <- system.file('input',package = 'hector')
hcores <- setup_cores(file.path(indir, 'hector_rcp26.ini'), n = nodes)
rslt2  <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), param_vals = emission_params)
saveRDS(object = rslt2, file = file.path(output_dir, 'emission_rcp26.rds'))


stop('end')

# The concentration driven runs
# RCP 85
# Right now lets start off with 5,000 emission driven runs for a single rcp
indir  <- system.file('input',package = 'hector')
hcores <- setup_cores(file.path(indir, 'hector_rcp85_constrained.ini'), n = nodes)
rslt1  <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), priorDist = concentration_dists)
saveRDS(object = rslt1, file = file.path(output_dir, 'concen_rcp85.rds'))

# Format the parmaters back into a list to reuse in the hectorSample function
concen_params <- list('S' = rslt1$params[,2],
                        'alpha' = rslt1$params[,3],
                        'diff' = rslt1$params[,4])

# RCP 60
# Right now lets start off with 5,000 emission driven runs for a single rcp
indir  <- system.file('input',package = 'hector')
hcores <- setup_cores(file.path(indir, 'hector_rcp60_constrained.ini'), n = nodes)
rslt1  <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), priorDist = concentration_dists)
saveRDS(object = rslt1, file = file.path(output_dir, 'concen_rcp60.rds'))

# RCP 45
# Right now lets start off with 5,000 emission driven runs for a single rcp
indir  <- system.file('input',package = 'hector')
hcores <- setup_cores(file.path(indir, 'hector_rcp45_constrained.ini'), n = nodes)
rslt1  <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), priorDist = concentration_dists)
saveRDS(object = rslt1, file = file.path(output_dir, 'concen_rcp45.rds'))

# RCP 26
# Right now lets start off with 5,000 emission driven runs for a single rcp
indir  <- system.file('input',package = 'hector')
hcores <- setup_cores(file.path(indir, 'hector_rcp26_constrained.ini'), n = nodes)
rslt1  <- hectorSample(n_runs, hcores, keeptime = years, keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), priorDist = concentration_dists)
saveRDS(object = rslt1, file = file.path(output_dir, 'concen_rcp26.rds'))




Sys.time()


