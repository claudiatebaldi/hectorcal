# This script creates an ensemble of Hector runs that will be for the PCA.
# Right now this is set up to run on pic.
#
# TODO a possible enchancemnet would be to make the number of cores in each element of the
# hector core list be flexibe instead of 3 but when I tried to use an lapply of for each
# statement to itterate over the hector list the computation time increased too much.
# be felxible instead of must each, this script returns a better data scrture with less calls
# but is much slower

# 0. Set Up --------------------------------------------------------------------------------
Sys.time()

pic_dir <- '.'

# Load the required pacakges
library(hector)
devtools::load_all('../')
library(doParallel)
library(dplyr)


# Define the prior hyperparameters, these should be the same as the hyperparameters defined in
# the liklihood function. These will be used as default inputs into functions defined
# in section 1 to create the pirior distirbtuion.
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


# Set the seed to ensure reproduceable results
set.seed(867-5309)

# This should be equal to the number of nodes selected in the
nodes <- 20

# The number of Hector runs per rcp scenario
n_runs <- 1200

# Years to extract the Hector output for
years  <- 1950:2100

# 1. Define functions ---------------------------------------------------------------------

# Randomly generate carbon values based on the prior distribtuion. These are the
# parameters used in the emission driven runs. Will be used in the hector_sample function.
#
# Args
#   ecsmu: the climate sensitivity mean, default is set to the valuein the hyper paramter defined in section 0
#   ecssig: the climate sensitivity  standard deviation
#   aeromu: the aerosol scalar mean
#   aerosig: the aerosol standard deviation
# Returns: a list of functions that will sample the climate senstivity, aero scaler, and ocean diffustvity
# parameter spaces. Each function in the list will sample the prior distrubtion n times.
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

# Randomly generate carbon cycle and climate parameter values based on the prior distribtuion. These are the
# parameters used in the emission driven runs. Will be used in the hector_sample function.
#
# Args
#   ecsmu: the climate sensitivity mean, default is set to the valuein the hyper paramter defined in section 0
#   ecssig: the climate sensitivity  standard deviation
#   aeromu: the aerosol scalar mean
#   aerosig: the aerosol standard deviation
#   c0mu: mean preindustiral co2
#   c0sig: sd peridustrial co2
#   q10mu: mean temperature effect on heterorophic respieration
#   q10sig: sd effect on heterorophic respiration
#   betamu: mean co2 fertlization
#   betasig: sd co2 fertlization
# Returns: a list of functions that will sample the climate senstivity, aero scaler, ocean diffustvity
# beta, q10, and preindustrical CO2. Each function in the list will sample the prior distrubtion n times.
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
#       the hector runs are going to parallized over
# Returns: the a list of hector cores to parallize over
setup_cores <- function(inifile, names, n) {

    stopifnot(length(inifile) == length(names))

    lapply(1:n, function(i) {
        lapply(1:length(inifile), function(k){
            newcore(inifile = inifile[k], suppresslogging = FALSE, name = names[k])
        })
    })
}


# run Hector once using a sampled paramter set
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
        out <- rep(NA, length(keeptimes)*length(keepvars) * length(core))

    }
    else {
        # Pull out the Hector results and format into a row of data
        rslt <- fetchvars(core, dates = keeptimes, vars = keepvars) %>%
            mutate(col_name = paste0('X', year, '_', scenario, '_', variable))

        out        <- rslt$value
        names(out) <- rslt$col_name
        out
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
hectorSample <- function(n, hcores, keeptime =1850:2100, keepvars =c(GLOBAL_TEMP()), priorDist) {

    # TODO change this requirement
    if(length(hcores[[1]]) != 4) stop('each element in hecores must equal length of 4')

    # Generate the paramter values
    param_vals <- lapply(priorDist, function(f){f(n)})


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
    runid <- as.integer(1:n)

    # For each of the compelte batches parallelize the Hector runs
    batch_rslts <- foreach(k=1:nbatch, .combine=rbind) %do% {
        ## Run a parallel batch
        foreach(i=1:ncore, .combine=rbind) %dopar% {
            core <- hcores[[i]]
            idx <- (k-1)*ncore + i

            # TODO make this part flexible
            c('runid' = runid[idx],
              c(run_sample(core[[1]], idx, param_vals, keeptime, keepvars, runid),
                run_sample(core[[2]], idx, param_vals, keeptime, keepvars, runid),
                run_sample(core[[3]], idx, param_vals, keeptime, keepvars, runid),
                run_sample(core[[4]], idx, param_vals, keeptime, keepvars, runid)))


        } # close dopar
    } # close do

    # If there are any runs that did not fit into the complete batches
    # run them now.
    if(nextra > 0) {
        extra_rslts <- foreach(i=1:nextra, .combine=rbind) %dopar% {
            core <- hcores[[i]]
            idx <- nbatch*ncore + i

            # TODO make this part flexible
            c('runid' = runid[idx],
              c(run_sample(core[[1]], idx, param_vals, keeptime, keepvars, runid),
                run_sample(core[[2]], idx, param_vals, keeptime, keepvars, runid),
                run_sample(core[[3]], idx, param_vals, keeptime, keepvars, runid),
                run_sample(core[[4]], idx, param_vals, keeptime, keepvars, runid)))


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
    NA_rows <- unique(which(is.na(rslt), arr.ind = TRUE)[ , 1])
    out     <- rslt[-NA_rows, ]

    # Determine the status of the run, if it failed or not
    param_df     <- data.frame(runid = runid , do.call(cbind, param_vals))
    run_complete          <- rep(TRUE, nrow(param_df))
    run_complete[NA_rows] <- FALSE
    param_df$run_complete <- run_complete

    list(rslt=rslt, params=param_df)

}


# 2. Runs ------------------------------------------------------------------------------------------------
# Make the output directory
output_dir <- file.path(pic_dir, 'pic_results')
dir.create(output_dir)

# Make the concentration and the emission prior distriubtuion functions.
concentration_dists <- make_concen_dist()
emission_dists      <- make_emiss_dist()

# Find the ini files and set up the ini vectors to used in the runs
all_inifiles <- list.files(system.file('input', package = 'hector'), full.names = TRUE)
emission_ini <- all_inifiles[grepl(pattern = 'rcp[0-9]{2}.ini', all_inifiles)]
concen_ini   <- all_inifiles[grepl(pattern = 'rcp[0-9]{2}_constrained', all_inifiles)]

# The emission driven runs with the carbon cycle paramters varying
emission_cores      <- setup_cores(emission_ini,
                                   names = c('emissRCP26', 'emissRCP45', 'emissRCP60', 'emiss85'), nodes)
message('starting the emissions cc run')
system.time(emission_results_CC <- hectorSample(n_runs, emission_cores, keeptime = years,
                                                keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), priorDist = emission_dists))
message('the emissions cc run complete')
saveRDS(emission_results_CC, file = file.path(output_dir, 'emission_CCparams.rds'))

# The emission driven runs with only the climate paramters varying, to keep the carbon paramters
# constant use the concentration distribution functions.
emission_results <- hectorSample(n_runs, emission_cores, keeptime = years,
                                 keepvars = c(GLOBAL_TEMP(), ATMOSPHERIC_C()), priorDist = concentration_dists)
saveRDS(emission_results, file = file.path(output_dir, 'emission.rds'))

# The concentration driven runs
concen_cores <- setup_cores(concen_ini, names = c('concenRCP26', 'concenRCP45', 'concenRCP60', 'concen85'), nodes)
concen_results <- hectorSample(n_runs, concen_cores, keeptime = years, keepvars = GLOBAL_TEMP(), priorDist = concentration_dists)
saveRDS(concen_results, file = file.path(output_dir, 'concentration.rds'))


Sys.time()



