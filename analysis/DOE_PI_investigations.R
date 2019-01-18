# MCMC analysis for the DOE PI meeting
# 2018-10-26 


# 0. Set Up ---------- 
#BASE <- '/pic/projects/GCAM/Dorheim/hector_calibration'; stopifnot(dir.exists(BASE))
#BASE <- getwd() # should be project location
#full_compdata <- read.csv(file.path(BASE, 'L1-output', 'esm_comparison_data.csv'), stringsAsFactors = FALSE)

# Source the functions 
#source(file.path(BASE, 'code', 'L02A.mcmc_functions.R'))

# This script is set up to run on pic, the first part should be done interactively as we try to figure out scale size ect. 
# After the metrop is set up correctly this section can be commented out. 

# Carbon cycle turned on. 
##full_compdata %>%  
##  filter(experiment == 'esmrcp85') -> 
##  comp_esmrcp85
comp_esmrcp85 <- filter(esm_comparison, experiment=='esmrcp85')


# Define the inifiles for hector and create a new core 
esmrcp85_ini <- system.file("input/hector_rcp85.ini", package = "hector")

esmrcp85_co2_cal_lpost <- make_log_post_function(esmrcp85_ini, comp_esmrcp85, carbon_cycle = TRUE, calibration = 'consensus', step = 19, showMessages = TRUE)
test1 <- metrop(obj = esmrcp85_co2_cal_lpost, initial = c(2, 1, 2.05, 270, 2.6, 1), nbatch = 1, scale = .5)
  

esmrcp85_co2_mean_lpost <- make_log_post_function(esmrcp85_ini, comp_esmrcp85, carbon_cycle = TRUE, calibration = 'mean', step = 19, showMessages = TRUE)
test2 <- metrop(obj = esmrcp85_co2_mean_lpost, initial = c(2, 1, 2.05, 270, 2.6, 0.5, 1, 5), nbatch = 50, scale = c(1, 1, 1, 1, 1, 0.03, 1, 1))
test2 <- metrop(test2, nbatch = 100, scale = 0.01)
