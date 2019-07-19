## Hector Calibration Paper 1
## This script produces the figures and the summary results for the first Hector calibration paper that
## looks uses a nonlinear optimization routine to tune Hector to act as a CMIP5 emulator.
## Working directory should be set to the analysis/paper_1 directory.
## 0. Set Up -----------------------------------------------------------------------------------------------------
# The required packages
library(ggplot2);  library(ggthemes)
library(dplyr);    library(tidyr)
library(hector);   library(hectorcal)
library(cowplot);  library(knitr)
library(kableExtra)
library(assertthat)

# Define the directories
BASE_DIR   <- getwd()
OUTPUT_DIR <- file.path(BASE_DIR, 'output')
FIGS_DIR   <- file.path(BASE_DIR, 'figs')

# Define the default script plotting options.
COLOR_THEME  <- c('cmip' = 'grey', 'hector' = "#56B4E9")
SCRIPT_THEME <- theme_bw(base_size = 10) +
    theme(legend.title = element_blank())
FIG_WIDTH  <- 8   #in
FIG_DPI    <- 300
FIG_RATIO  <- 3


#################################################################################################################
### A. Concentration Results
#################################################################################################################
## 1. Import Fits -----------------------------------------------------------------------------------------------------

# Import the fits from the temperature-heat flux calibration experiment.
list.files(file.path(OUTPUT_DIR, 'conc_temp_heatflux'), 'fit_summary_table.csv', full.names = TRUE) %>%
    read.csv(stringsAsFactors = FALSE) %>%
    mutate(comp_data = 'temp heatflux',
           method = 'temp heatflux') ->
    fits_heatflux

# Import the fits from the temperature-heat flux range calibration experiment but only for the models that
# were missing heat flux data.
list.files(file.path(OUTPUT_DIR, 'conc_temp_heatfluxRange'), 'fit_summary_table.csv', full.names = TRUE) %>%
    read.csv(stringsAsFactors = FALSE) %>%
    filter(!model %in% fits_heatflux$model) %>%
    mutate(comp_data = 'temp heatflux Range',
           method = 'temp heatflux') ->
    fits_heatfluxRange

# Import the fits from the temperature-only calibration experiment.
list.files(file.path(OUTPUT_DIR, 'conc_temp'), 'fit_summary_table.csv', full.names = TRUE) %>%
    read.csv(stringsAsFactors = FALSE) %>%
    filter(model %in% c(fits_heatflux$model, fits_heatfluxRange$model)) %>%
    mutate(comp_data = 'temp',
           method = 'temp') ->
    fits_temp

# Check to make sure that there are no dpublicates in the temperature heat flux runs
assert_that(length(c(intersect(fits_heatfluxRange$model, fits_heatflux$model),
                     intersect(fits_heatflux$model, fits_heatfluxRange$model))) == 0)

# I am not sure why this test is failing but it should not be right now.
assert_that(length(c(setdiff(c(fits_heatfluxRange$model, fits_heatflux$model), fits_temp$model),
                     setdiff(fits_temp$model, c(fits_heatfluxRange$model, fits_heatflux$model)))) == 0)

# Save the data in all of fitted results in a single data frame.
fits <- bind_rows(fits_temp, fits_heatflux, fits_heatfluxRange) %>%
    filter(!model %in% c('MRI-ESM1', 'CNRM-CM5-2'))


## 2. Summary Results -------------------------------------------------------------------------------------------
# Temperature-only Calibration Results
# Figure out how many fall into each category, the sum of the temp_hf and the temp_hfRange
# should equal the count of temp
fits %>%
    group_by(comp_data) %>%
    summarise(count = n_distinct(model))

# Determine the RMSE for the temperature-only calibration fits
fits %>%
    filter(comp_data == 'temp') %>%
    mutate(min_value = min_value^(1/2)) %>%
    pull(min_value) ->
    MSE_temp

mean(MSE_temp)
sd(MSE_temp)

# Summary stats of the temperature only constrained calibration parameter values.
summary(fits_temp$S)
summary(fits_temp$diff)
summary(fits_temp$alpha)
summary(fits_temp$volscl)


# Temperature-heat flux calibration results
# Determine the RMSE for the temperature - heat flux calibration fits
fits %>%
    filter(comp_data != 'temp') %>%
    mutate(min_value = min_value^(1/2)) %>%
    pull(min_value) ->
    MSE_TH

summary(MSE_TH)
sd(MSE_TH)

fits %>%
    filter(comp_data != 'temp') ->
    temp_hf_fits

summary(temp_hf_fits$S)
summary(temp_hf_fits$diff)
summary(temp_hf_fits$alpha)
summary(temp_hf_fits$volscl)

## 3. Parameter Dot Plots -------------------------------------------------------------------------------------------
# Make the multi panneled dot plots of the different parameter fits from the concentration driven experiments.

# Temperature - only calibration fits
fits %>%
    filter(S <= 30 & method == 'temp') %>%
    ggplot(aes(S)) +
    geom_dotplot() +
    labs(x = expression(~italic(S)~' ('~degree~'C)'),
         y = 'Frequency') +
    SCRIPT_THEME ->
    S_temp_plot

fits %>%
    filter(diff <= 30  & method == 'temp') %>%
    ggplot(aes(diff)) +
    geom_dotplot() +
    labs(x = expression(kappa~' ('~cm^2/'m)'),
         y = 'Frequency') +
    SCRIPT_THEME ->
    diff_temp_plot

fits %>%
    filter(method == 'temp') %>%
    ggplot(aes(alpha)) +
    geom_dotplot() +
    labs(x = expression(alpha[a]),
         y = 'Frequency') +
    SCRIPT_THEME ->
    alpha_temp_plot

fits %>%
    filter(method == 'temp') %>%
    ggplot(aes(volscl)) +
    geom_dotplot() +
    labs(x = expression(alpha[v]),
         y = 'Frequency') +
    SCRIPT_THEME ->
    volscl_temp_plot

# Format the four climate paramter dot pltos into a single plot with mulitple pannels.
temp_param_dotplot <- plot_grid(alpha_temp_plot, volscl_temp_plot, S_temp_plot, diff_temp_plot,
                              labels = c('A', 'B', 'C', 'D'),
                              label_size = 10)

# Save the figure
ggsave(temp_param_dotplot,
       filename = file.path(FIGS_DIR, 'conc_temp_param_dotplots.pdf'),
       device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO * 2,
       dpi = FIG_DPI * 2, units = 'in')


# Temperature-heat flux dot plots
fits %>%
    filter(method != 'temp') %>%
    ggplot(aes(S)) +
    geom_dotplot() +
    labs(x = expression(~italic(S)~' ('~degree~'C)'),
         y = 'Frequency') +
    SCRIPT_THEME ->
    S_tempHF_plot

fits %>%
    filter( method != 'temp') %>%
    ggplot(aes(diff)) +
    geom_dotplot() +
    labs(x = expression(kappa~' ('~cm^2/'m)'),
         y = 'Frequency') +
    SCRIPT_THEME ->
    diff_tempHF_plot

fits %>%
    filter(method != 'temp') %>%
    ggplot(aes(alpha)) +
    geom_dotplot() +
    labs(x = expression(alpha[a]),
         y = 'Frequency') +
    SCRIPT_THEME ->
    alpha_tempHF_plot

fits %>%
    filter(method != 'temp') %>%
    ggplot(aes(volscl)) +
    geom_dotplot() +
    labs(x = expression(alpha[v]),
         y = 'Frequency') +
    SCRIPT_THEME ->
    volscl_tempHF_plot


tempHF_param_plots <- plot_grid(alpha_tempHF_plot,
                              volscl_tempHF_plot,
                              S_tempHF_plot,
                              diff_tempHF_plot, labels = c('A', 'B', 'C', 'D'),
                              label_size = 10)

ggsave(tempHF_param_plots,
       filename = file.path(FIGS_DIR, 'conc_tempHF_param_dotplots.pdf'),
       device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO * 2,
       dpi = FIG_DPI * 2, units = 'in')


## 4. Fitted Parameter Table -------------------------------------------------------------------------------------------
# Create the bones for the table of the fitted paratmers for latex, some additional editing will have to be done in
# latex but this makes something that can be copied and pasted into the LaTex document.
fits %>%
    gather(param, value, S, diff, alpha, volscl) %>%
    select(model, method, param, value) %>%
    mutate(param = if_else(param == 'alpha', '1', param),
           param = if_else(param == 'volscl', '2', param),
           param = if_else(param == 'S', '3', param),
           param = if_else(param == 'diff', '4', param),
           param = as.integer(param)) %>%
    mutate(method = if_else(method == 'temp', 1, 2)) %>%
    mutate(param = paste0(method, param)) %>%
    select(model, param, value) %>%
    spread(param, value) %>%
    na.omit %>%
    kable(format = 'latex', digits = 2, align = 'c',
          col.names = c('Model',  'alpha', 'volscl', 'S', 'diff',
                        'alpha', 'volscl', 'S', 'diff')) %>%
    kable_styling(c("striped"), full_width = F) %>%
    add_header_above(c(" " = 1, "Temperature" = 4, "Temperature & Heat Flux" = 4)) %>%
    save_kable(file = "conc_fitted_params.pdf")

## 5. Selected Fits Plot ------------------------------------------------------------------------------------------------
# Select three fits and plot them against the ESM comparison data. This will be used to illsutrate Hector's performance
# as an emulator.

# Select the min, median, and max optmized fits, these correspond to the best, worst, and median emulator fits.
fits %>%
    filter(comp_data == 'temp') %>%
    mutate(dif_mean = abs(min_value - mean(min_value))) %>%
    filter(min_value == min(min_value) | min_value == max(min_value) | dif_mean == min(dif_mean)) ->
    extreme_temp_fits

# Create a tibble of the ini files / scenario names to run.
tibble(ini = c(system.file('input/hector_rcp26_constrained.ini', package = 'hector'),
               system.file('input/hector_rcp45_constrained.ini', package = 'hector'),
               system.file('input/hector_rcp60_constrained.ini', package = 'hector'),
               system.file('input/hector_rcp85_constrained.ini', package = 'hector')),
       name = c('rcp26', 'rcp45', 'rcp60', 'rcp85')) ->
    hector_ini_tib

# Create a list of the Hector cores that will be run using the temperature-only calibration fits
# selected above.
core_list <- mapply(newcore, inifile = hector_ini_tib$ini, name = hector_ini_tib$name)

# Run Hector for the historic and future scenarios and join with comparison data.
lapply(X = split(extreme_temp_fits, extreme_temp_fits$model), FUN = function(input = X){

    # Select the paramters out of the input vector.
    param_names   <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE())
    params        <- input[names(input) %in% param_names]
    names(params) <- param_names

    # Parmeterize and run the Hector cores.
    lapply(core_list, FUN = parameterize_core, params = params)
    lapply(core_list, reset)
    lapply(core_list, run)

    # Extract the temperature results and add the other input information to
    # the returned data frame.
    lapply(core_list, fetchvars, dates = 1850:2100, vars = GLOBAL_TEMP()) %>%
        bind_rows() %>%
        mutate(join = 1) %>%
        left_join(input %>% mutate(join = 1), by = 'join') %>%
        select(-join)

}) %>%
    # Format the output into a single data frame with the appropriate years for the  different experiments.
    bind_rows() %>%
    rename(experiment = scenario) %>%
    mutate(experiment = if_else(year <= 2005, 'historical', experiment)) %>%
    distinct ->
    hector_temp_output

# Add the cmip comparison data to the temperature data frame.
hector_temp_output %>%
    rename(hector = value) %>%
    mutate(variable = 'tas') %>%
    inner_join(cmip_individual %>%
                   select(year, model, ensemble, variable, experiment, cmip = value)) %>%
    gather(source, value, hector, cmip) ->
    hector_cmip

ggplot(data = hector_cmip) +
    geom_line(data = hector_cmip %>% filter(source == 'cmip'),
              aes(year, value, color = source, group = interaction(model, experiment, ensemble, source))) +
    geom_line(data = hector_cmip %>% filter(source == 'hector'),
              aes(year, value, color = source, group = interaction(model, experiment, ensemble, source))) +
    facet_wrap('model') +
    scale_color_manual(values = COLOR_THEME) +
    SCRIPT_THEME +
    labs(x = 'Year',
         y = expression('Temperature Anomaly '~degree~'C')) +
    ggsave(filename = file.path(FIGS_DIR, 'conc_selected_temp_comparison.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO,
           dpi = FIG_DPI, units = 'in')

## 6. How Did the Ocean Heat Flux Constraint Change the Temperature Output? -----------------------------------------------
# When we compare the Hector output from the temperature-only and temperature heat flux calibration exercises how much of
# a difference does the calibration method make on the temperature values returned? We would expect the temperature
# values to be fairly similar to one another despite having very different parameterizations.

tibble(ini = c(system.file('input/hector_rcp26_constrained.ini', package = 'hector'),
               system.file('input/hector_rcp26_constrained.ini', package = 'hector'),
               system.file('input/hector_rcp45_constrained.ini', package = 'hector'),
               system.file('input/hector_rcp60_constrained.ini', package = 'hector'),
               system.file('input/hector_rcp85_constrained.ini', package = 'hector')),
       name = c('historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85')) ->
    hector_ini_tib

core_list <- mapply(newcore, inifile = hector_ini_tib$ini, name = hector_ini_tib$name)

# Run Hector for all of the fits but make sure that when we run the Hector core that the
# we use a try catch so that if an error occurs it does not break everything.
lapply(X = split(fits, interaction(fits$model, fits$method, drop = TRUE)), FUN = function(input = X){

    # Select the paramters out of the input vector.
    param_names   <- c(ECS(), DIFFUSIVITY(), AERO_SCALE(), VOLCANIC_SCALE())
    params        <- input[names(input) %in% param_names]
    names(params) <- param_names
print(input)
print('1')
    lapply(core_list, FUN = parameterize_core, params = params)
    print('2')
    lapply(core_list, reset)

    tryCatch({

        lapply(core_list, run)
        lapply(core_list, fetchvars, dates = 1850:2100, vars = GLOBAL_TEMP()) %>%
            bind_rows() %>%
            mutate(join = 1) %>%
            left_join(input %>% mutate(join = 1), by = 'join') %>%
            select(-join)
    }, error = function(e){NULL})


}) %>%
    bind_rows() ->
    hector_rslts

hector_rslts %>%
    rename(experiment = scenario) %>%
    filter(c(experiment == 'hisotrical' & year <= 2005) | c(experiment != 'historical' & year > 2005)) %>%
    distinct ->
    hector_temp_output_all_fits


# The mean squared error between the hector results, I am not sure if this is actually the metric  we
# will want to report or if we  want it to be something else.
hector_temp_output_all_fits %>%
    select(experiment, year, variable, value, model, method) %>%
    spread(method, value) %>%
    na.omit() ->
    Hector_rslts_temp_tempHF

Hector_rslts_temp_tempHF %>%
    mutate(dif = (`temp heatflux` - temp)^2) %>%
    pull(dif) %>%
    mean ->
    MSE

# The RMSE between the temp-only and the temp-heatflux calibrated Hector results
MSE^(1/2)

# The percent difference in the results
Hector_rslts_temp_tempHF %>%
    mutate(dif = 100 * (`temp heatflux` - temp)/ temp) ->
    percent_change

summary(percent_change$dif)

# The alternative is to look at how much the what the min optmiized value for the temperature is.
fits %>%
    filter(comp_data == 'temp heatflux') %>%
    pull(model) %>%
    paste(collapse = '|') ->
    heatflux_models

fits %>%
    filter(comp_data == 'temp heatflux Range') %>%
    pull(model) %>%
    paste(collapse = '|') ->
    heatfluxR_models

heatflux_models_paths  <- list.files(path = file.path(OUTPUT_DIR, 'conc_temp_heatflux'), pattern = heatflux_models, full.names = TRUE)
heatfluxR_models_paths <- list.files(path = file.path(OUTPUT_DIR, 'conc_temp_heatfluxRange'), pattern = heatfluxR_models, full.names = TRUE)

lapply(append(heatflux_models_paths, heatfluxR_models_paths), function(input){

    object <- readRDS(input)

    if(object$optim_rslt$convergence == 0){

     object$MSE$method <- 'temp heatflux'
     object$MSE

    }
    }) %>%
    bind_rows() ->
    heatflux_MSE

# Calculate the temperature MSE
heatflux_MSE %>%
    filter(variable == GLOBAL_TEMP()) %>%
    group_by(model, experiment, variable) %>%
    summarise(value = mean(value)) %>%
    ungroup %>%
    group_by(model) %>%
    summarise(min_temp = mean(value)) %>%
    ungroup %>%
    mutate(method = 'temp heatflux') ->
    minvalue_tempheatflux

# How does this compare to the MSE from the temperature-only fits?
fits %>%
    select(model, method, min_value) %>%
    filter(method == 'temp') %>%
    select(-method) %>%
    left_join(minvalue_tempheatflux) %>%
    mutate(dif = (min_value - min_temp))


## 7. Symmetry Between S and Kappa -------------------------------------------------------------------------------
# The reason why the single calibration results contain unreasonable fitted parameters is because of the symmetry
# that exists in the simplified climate system. Here we visualize results from the conc_symmetry_runs.R that
# illustrate how S and Kappa can trade off with one another during the temperature-only calibration experiment
# but that this relationship is limited when heat flux output is added as an additional constraint.

S_kappa_temp   <- read.csv(list.files(path = file.path(OUTPUT_DIR, 'conc_S_kappa_temp'), '.csv', full.names = TRUE), stringsAsFactors = FALSE)
S_kappa_tempHF <- read.csv(list.files(path = file.path(OUTPUT_DIR, 'conc_S_kappa_tempHF'), '.csv', full.names = TRUE), stringsAsFactors = FALSE)
S_kappa        <- bind_rows(S_kappa_temp, S_kappa_fits_tempHF)

S_kappa_temp %>%
    bind_rows(S_kappa_fits_tempHF) %>%
    mutate(comp_data = if_else(comp_data == 'temp', 'Temp-Only', 'Temp-Heat Flux')) %>%
    filter(comp_data =='Temp-Only') ->
    S_kappa

ggplot(data = S_kappa, aes(kappa, min, color = model, shape = comp_data)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression('Temp MSE ('~degree~'C'^2~')'),
         x = expression(kappa~'('~cm^2*s^-1~')')) +
    coord_cartesian(ylim = c(0, 0.4)) +
    guides(shape = FALSE) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color = 'black')) ->
    conc_S_kappa_sym_plot


ggplot(data = S_kappa, aes(kappa, S, color = model, shape = comp_data)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression(~italic(S)~' ('~degree~'C)'),
         x = expression(kappa~'('~cm^2*s^-1~')')) +
    guides(shape = FALSE) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color = 'black')) ->
    conc_S_kappa_trade_off


conc_S_kapp_sym_plot_final <- plot_grid(conc_S_kappa_sym_plot,
                                conc_S_kappa_trade_off, labels = c('A', 'B'),
                                label_size = 10)

ggsave(conc_S_kapp_sym_plot_final,
       filename = file.path(FIGS_DIR, 'conc_S_kappa_sym_plot.pdf'),
       device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO * 2,
       dpi = FIG_DPI * 2, units = 'in')

## 8. Log Mesa Plot ---------------------------------------------------------------------------------------------
# For the section of the paper that introduces the idea of the mesa function we would like to include a figure
# to help facilitate our explanation of what is going on.

## Function for generating a log-mesa function for use in plotting functions
logmesaplt <- function(a,b,sig) {
    function(x) {
        -log(mesa(x, a, b, sig))
    }
}

ggplot(data=data.frame(x=c(-1.5, 1.5)),
       aes(x)) + stat_function(fun=logmesaplt(-1,1,0.15), size=1.1) +
    ylab('Log(x)') +
    SCRIPT_THEME ->
    mesaplt

ggsave(mesaplt,
       filename = file.path(FIGS_DIR, 'log_mesa_plot.pdf'),
       device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO * 2,
       dpi = FIG_DPI * 2, units = 'in')
