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
library(kableExtra); library(scales)
library(assertthat)
library(GGally)

# Define the directories
BASE_DIR   <- getwd()
INPUT_DIR  <- file.path(BASE_DIR, 'output')
FIGS_DIR   <- file.path(BASE_DIR, 'figs')

# Define the default script plotting options, based on the color blind palette with grey.
COLOR_THEME  <- c('CMIP' = 'grey',
                  'Hector' = "#56B4E9",
                  "#E69F00",
                  "#56B4E9",
                  "#009E73",
                  "#56B4E9",
                  "#0072B2",
                  "#D55E00",
                  "#CC79A7")
SCRIPT_THEME <- theme_bw(base_size = 10) +
    theme(legend.title = element_blank())
FIG_WIDTH  <- 8   #in
FIG_DPI    <- 300
FIG_RATIO  <- 3


### A. Import and Format Results ###############################################################
# Because the calibration runs themselves are time intensive the calibration runs are launched by the
# 1A. scripts. Import all of the fitted calibration results. And subset the calibration results for
# so that it contains the results we are interested in. For example for this paper we perfer the
# results from the temp and heat flux calibration runs over the fits from the temp and heat flux range
# results.
## Concentration Driven -------------------------------------------------------------------------------
read.csv(list.files(INPUT_DIR, 'calibration_results.csv', full.names = TRUE), stringsAsFactors = FALSE) %>%
    filter(grepl('conc_', model)) %>%
    mutate(comp_data = gsub(pattern = 'conc_', replacement = '', x = method)) %>%
    mutate(method = if_else(grepl('heat' , comp_data), 'temp heatflux', 'temp')) %>%
    split(.$comp_data) ->
    fit_list

# Determine which models are missing from the temp - heatflux calibration experiments, this is
# because the comparison data set must be missing specific heat flux data so we will want to
# use the temp - heatflux Range results.
temp_heatfxRange_models <- setdiff(fit_list$`temp-heatfluxRange`$model,
                                   fit_list$`temp-heatflux`$model)

fit_list$`temp-heatfluxRange` %>%
    filter(model %in% temp_heatfxRange_models) %>%
    bind_rows(fit_list$`temp-heatflux`) ->
    fits_temp_heatflux

# Determine if there are any missing model fits from the two different
# calibration methods. There should be none.
missing_fits <- c(setdiff(fit_list$temp$model, fits_temp_heatflux$model),
                  setdiff(fits_temp_heatflux$model, fit_list$temp$model))

stopifnot(length(missing_fits) == 0)

conc_fits <- bind_rows(fits_temp_heatflux, fit_list$temp)

## Subset the emulated results to reflect the final conc_fits
read.csv(list.files(INPUT_DIR, 'emulated_Hector_output.csv', full.names = TRUE),
         stringsAsFactors = FALSE) %>%
    mutate(comp_data = gsub(pattern = 'conc_', replacement = '', x = method)) %>%
    select(-method) %>%
    inner_join(conc_fits %>%
                   select(model, comp_data, method), by = c("model", "comp_data")) ->
    emulated_conc


# Import the concentration driven symmetry test results.
read.csv(list.files(INPUT_DIR, 'sym_results.csv', full.names = TRUE), stringsAsFactors = FALSE) %>%
    filter(driven == 'conc') ->
    conc_sym_rslts

## Emission Driven --------------------------------------------------------------------------------------
# Import emulated emission driven results
read.csv(list.files(INPUT_DIR, 'emulated_Hector_output.csv', full.names = TRUE), stringsAsFactors = FALSE) %>%
    filter(grepl(pattern = 'emiss_', x = model)) %>%
    mutate(model = gsub(pattern = 'emiss_', replacement = '', x = model),
           penalty = if_else(grepl(pattern = 'penalty', x = tolower(method)), 'with', 'without')) %>%
    mutate(variable = if_else(variable == 'Ca', 'co2', variable)) ->
    emulated_emiss

# Import the emission driven Hector fit results.
list.files(INPUT_DIR, 'calibration_results.csv', full.names = TRUE) %>%
    read.csv(stringsAsFactors = FALSE) %>%
    filter(grepl(pattern = 'emiss_', x = model)) %>%
    mutate(model = gsub(pattern = 'emiss_', replacement = '', x = model),
           penalty = if_else(grepl(pattern = 'penalty', x = tolower(method)), 'with', 'without')) ->
    emiss_fits

# Import the emission driven symetry runs
list.files(INPUT_DIR, "sym_results.csv", full.names = TRUE)  %>%
    read.csv(stringsAsFactors = FALSE) %>%
    filter(driven == 'emiss') %>%
    filter(!signif(beta, 7) %in% c(signif(0.6755556, 7), signif(0.0100000, 7))) %>%
    filter(beta <= 1) ->
    emiss_sym




### A. Concentration Results ###############################################################
## Select Example Models
example_models <- sort(c('CCSM4', 'GFDL-CM3', 'MRI-CGCM3'))
example_models_COLORS <- COLOR_THEME[3:7]
names(example_models_COLORS) <- example_models

## 1. Results Summary -------------------------------------------------------------------------------------------
# Temperature-only Calibration Results
# Figure out how many fall into each category, the sum of the temp_hf and the temp_hfRange
# should equal the count of temp
conc_fits %>%
    group_by(comp_data) %>%
    summarise(count = n_distinct(model))

# Determine the RMSE for the temperature-only calibration fits
conc_fits %>%
    filter(comp_data == 'temp') %>%
    mutate(RMSE = optmized_value^(1/2)) ->
    fits_temp

# RMSE of the temp-only fits
mean(fits_temp$RMSE)
sd(fits_temp$RMSE)

# Summary stats of the temperature only constrained calibration parameter values.
summary(fits_temp$S)
summary(fits_temp$diff)
summary(fits_temp$alpha)
summary(fits_temp$volscl)

fits_temp %>%
    mutate(good = FALSE,
           good = if_else(S >= 1 & S <= 6, TRUE, good)) %>%
   #  mutate(good = if_else(diff > 1e-6, TRUE, good)) %>%
    filter(diff > 1e-6) %>%
    filter(good) %>%
    group_by(good) %>%
    summarise(n = n_distinct(model))



# Summary about the fits for the temp-heatflux calibration experiments
conc_fits %>%
    filter(comp_data != 'temp') ->
    fits_tempHF

summary(fits_tempHF$S)
summary(fits_tempHF$diff)
summary(fits_tempHF$alpha)
summary(fits_tempHF$volscl)


## 2. Parameter Dot Plots -------------------------------------------------------------------------------------------
# Make the multi panneled dot plots of the different parameter fits from the concentration driven experiments.
mround <- function(x,base){
    base*round(x/base)
}
mround(14,5)


# Temperature - only calibration fits
conc_fits %>%
    filter(method == 'temp') %>%
    ggplot(aes(S)) +
    geom_rect(aes(xmin=1.5, xmax=4.5, ymin=0, ymax=Inf), fill = '#D3D3D3') +
    geom_dotplot() +
    labs(x = expression(~italic(S)~' ('~degree~'C)'),
         y = 'Frequency') +
    SCRIPT_THEME +
    scale_x_continuous(trans = 'log10',
                       breaks = trans_breaks("log10", function(x) signif(10^x, digits = 1))) +
    annotation_logticks(sides = 'b') ->
    S_temp_plot

conc_fits %>%
    filter(method == 'temp') %>%
    ggplot(aes(diff)) +
    geom_rect(aes(xmin=1e-6, xmax=Inf, ymin=0, ymax=Inf), fill = '#D3D3D3') +
    geom_dotplot() +
    labs(x = expression(kappa~' ('~cm^2/'s)'),
         y = 'Frequency') +
    SCRIPT_THEME +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x){10^x})) ->
    diff_temp_plot

conc_fits %>%
    filter(method == 'temp') %>%
    ggplot(aes(alpha)) +
    geom_dotplot() +
    labs(x = expression(alpha[a]),
         y = 'Frequency') +
    SCRIPT_THEME ->
    alpha_temp_plot

conc_fits %>%
    filter(method == 'temp') %>%
    ggplot(aes(volscl)) +
    geom_dotplot() +
    labs(x = expression(alpha[v]),
         y = 'Frequency') +
    SCRIPT_THEME ->
    volscl_temp_plot

# Format the four climate paramter dot pltos into a single plot with mulitple pannels.
temp_param_dotplot <- plot_grid(S_temp_plot, diff_temp_plot,
                                labels = c('A', 'B'),
                                label_size = 10)

# Save the figure
ggsave(temp_param_dotplot,
       filename = file.path(FIGS_DIR, 'conc_temp_param_dotplots.pdf'),
       device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO,
       dpi = FIG_DPI * 2, units = 'in')


# Temperature-heat flux dot plots
conc_fits %>%
    filter(method != 'temp') %>%
    ggplot(aes(S)) +
    geom_dotplot() +
    labs(x = expression(~italic(S)~' ('~degree~'C)'),
         y = 'Frequency') +
    SCRIPT_THEME ->
    S_tempHF_plot

conc_fits %>%
    filter( method != 'temp') %>%
    ggplot(aes(diff)) +
    geom_dotplot() +
    labs(x = expression(kappa~' ('~cm^2/'m)'),
         y = 'Frequency') +
    SCRIPT_THEME ->
    diff_tempHF_plot

conc_fits %>%
    filter(method != 'temp') %>%
    ggplot(aes(alpha)) +
    geom_dotplot() +
    labs(x = expression(alpha[a]),
         y = 'Frequency') +
    SCRIPT_THEME ->
    alpha_tempHF_plot

conc_fits %>%
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

## 3. Change in Parameter Plots --------------------------------------------------------------------
# For the change in parameter value plots use the colors from the script theme.
CHANGE_PLOT_COLORS <- c('temp' = "#D3D3D3",
                        'temp heatflux' = "#4169E1")


# Edit the model column to remove the conc_ string of the name.
conc_fits  <-  mutate(conc_fits, model = gsub(pattern = 'conc_', replacement = '', x = model))



# Add a factor to the fits data frame so that the models will be plotted in alphabetical order.
conc_fits$model <- factor(x = conc_fits$model, levels = sort(unique(conc_fits$model), decreasing = TRUE),
                     ordered = TRUE)

ggplot(data = conc_fits) +
    geom_point(aes(y = model, x = S, color = method),
               size = 2.5) +
    geom_path(aes(y = model, x = S),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'first')) +
    labs(y = NULL)  +
    coord_cartesian(xlim = c(0, 15)) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    labs(x = expression(~italic(S)~' ('~degree~'C)')) ->
    change_fitted_S

ggplot(data = conc_fits) +
    geom_point(aes(y = model, x = diff, color = method),
               size = 2.5) +
    geom_path(aes(y = model, x = diff),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'first')) +
    labs(y = NULL) +
    coord_cartesian(xlim = c(-1, 20)) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    labs(x = expression(kappa~' ('~cm^2/'s)')) ->
    change_fitted_diff

ggplot(data = conc_fits) +
    geom_point(aes(y = model, x = alpha, color = method),
               size = 2.5) +
    geom_path(aes(y = model, x = alpha),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'first')) +
    labs(y = NULL) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS)+
    coord_cartesian(xlim = c(0, 2.25)) +
    labs(x = expression(alpha[a])) ->
    change_fitted_alpha

ggplot(data = conc_fits) +
    geom_point(aes(y = model, x = volscl, color = method),
               size = 2.5) +
    geom_path(aes(y = model, x = volscl),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'first')) +
    labs(y = NULL) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    coord_cartesian(xlim = c(-.5, 2)) +
    labs(x = expression(alpha[v])) ->
    change_fitted_volscl


change_param_plot <- plot_grid(change_fitted_S + theme(legend.position = 'none'),
                               change_fitted_diff + theme(legend.position = 'none'),
                               change_fitted_alpha + theme(legend.position = 'none'),
                               change_fitted_volscl + theme(legend.position = 'none'),
                               labels = c('A', 'B', 'C', 'D'),
                               label_size = 10)

legend <- get_legend(change_fitted_S + theme(legend.box.margin = margin(0, 0, 0, 12)))

final_change_param_plot <- plot_grid(change_param_plot, legend, ncol = 1, rel_heights = c(4, .25))
# Save the figure
ggsave(final_change_param_plot,
       filename = file.path(FIGS_DIR, 'change_fitted_param.pdf'),
       device = 'pdf', width = FIG_WIDTH , height = FIG_WIDTH ,
       dpi = FIG_DPI * 3, units = 'in')


## 4. Fitted Parameter Table -------------------------------------------------------------------------------------------
# Create the bones for the table of the fitted paratmers for latex, some additional editing will have to be done in
# latex but this makes something that can be copied and pasted into the LaTex document.
conc_fits %>%
    gather(param, value, S, diff, alpha, volscl) %>%
    select(model, method, param, value) %>%
    mutate(model = gsub('conc_', '', model)) %>%
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
    add_header_above(c(" " = 1, "Temperature" = 4, "Temperature & Heat Flux" = 4))  %>%
    save_kable(file = file.path(FIGS_DIR, "conc_fitted_params.pdf"))

## 5. Selected Fits Plot ------------------------------------------------------------------------------------------------
# Subset the emulated result to only include data for the
emulated_conc %>%
    mutate(model = gsub(pattern = 'conc_', replacement = '', x = model)) %>%
    filter(model %in% example_models) %>%
    rename(Hector = value) %>%
    inner_join(cmip_individual %>%
                  select(year, model, variable, experiment, ensemble, CMIP = value),
              by = c("model", "experiment", "year", "variable")) %>%
    filter(variable == 'tas' & comp_data == 'temp') %>%
    gather(source, value, Hector, CMIP) ->
    hector_cmip

# Plot each of the Hector vs CMIP ESM output data separately.
split(hector_cmip, hector_cmip$model, drop = TRUE) %>%
    lapply(function(data){

        ggplot() +
            geom_line(data = data %>% filter(source == 'CMIP'),
                      aes(year, value, color = source, group = interaction(model, experiment, ensemble, source)),
                      size = 0.25,
                      alpha = 0.75) +
            geom_line(data = data %>% filter(source == 'Hector'),
                      aes(year, value, color = source, group = interaction(model, experiment, ensemble, source)),
                      size = 0.75) +
            scale_color_manual(values = COLOR_THEME) +
            SCRIPT_THEME +
            labs(x = 'Year',
                 y = expression('Temperature Anomaly '~degree~'C')) +
            geom_text(data = tibble(x=min(data$year),
                                    y = max(data$value),
                                    label = unique(data$model)),
                      aes(x, y, label = label, hjust = 'left'))
    }) ->
    plots


# Arrange the plots in a list in a grid based on the MSE value.
arranged_plot_list <- list()
for(i in 1:length(example_models)){

    name <- names(plots)[i]
    arranged_plot_list[[name]] <- plots[[name]] + theme(legend.position = 'none')

    if(i > 1){

        arranged_plot_list[[name]] <- arranged_plot_list[[name]] + labs(y = NULL)

    }


}
arranged_plots <- plot_grid(plotlist = arranged_plot_list, nrow = 1)

# Extract the legend from one of the plots, so that there is a universal plot.
legend <- get_legend(plots[[1]]+ theme(legend.box.margin = margin(0, 0, 0, 12)))

# Add the legend to the row we made earlier and save.
plot_grid(arranged_plots, legend, rel_widths = c(4, .5))+
    ggsave(filename = file.path(FIGS_DIR, 'conc_selected_temp_comparison.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO,
           dpi = FIG_DPI, units = 'in')


## 6. How Did the Calibration Method Change the Temperature Results? ---------------------------------------------
# Calculate the RMSE between the temp-only and temp-heatflux data sets.
emulated_conc %>%
    filter(variable == 'tas') %>%
    select(model, method, experiment, year, value) %>%
    distinct %>%
    spread(method, value) ->
    wide_emualted_conc

assert_that(all(!is.na(wide_emualted_conc)))

wide_emualted_conc %>%
    mutate(SE = (temp - `temp heatflux`) ^ 2) %>%
    group_by(model, experiment) %>%
    summarise(RMSE = mean(SE) ^ 1/2 )%>%
    ungroup %>%
    group_by(model) %>%
    summarise(modelRMSE = mean(RMSE)) %>%
    ungroup %>%
    pull(modelRMSE) %>%
    mean(na.rm = TRUE)



## 7. Symmetry Between S and Kappa -------------------------------------------------------------------------------
# The reason why the single calibration results contain unreasonable fitted parameters is because of the symmetry
# that exists in the simplified climate system. Here we visualize results from the conc_symmetry_runs.R that
# illustrate how S and Kappa can trade off with one another during the temperature-only calibration experiment
# but that this relationship is limited when heat flux output is added as an additional constraint.

conc_sym_rslts %>%
    mutate(comp_data = if_else(grepl(pattern = 'HF', x = sym_method), 'Temp-Heat Flux', 'Temp-Only')) %>%
    filter(comp_data =='Temp-Only' & model %in% gsub(pattern = 'conc_', replacement = '', x = example_models)) %>%
    rename(kappa = diff) ->
    S_kappa

names(COLOR_THEME)[4:6] <- example_models

ggplot(data = S_kappa, aes(kappa, optmized_value, color = model, shape = comp_data)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression('Optmial Value'),
         x = expression(kappa~'('~cm^2*s^-1~')')) +
    coord_cartesian(ylim = c(0, 0.4)) +
    guides(shape = FALSE) +
    scale_color_manual(values = example_models_COLORS) +
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
    scale_color_manual(values = example_models_COLORS) +
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

## Now make the symetry plot with the temperature heat flux data
conc_sym_rslts %>%
    mutate(comp_data = if_else(grepl(pattern = 'tempHF', x = sym_method), 'Temp-Heat Flux', 'Temp-Only')) %>%
    filter(comp_data == 'Temp-Heat Flux') %>%
    filter(model %in% example_models) %>%
    rename(kappa = diff) ->
    S_kappa_tempHF

ggplot(data = S_kappa_tempHF, aes(kappa, optmized_value, color = model, shape = comp_data)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression('Optmial Value'),
         x = expression(kappa~'('~cm^2*s^-1~')')) +
    coord_cartesian(ylim = c(0, 1)) +
    guides(shape = FALSE) +
   scale_color_manual(values = example_models_COLORS) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color = 'black')) ->
    conc_S_kappa_sym_plot_HF


ggplot(data = S_kappa_tempHF, aes(kappa, S, color = model, shape = comp_data)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression(~italic(S)~' ('~degree~'C)'),
         x = expression(kappa~'('~cm^2*s^-1~')')) +
    guides(shape = FALSE) +
    scale_color_manual(values = example_models_COLORS) +
    coord_cartesian(ylim = c(0, 7)) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color = 'black')) ->
    conc_S_kappa_trade_off_HF


conc_S_kapp_sym_plot_final_HF <- plot_grid(conc_S_kappa_sym_plot_HF,
                                           conc_S_kappa_trade_off_HF, labels = c('A', 'B'),
                                           label_size = 10)

ggsave(conc_S_kapp_sym_plot_final_HF,
       filename = file.path(FIGS_DIR, 'conc_S_kappa_sym_plot_HF.pdf'),
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
       aes(x)) + stat_function(fun=logmesaplt(-1,1,0.15), size=0.50) +
    ylab('L(x)') +
    SCRIPT_THEME ->
    mesaplt

ggsave(mesaplt,
       filename = file.path(FIGS_DIR, 'log_mesa_plot.pdf'),
       device = 'pdf', width = FIG_WIDTH / 2 , height = (FIG_WIDTH / FIG_RATIO)  ,
       dpi = FIG_DPI * 2, units = 'in')

## 9. Supplemental Figures -------------------------------------------------------------------------------------
## Temp-Only calibration results
# Make a plot comparing the Hector and the ESM comparison data.
emulated_conc %>%
    filter(method == 'temp' & variable == 'tas') %>%
    rename(Hector = value) %>%
    mutate(model = gsub(pattern = 'conc_', replacement = '', x = model)) %>%
    left_join(cmip_individual %>%
                  select(model, experiment, year, variable, ensemble, CMIP = value),
              by = c("model", "experiment", "year", "variable")) %>%
    na.omit %>%
    split(.$model, drop = TRUE) %>%
    lapply(function(input){
        ggplot(data = input) +
            geom_line(aes(year, CMIP, group = interaction(model, ensemble, experiment),
                          color = 'CMIP'),
                      size = 1,
                      alpha = 0.5) +
            geom_line(aes(year, Hector, color = 'Hector',
                          group = interaction(model, ensemble, experiment)),
                      size = 0.5) +
            scale_color_manual(values = COLOR_THEME) +
            SCRIPT_THEME +
            labs(x = 'Year',
                 y = unique(input$units)) +
            geom_text(data = tibble(x=1850,
                                    y = max(input$Hector, input$CMIP),
                                    label = unique(input$model)),
                      size = 3,
                      aes(x, y, label = label, hjust = 'left')) +
            coord_cartesian(xlim = c(1850, 2100))

    }) ->
    plot_list

# Arrange the plots in a list in a grid based on the MSE value.
arranged_plot_list <- list()
for(i in 1:length(names(plot_list))){

    name <- names(plot_list)[i]
    arranged_plot_list[[name]] <- plot_list[[name]] + theme(legend.position = 'none')


        arranged_plot_list[[name]] <- arranged_plot_list[[name]] + labs(y = expression(degree~'C'))
}

# Extract the legend from one of the plots, so that there is a universal plot.
arranged_plots <- plot_grid(plotlist = arranged_plot_list[1:16], nrow = 4, ncol = 4)
legend         <- get_legend(plot_list$`ACCESS1-3` + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_conc_omparison_tempOnly_temp1.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')

# Extract the legend from one of the plots, so that there is a universal plot.
arranged_plots <- plot_grid(plotlist = arranged_plot_list[17:32], nrow = 4, ncol = 4)
legend <- get_legend(plot_list$`CMCC-CESM` + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_conc_comparison_tempOnly_temp2.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')


## The Temp - HF Results
# Make a plot comparing the Hector and the ESM comparison data.
emulated_conc %>%
    rename(Hector = value) %>%
    mutate(model = gsub(pattern = 'conc_', replacement = '', x = model)) %>%
    filter(method == 'temp heatflux') %>%
    left_join(cmip_individual %>%
                  select(model, experiment, year, variable, ensemble, CMIP = value),
              by = c("model", "experiment", "year", "variable")) %>%
    na.omit %>%
    filter(variable == 'tas') %>%
    split(.$model, drop = TRUE) %>%
    lapply(function(input){

        ggplot(data = input) +
            geom_line(aes(year, CMIP, group = interaction(model, ensemble, experiment), color = 'CMIP'),
                      size = 0.75,
                      alpha = 0.5) +
            geom_line(aes(year, Hector, color = 'Hector',
                          group = interaction(model, ensemble, experiment)),  size = 0.75) +
            scale_color_manual(values = COLOR_THEME) +
            SCRIPT_THEME +
            labs(x = 'Year',
                 y = expression(degree~'C')) +
            geom_text(data = tibble(x=1850,
                                    y = max(input$Hector, input$CMIP),
                                    label = unique(input$model)), size = 3,
                      aes(x, y, label = label, hjust = 'left')) +
            coord_cartesian(xlim = c(1850, 2100))

    }) ->
    plot_list_tempHF_temp

# Arrange the plots in a list in a grid based on the MSE value.
arranged_plot_list <- list()
for(i in 1:length(names(plot_list_tempHF_temp))){

    name <- names(plot_list)[i]
    arranged_plot_list[[name]] <- plot_list_tempHF_temp[[name]] + theme(legend.position = 'none')
}


# Extract the legend from one of the plots, so that there is a universal plot.
arranged_plots <- plot_grid(plotlist = arranged_plot_list[1:17], nrow = 4, ncol = 4)
legend         <- get_legend(plot_list_tempHF_temp$`CMCC-CESM` + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_conc_omparison_tempHF_temp1.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')

# Extract the legend from one of the plots, so that there is a universal plot.
arranged_plots <- plot_grid(plotlist = arranged_plot_list[17:32], nrow = 4, ncol = 4)
legend <- get_legend(plot_list$`CMCC-CESM` + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_conc_comparison_tempHF_temp2.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')


# The Heat Flux Plots are going to have to be boken up into the the following categories,
# the runs that specfic heat flux comparison data and the runs that used the cmip heat flux range.
# First we will plot all of the runs that used the specific heat flux data.
emulated_conc %>%
    rename(Hector = value) %>%
    mutate(model = gsub(pattern = 'conc_', replacement = '', x = model)) %>%
    filter(method == 'temp heatflux') %>%
    filter(comp_data == 'temp-heatflux' & variable == 'heatflux') %>%
    filter(experiment != 'historical') %>%
    left_join(cmip_individual %>%
                  select(model, experiment, year, variable, ensemble, CMIP = value),
              by = c("model", "experiment", "year", "variable")) %>%
    na.omit %>%
    split(.$model, drop = TRUE) %>%
    lapply(function(input){

        ggplot(data = na.omit(input)) +
            geom_line(aes(year, CMIP, group = interaction(experiment, ensemble), color = 'CMIP'), alpha = 0.5, size = 0.75) +
            geom_line( aes(year, Hector, group = experiment, color = 'Hector')) +
            scale_color_manual(values = COLOR_THEME) +
            SCRIPT_THEME +
            labs(x = 'Year',
                 y = unique(input$units)) +
            geom_text(data = tibble(x=2005,
                                    y = 8.5,
                                    label = unique(input$model)), size = 3,
                      aes(x, y, label = label, hjust = 'left')) +
            coord_cartesian(xlim = c(2005, 2100),
                            ylim = c(-0.15, 8.75))

    }) ->
    specfic_HF


## The pltos that compare the heat flux range values.
emulated_conc %>%
    rename(Hector = value) %>%
    mutate(model = gsub(pattern = 'conc_', replacement = '', x = model)) %>%
    filter(method == 'temp heatflux') %>%
    filter(comp_data == 'temp-heatfluxRange') %>%
    filter(experiment != 'historical') %>%
    select(model, year, Hector, experiment, variable) %>%
    left_join(esm_comparison, by = c("year", "experiment", "variable")) %>%
    split(.$model, drop = TRUE)  %>%
    lapply(function(input){


        cmip_individual %>%
            filter(variable == 'tas' & experiment != 'historical') %>%
            filter(model == unique(input$model)) %>%
            pull(experiment) %>%
            unique() ->
            tas_exp_list

        input %>%
            filter(experiment %in% tas_exp_list &
                   experiment != 'historical' &
                       variable == 'heatflux') ->
            cmip_range

        ggplot() +
            geom_ribbon(data = cmip_range, aes(year, ymin = mina, ymax = maxb, group = experiment, fill = 'CMIP'),
                        alpha = 0.25) +
            geom_line(data = cmip_range, aes(year, Hector, group = experiment, color = 'Hector')) +
            scale_color_manual(values = COLOR_THEME) +
            scale_fill_manual(values = COLOR_THEME) +
            SCRIPT_THEME +
            labs(x = 'Year',
                 y = unique(input$units)) +
            geom_text(data = tibble(x=2005,
                                    y = 8.5,
                                    label = unique(input$model)), size = 3,
                      aes(x, y, label = label, hjust = 'left')) +
            coord_cartesian(xlim = c(2005, 2100),
                            ylim = c(-0.15, 8.75))

    }) ->
    range_HF

# Combine the heat flux runs into a list
heat_flux_comparison_list_unordered <- append(specfic_HF, range_HF)
heat_flux_comparison_list_ordered  <- heat_flux_comparison_list_unordered[sort(names(heat_flux_comparison_list_unordered))]


# Arrange the plots in a list in a grid based on the MSE value.
arranged_plot_list <- list()
for(i in 1:length(names(heat_flux_comparison_list_ordered))){

    name <- names(heat_flux_comparison_list_ordered)[i]
    arranged_plot_list[[name]] <- heat_flux_comparison_list_ordered[[name]] +
        theme(legend.position = 'none') +
        labs(y = expression(Wm^-2))
}


# Make a quick plot with the characteristics for the legend.
tibble(Hector = 1:10,
       CMIP = 1:10,
       source = CMIP,
       lower = 1:10,
       upper = 2:11,
       time = 1:10) %>%
    ggplot() +
    geom_line(aes(time, Hector, color = 'Hector')) +
    geom_line(aes(time, CMIP, color = 'CMIP')) +
    geom_ribbon(aes(time, ymin = lower, ymax = upper, fill = 'CMIP Range'), alpha = 0.25) +
    scale_color_manual(values = COLOR_THEME) +
    scale_fill_manual(values = c('CMIP Range' = 'grey')) +
    SCRIPT_THEME ->
    legend_plot

legend <- get_legend(legend_plot)

# Extract the legend from one of the plots, so that there is a universal plot.
arranged_plots <- plot_grid(plotlist = arranged_plot_list[1:16], nrow = 4, ncol = 4)
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_conc_summpelmental_comparison_HFtempHF1.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')

# Extract the legend from one of the plots, so that there is a universal plot.
arranged_plots <- plot_grid(plotlist = arranged_plot_list[17:32], nrow = 4, ncol = 4)
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_conc_summpelmental_comparison_HFtempHF2.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')





### B. Emission Driven Results ##################################################################################
## Select Example Models
# Save a vector of the example emission driven models to plot in manuscript figures, the other models will be
# visualized in the supplemental materials.
example_emiss_models <- c('GFDL-ESM2G', 'CESM1-BGC', 'CanESM2')
EMISS_COLORS         <- COLOR_THEME[5:7]
names(EMISS_COLORS)  <- example_emiss_models

## 1. Results Summary ----------------------------------------------------------------------------------------

# Summarise the paramter fits for the two different emission driven calibration protocols.
emiss_fits %>%
    gather(param, value, S, diff, alpha, volscl, C0, beta, q10_rh) %>%
    filter(penalty == 'without') %>%
    group_by(param, penalty) %>%
    summarise(min = min(value),
              mean = mean(value),
              max = max(value))

## 2. Change in Dot Plots -------------------------------------------------------------------------------------
# The color codes for the

names(CHANGE_PLOT_COLORS) <- c('without', 'with')

# Add a factor to the fits data frame so that the models will be plotted in alphabetical order.
emiss_fits$model <- factor(x = emiss_fits$model, levels = sort(unique(emiss_fits$model), decreasing = TRUE), ordered = TRUE)
emiss_fits$penalty <- factor(x = emiss_fits$penalty, levels = c('without', 'with'), ordered = TRUE)
emiss_fits <- arrange(emiss_fits, model, penalty)


ggplot(data = emiss_fits) +
    geom_point(aes(y = model, x = S, color = penalty),
               size = 2.5) +
    geom_path(aes(y = model, x = S),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'last')) +
    labs(y = NULL)  +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    labs(x = expression(~italic(S)~' ('~degree~'C)')) ->
    change_fitted_S

ggplot(data = emiss_fits) +
    geom_point(aes(y = model, x = diff, color = penalty),
               size = 2.5) +
    geom_path(aes(y = model, x = diff),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'last')) +
    labs(y = NULL) +
    coord_cartesian(xlim = c(-1, 20)) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    labs(x = expression(kappa~' ('~cm^2/'s)')) ->
    change_fitted_diff

ggplot(data = emiss_fits) +
    geom_point(aes(y = model, x = alpha, color = penalty),
               size = 2.5) +
    geom_path(aes(y = model, x = alpha),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'last')) +
    labs(y = NULL) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS)+
    labs(x = expression(alpha[a])) ->
    change_fitted_alpha

ggplot(data = emiss_fits) +
    geom_point(aes(y = model, x = volscl, color = penalty),
               size = 2.5) +
    geom_path(aes(y = model, x = volscl),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'last')) +
    labs(y = NULL) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    labs(x = expression(alpha[v])) ->
    change_fitted_volscl

ggplot(data = emiss_fits) +
    geom_point(aes(y = model, x = C0, color = penalty),
               size = 2.5) +
    geom_path(aes(y = model, x = C0),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'last')) +
    labs(y = NULL) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    labs(x = expression(C[0]~'(ppmv '~CO[2]~')')) ->
    change_fitted_C0


ggplot(data = emiss_fits) +
    geom_point(aes(y = model, x = beta, color = penalty),
               size = 2.5) +
    geom_path(aes(y = model, x = beta),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'last')) +
    labs(y = NULL) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    labs(x = expression(beta)) ->
    change_fitted_beta


ggplot(data = emiss_fits) +
    geom_point(aes(y = model, x = q10_rh, color = penalty),
               size = 2.5) +
    geom_path(aes(y = model, x = q10_rh),
              arrow = arrow(length=unit(0.20,"cm"), ends = 'last')) +
    labs(y = NULL) +
    SCRIPT_THEME +
    scale_color_manual(values = CHANGE_PLOT_COLORS) +
    labs(x = expression(Q[10])) ->
    change_fitted_q10_rh

change_param_plot <- plot_grid(change_fitted_S + theme(legend.position = 'none'),
                               change_fitted_diff + theme(legend.position = 'none'),
                               change_fitted_alpha + theme(legend.position = 'none'),
                               change_fitted_volscl + theme(legend.position = 'none'),
                               change_fitted_C0 + theme(legend.position = 'none'),
                               change_fitted_beta + theme(legend.position = 'none'),
                               change_fitted_q10_rh + theme(legend.position = 'none'),
                               labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'),
                               label_size = 10, ncol = 2)

legend <- get_legend(change_fitted_S + theme(legend.box.margin = margin(0, 0, 0, 12)))

final_change_param_plot <- plot_grid(change_param_plot, legend, ncol = 1, rel_heights = c(4, .25))
# Save the figure
ggsave(final_change_param_plot,
       filename = file.path(FIGS_DIR, 'emiss_change_fitted_param.pdf'),
       device = 'pdf', width = FIG_WIDTH , height = FIG_WIDTH ,
       dpi = FIG_DPI * 3, units = 'in')





## 3. Fitted Paramter Table ----------------------------------------------------------------------------------

emiss_fits %>%
    gather(param, value, S, diff, alpha, volscl, C0, beta, q10_rh) %>%
    select(model, penalty, param, value) %>%
    mutate(param = if_else(param == 'alpha', '1', param),
           param = if_else(param == 'volscl', '2', param),
           param = if_else(param == 'S', '3', param),
           param = if_else(param == 'diff', '4', param),
           param = if_else(param == 'C0', '5', param),
           param = if_else(param == 'beta', '6', param),
           param = if_else(param == 'q10_rh', '7', param),
           param = as.integer(param)) %>%
    mutate(method = if_else(penalty == 'without', 1, 2)) %>%
    mutate(param = paste0(method, param)) %>%
    select(model, param, value) %>%
    spread(param, value) %>%
    na.omit %>%
    kable(format = 'latex', digits = 2, align = 'c',
          col.names = c('Model', 'alpha', 'volscl', 'S', 'diff', 'C0', 'beta', 'Q10',
                        'alpha', 'volscl', 'S', 'diff', 'C0', 'beta', 'Q10')) %>%
    kable_styling(c("striped"), full_width = F) %>%
    add_header_above(c(" " = 1, "Without Penalty" = 7, "With Q10 Penalty" = 7))  %>%
    save_kable(file = file.path(FIGS_DIR, "emiss_fitted_params.pdf"))


## 4. How did the Calibration Method Change the Emulation Results? -------------------------------------------
# Caluclate the difference between the output variables when the calibration protocol changed.
emulated_emiss %>%
    filter(model %in% example_emiss_models) %>%
    select(model, penalty, value, variable, year) %>%
    spread(penalty, value) %>%
    mutate(SE = (with - without)^2) %>%
    group_by(model, variable) %>%
    summarise(RMSE = mean(SE)) %>%
    ungroup %>%
    na.omit() %>%
    group_by(variable) %>%
    summarise(RMSE_SD = sd(RMSE^1/2),
              RMSE = mean(RMSE^1/2))


## 5. Symmetry Between Q10 and Beta ---------------------------------------------------------------------------

# Without the penalty beta and q10 symmetry plots.
emiss_sym %>%
    filter(sym_method == 'emiss_carbon') %>%
    filter(model %in% example_emiss_models) ->
    beta_q10

ggplot(data = beta_q10, aes(beta, optmized_value, color = model)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression('Optmial Value'),
         x = expression(beta)) +
    coord_cartesian(ylim = c(0, 0.7)) +
    guides(shape = FALSE) +
    scale_color_manual(values = EMISS_COLORS) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color = 'black')) ->
    emiss_beta_sym_plot

ggplot(data = beta_q10, aes(beta, q10_rh, color = model)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression(Q[10]),
         x = expression(beta)) +
    guides(shape = FALSE) +
    scale_color_manual(values = EMISS_COLORS) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color = 'black')) ->
    emiss_beta_tradeoff

emiss_beta_q10_sym_plot_final <- plot_grid(emiss_beta_sym_plot,
                                           emiss_beta_tradeoff, labels = c('A', 'B'),
                                        label_size = 10)

ggsave(emiss_beta_q10_sym_plot_final,
       filename = file.path(FIGS_DIR, 'emiss_beta_q10_sym_plot.pdf'),
       device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO * 2,
       dpi = FIG_DPI * 2, units = 'in')


# With the penalty beta and q10 symmetry plots.
emiss_sym %>%
    filter(sym_method != 'emiss_carbon') %>%
    filter(model %in% example_emiss_models) ->
    beta_q10_penalty

ggplot(data = beta_q10_penalty, aes(beta, optmized_value, color = model)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression('Optmial Value'),
         x = expression(beta)) +
    coord_cartesian(ylim = c(0, 0.7)) +
    guides(shape = FALSE) +
    scale_color_manual(values = EMISS_COLORS) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color = 'black')) ->
    emiss_beta_sym_plot

ggplot(data = beta_q10_penalty, aes(beta, q10_rh, color = model)) +
    geom_line() +
    geom_point(size = 2.5) +
    SCRIPT_THEME +
    labs(y = expression(Q[10]),
         x = expression(beta)) +
    guides(shape = FALSE) +
    coord_cartesian(ylim = c(1.6, 2)) +
    scale_color_manual(values = EMISS_COLORS) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color = 'black')) ->
    emiss_beta_tradeoff

emiss_beta_q10_sym_plot_final <- plot_grid(emiss_beta_sym_plot,
                                           emiss_beta_tradeoff, labels = c('A', 'B'),
                                           label_size = 10)

ggsave(emiss_beta_q10_sym_plot_final,
       filename = file.path(FIGS_DIR, 'emiss_beta_q10_penalty_sym_plot.pdf'),
       device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH / FIG_RATIO * 2,
       dpi = FIG_DPI * 2, units = 'in')






## 6. Plot Hector vs Comp Data ---------------------------------------------------------

# This is just for the first calibration protocol.
emulated_emiss %>%
    filter(model %in% example_emiss_models) %>%
    filter(method == 'emiss_temp-heatflux-co2') %>%
    filter(!(experiment == 'esmHistorical' & variable == 'heatflux')) %>%
    rename(Hector = value) %>%
    left_join(cmip_individual %>%
                  select(model, experiment, year, variable, CMIP = value, ensemble) %>%
                  filter(grepl('esm', experiment)),
              by = c("model", "year", "variable", "experiment")) %>%
    na.omit() %>%
    split(interaction(.$model, .$variable), drop = TRUE) %>%
    lapply(function(input){

        ggplot(data = input) +
            geom_line(aes(year, CMIP, group = interaction(model, ensemble), color = 'CMIP'),
                      size = 0.75,
                      alpha = 0.5) +
            geom_line(aes(year, Hector, color = 'Hector'),  size = 1, linetype = 2) +
            scale_color_manual(values = COLOR_THEME) +
            SCRIPT_THEME +
            labs(x = 'Year',
                 y = unique(input$units)) +
            geom_text(data = tibble(x=1850,
                                    y = max(input$Hector, input$CMIP),
                                    label = unique(input$model)),
                      aes(x, y, label = label, hjust = 'left')) +
            coord_cartesian(xlim = c(1850, 2100))

    }) ->
    comparison_plots


plot_list <- list(comparison_plots$CanESM2.co2 +
                      theme(legend.position = 'none') + labs(y = expression('ppmv CO'[2])),
                  comparison_plots$`CESM1-BGC.co2` +
                      theme(legend.position = 'none') + labs(y = expression('ppmv CO'[2])),
                  comparison_plots$`GFDL-ESM2G.co2`+
                      theme(legend.position = 'none') + labs(y = expression('ppmv CO'[2])),
                  comparison_plots$CanESM2.heatflux +
                      theme(legend.position = 'none') + labs(y = expression(Wm^-2)),
                  comparison_plots$`CESM1-BGC.heatflux` +
                      theme(legend.position = 'none') + labs(y = expression(Wm^-2)),
                  comparison_plots$`GFDL-ESM2G.heatflux` +
                      theme(legend.position = 'none') + labs(y = expression(Wm^-2)),
                  comparison_plots$CanESM2.tas +
                      theme(legend.position = 'none') + labs(y = expression(degree~C)),
                  comparison_plots$`CESM1-BGC.tas`+
                      theme(legend.position = 'none') + labs(y = expression(degree~C)),
                  comparison_plots$`GFDL-ESM2G.tas`+
                      theme(legend.position = 'none') + labs(y = expression(degree~C)))


arranged_plots <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 3, labels = LETTERS[1:length(plot_list)])

# Extract the legend from one of the plots, so that there is a universal plot.
legend <- get_legend(comparison_plots$`GFDL-ESM2G.tas` + theme(legend.box.margin = margin(0, 0, 0, 12)))

# Add the legend to the row we made earlier and save.
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'emiss_hector_cmip_comparison.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')


## 6. Supplemental Figures ---------------------------------------------------------------------------------

# The figures the emission driven runs. Without the Q10 penalty.
emulated_emiss %>%
    filter(penalty == 'without') %>%
    mutate(keep = FALSE,
           keep = if_else(variable != 'heatflux', TRUE, keep),
           keep = if_else(variable == 'heatflux' & year > 2005, TRUE, keep)) %>%
    filter(keep) %>%
    select(-keep) %>%
    rename(Hector = value) %>%
    left_join(cmip_individual %>%
                  select(model, year, CMIP = value, variable, ensemble, experiment),
                         by = c("model", "year", "variable", "experiment")) %>%
    # Modify the ensemble name for the
    mutate(ensemble = if_else(model == 'CanESM2' & variable == 'heatflux' & year > 2055 &
                                  ensemble == 'r1i1p1', 'new', ensemble)) %>%
na.omit ->
     to_plot

 to_plot %>%
    split(., interaction(.$variable, .$model)) %>%
    lapply(function(input){

        ggplot(data = input) +
            geom_line(aes(year, CMIP, group = interaction(model, ensemble, experiment, drop = TRUE),
                          color = 'CMIP'),
                      size = 1,
                      alpha = 0.5, na.rm = FALSE) +
            geom_line(aes(year, Hector, color = 'Hector'),
                      size = 0.5) +
            scale_color_manual(values = COLOR_THEME) +
            SCRIPT_THEME +
            labs(x = 'Year',
                 y = unique(input$units)) +
            geom_text(data = tibble(x= min(input$year),
                                    y = max(input$Hector, input$CMIP, na.rm = TRUE),
                                    label = unique(input$model)), size = 3,
                      aes(x, y, label = label, hjust = 'left'))
    }) ->
    plot_list

 # Because MIRCO-ESM heat flux data stopsin 2056 we had to use the heat flux range to constrain the model,
 # use the range to constrain the model so the plot acutally needs to be of the CMIP range.
 emulated_emiss %>%
     filter(model == 'MIROC-ESM' & variable == 'heatflux' & experiment == 'esmrcp85') %>%
     select(model, experiment, year, variable, Hector = value, penalty) %>%
     left_join(esm_comparison, by = c("experiment", "year", "variable")) %>%
     split(interaction(.$model, .$penalty), drop = TRUE) %>%
     lapply(function(input){

         ggplot(data = input) +
             geom_ribbon(aes(year, ymin = mina, ymax = maxb, group = experiment, fill = 'CMIP'),
                         alpha = 0.25) +
             geom_line(aes(year, Hector, group = experiment, color = 'Hector')) +
             scale_color_manual(values = COLOR_THEME) +
             scale_fill_manual(values = COLOR_THEME) +
             SCRIPT_THEME +
             labs(x = 'Year',
                  y = unique(input$units)) +
             geom_text(data = tibble(x=2005,
                                     y = 8.5,
                                     label = unique(input$model)), size = 3,
                       aes(x, y, label = label, hjust = 'left')) +
             coord_cartesian(xlim = c(2005, 2100),
                             ylim = c(-0.15, 8.75))

     }) ->
     emiss_heatflux_range


co2_plots <- plot_list[grepl(pattern = 'co2.', names(plot_list))]
hf_plots  <- plot_list[grepl(pattern = 'heatflux.', names(plot_list))]
hf_plots$`heatflux.MIROC-ESM` <- emiss_heatflux_range$`MIROC-ESM.without`

tas_plots <- plot_list[grepl(pattern = 'tas.', names(plot_list))]




# Arrange the co2 plots in a list that can be plotted as a grid without
# the legends.
co2_plot_list <- list()
for(i in 1:length(names(co2_plots))){

    name <- names(co2_plots)[i]
    co2_plot_list[[name]] <- co2_plots[[name]] +
        theme(legend.position = 'none') +
        labs(y = expression('ppmv CO'[2]))
}

legend <- get_legend(co2_plots$co2.CanESM2)

arranged_plots <- plot_grid(plotlist = co2_plot_list, ncol = 3)
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_emiss_withoutPen_co2.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')


# Arrange the temp plots in a list that can be plotted as a grid without
# the legends.
temp_plot_list <- list()
for(i in 1:length(names(tas_plots))){

    name <- names(tas_plots)[i]
    temp_plot_list[[name]] <- tas_plots[[name]] +
        theme(legend.position = 'none') +
        labs(y = expression(degree~C))
}

legend <- get_legend(tas_plots$tas.CanESM2)

arranged_plots <- plot_grid(plotlist = temp_plot_list, ncol = 3)
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_emiss_withoutPen_tas.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')



# Arrange the temp plots in a list that can be plotted as a grid without
# the legends.
hf_plot_list <- list()
for(i in 1:length(names(hf_plots))){

    name <- names(hf_plots)[i]
    hf_plot_list[[name]] <- hf_plots[[name]] +
        theme(legend.position = 'none') +
        labs(y = expression(Wm^-2))
}

legend <- get_legend(hf_plots$heatflux.CanESM2)

arranged_plots <- plot_grid(plotlist = hf_plot_list, ncol = 3)
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_emiss_withoutPen_hf.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')


# The figures the emission driven runs. With the Q10 penalty.
emulated_emiss %>%
    filter(penalty == 'with') %>%
    mutate(keep = FALSE,
           keep = if_else(variable != 'heatflux', TRUE, keep),
           keep = if_else(variable == 'heatflux' & year > 2005, TRUE, keep)) %>%
    filter(keep) %>%
    select(-keep) %>%
    rename(Hector = value) %>%
    left_join(cmip_individual %>%
                  select(model, year, CMIP = value, variable, ensemble, experiment),
              by = c("model", "year", "variable", "experiment")) %>%
    # Modify the ensemble name for the
    mutate(ensemble = if_else(model == 'CanESM2' & variable == 'heatflux' & year > 2055 &
                                  ensemble == 'r1i1p1', 'new', ensemble)) %>%
    split(., interaction(.$variable, .$model)) %>%
    lapply(function(input){

        ggplot(data = input) +
            geom_line(aes(year, CMIP, group = interaction(model, ensemble, experiment, drop = TRUE),
                          color = 'CMIP'),
                      size = 1,
                      alpha = 0.5, na.rm = FALSE) +
            geom_line(aes(year, Hector, color = 'Hector'),
                      size = 0.5) +
            scale_color_manual(values = COLOR_THEME) +
            SCRIPT_THEME +
            labs(x = 'Year',
                 y = unique(input$units)) +
            geom_text(data = tibble(x= min(input$year),
                                    y = max(input$Hector, input$CMIP, na.rm = TRUE),
                                    label = unique(input$model)), size = 3,
                      aes(x, y, label = label, hjust = 'left'))
    }) ->
    plot_list

co2_plots <- plot_list[grepl(pattern = 'co2.', names(plot_list))]
hf_plots  <- plot_list[grepl(pattern = 'heatflux.', names(plot_list))]
hf_plots$`heatflux.MIROC-ESM` <- emiss_heatflux_range$`MIROC-ESM.with`
tas_plots <- plot_list[grepl(pattern = 'tas.', names(plot_list))]


# Arrange the co2 plots in a list that can be plotted as a grid without
# the legends.
co2_plot_list <- list()
for(i in 1:length(names(co2_plots))){

    name <- names(co2_plots)[i]
    co2_plot_list[[name]] <- co2_plots[[name]] +
        theme(legend.position = 'none') +
        labs(y = expression('ppmv CO'[2]))
}

legend <- get_legend(co2_plots$co2.CanESM2)

arranged_plots <- plot_grid(plotlist = co2_plot_list, ncol = 3)
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_emiss_withPen_co2.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')


# Arrange the temp plots in a list that can be plotted as a grid without
# the legends.
temp_plot_list <- list()
for(i in 1:length(names(tas_plots))){

    name <- names(tas_plots)[i]
    temp_plot_list[[name]] <- tas_plots[[name]] +
        theme(legend.position = 'none') +
        labs(y = expression(degree~C))
}

legend <- get_legend(tas_plots$tas.CanESM2)

arranged_plots <- plot_grid(plotlist = temp_plot_list, ncol = 3)
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_emiss_withPen_tas.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')



# Arrange the temp plots in a list that can be plotted as a grid without
# the legends.
hf_plot_list <- list()
for(i in 1:length(names(hf_plots))){

    name <- names(hf_plots)[i]
    hf_plot_list[[name]] <- hf_plots[[name]] +
        theme(legend.position = 'none') +
        labs(y = expression(Wm^-2))
}

legend <- get_legend(hf_plots$heatflux.CanESM2)

arranged_plots <- plot_grid(plotlist = hf_plot_list, ncol = 3)
plot_grid(arranged_plots, legend, rel_widths = c(4, .5)) +
    ggsave(filename = file.path(FIGS_DIR, 'summpelmental_emiss_withPen_hf.pdf'),
           device = 'pdf', width = FIG_WIDTH, height = FIG_WIDTH,
           dpi = FIG_DPI, units = 'in')


