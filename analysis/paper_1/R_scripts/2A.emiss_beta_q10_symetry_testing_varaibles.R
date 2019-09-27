library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
#library(hector)
library(hectorcal)
library(cowplot)
devtools::load_all('~/Documents/hector')



BASE_DIR <- "/Users/dorh012/Documents/2019/hectorcal/analysis/paper_1"


beta_q10 <- read.csv(list.files(file.path(BASE_DIR, 'output', 'emiss_beta_q10'), '.csv', full.name = TRUE),
                     stringsAsFactors = FALSE)
beta_q10_penalty <- read.csv(list.files(file.path(BASE_DIR, 'output', 'emiss_beta_q10_penalty2'), '.csv',
                                        full.name = TRUE), stringsAsFactors = FALSE)

beta_q10$penalty         <- 'no penalty'
beta_q10_penalty$penalty <- 'penalty'
beta_q10_df <- bind_rows(beta_q10, beta_q10_penalty)

tibble::tibble( file = c(system.file('input/hector_rcp85.ini', package = 'hector'),
                         system.file('input/hector_rcp85.ini', package = 'hector')),
                experiment = c('esmHistorical', 'esmrcp85')) %>%
    select(ini_file = file, core_name = experiment, experiment) ->
    ini_files_tib

# Create the list of the Hector cores.
core_list <- mapply(newcore, inifile = ini_files_tib$ini_file, name = ini_files_tib$core_name )

split(beta_q10_df, interaction(beta_q10_df$model, beta_q10_df$penalty, beta_q10_df$beta, drop = TRUE)) %>%
lapply(function(input){

    input <- data.frame(input)
    param_names <- c(ECS(), BETA(), DIFFUSIVITY(), VOLCANIC_SCALE(), AERO_SCALE(), PREINDUSTRIAL_CO2(), Q10_RH())
    params      <- select(input, param_names)
    info        <- select(input, -param_names)

    lapply(core_list, parameterize_core, params = params)
    lapply(core_list, run, runtodate = 2100)
    lapply(core_list, fetchvars, dates = 1850:2100, vars = c(GLOBAL_TEMP(), HEAT_FLUX(), ATMOSPHERIC_CO2(),  NPP())) %>%
        bind_rows() %>%
        cbind(info) %>%
        cbind(params)

    }) %>%
    bind_rows() %>%
    filter(!(variable == "esmHistorical" & year <= 2005)) %>%
    filter(model  == 'CanESM2') ->
    df

ggplot(data = df) +
    geom_line(aes(year, value, group = interaction(variable, beta, q10_rh), color = penalty)) +
    facet_wrap('variable', scales = 'free')



split(df, df$variable) %>%
    lapply(function(input){

        ylab <- unique(input[['units']])

        ggplot(data = input) +
            geom_line(aes(year, value, group = interaction(beta, q10_rh))) +
            labs(y = ylab) +
            facet_wrap('penalty')
    }) ->
    list





