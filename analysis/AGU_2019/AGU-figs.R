## AGU-figs.R
## This script plots the results for the 2019 AGU Fall meeting.
# 0. Set Up -----------------------------------------------------------------------------
# Load the required packages.
library(dplyr)
library(tidyr)
library(ggplot2)
library(nationalparkcolors)
library(hectorcal)

# Define where to save the data.
OUTDIR <- here::here('analysis', 'AGU_2019')

# Define script figure aestheics.
THEME  <- theme_bw(base_size = 32)
COLORS <- park_palettes$Acadia

# 1. Kappa vs S Symmetry -----------------------------------------------------------------------------
# Load the summary results for when Hector is optimized at fixed kappa values. These results
# will be used to illsutrate how unless the cliamte system is properly constrained that it is
# not possible for the nonlinear optmization routine to identify a unique solution.
# Because the temp only and temp - heat flux optim values are going to be very different
# from one another the results will neeed to be centered around the mean so they can be
# reported on the same axis.
here::here('analysis', 'paper_1', 'output', '2A.symmetry_poster_results', 'conc_clim_temp') %>%
    list.files(pattern = 'summary_fits.csv', full.names = TRUE) %>%
    read.csv(stringsAsFactors = FALSE) %>%
    mutate(method = 'insufficently constrained') ->
    sym_restults_1

sym_restults_1 %>%
    group_by(model, method) %>%
    summarise(mean = mean(min)) %>%
    ungroup ->
    mean_df

sym_restults_1 %>%
    left_join(mean_df) %>%
    mutate(min = min - mean) ->
    sym_restults_1

here::here('analysis', 'paper_1', 'output', '2A.symmetry_poster_results', 'conc_clim_tempHF') %>%
    list.files(pattern = 'summary_fits.csv', full.names = TRUE) %>%
    read.csv(stringsAsFactors = FALSE) %>%
    mutate(method = 'sufficently constrained') %>%
    select(-value)  ->
    sym_restults_2

sym_restults_2  %>%
    group_by(model, method) %>%
    summarise(mean = mean(min)) %>%
    ungroup ->
    mean_df

sym_restults_2 %>%
    left_join(mean_df) %>%
    mutate(min = min - mean) %>%
    bind_rows(sym_restults_1) %>%
    filter(model == "GFDL-CM3") %>%
    ggplot(aes(kappa, min, color = method)) +
    geom_line(size = 3) +
    coord_cartesian(ylim = c(-0.2, 0.5)) +
    labs(x =  expression(kappa~cm^2~s^-1),
         y = 'Goodness of Fit') +
    scale_color_manual(values = c('insufficently constrained' = COLORS[[2]],
                                  'sufficently constrained' = COLORS[[3]])) +
    theme_bw(base_size = 26) +
    theme(legend.justification = c(1,1),
          legend.position = c(1,1),
          legend.background  = element_rect(color = 1),
          legend.title = element_blank()) +
    ggsave(filename = file.path(OUTDIR, 'goodness_of_fit.pdf'),
           device = 'pdf', width = 6, height = 6, units = 'in', dpi = 200)

# 2. Log Mesa Plot -----------------------------------------------------------------------------
# For the section of the paper that introduces the idea of the mesa function we would like to include a figure
# to help facilitate our explanation of what is going on.

## Function for generating a log-mesa function for use in plotting functions
logmesaplt <- function(a,b,sig) {
    function(x) { mesa(x, a, b, sig) }
}

ggplot(data=data.frame(x=c(-1.5, 1.5)),
       aes(x)) + stat_function(fun=logmesaplt(-1,1,0.15), size=3) +
    ylab('Likelihood') +
    labs(x = 'Value') +
    THEME ->
    mesaplt

ggsave(mesaplt,filename = file.path(OUTDIR, 'log_mesa_plot.pdf'), width = 11, height = 7)

# 3. Avaibale Data Table-----------------------------------------------------------------------------
# Make a table of the ESM data that is aviaible to use in calibration.

hectorcal::cmip_individual %>%
    filter(!grepl(pattern = 'esm', x = experiment)) %>%
    filter(variable %in% c('tas', 'heatflux')) %>%
    group_by(experiment, variable) %>%
    summarise('model count' = n_distinct(model)) %>%
    arrange(experiment, rev(variable)) %>%
    write.csv(file = file.path(OUTDIR, 'singel-esm-comp-data.csv'),  row.names = FALSE)

# 4. Bivariate S and Kappa -----------------------------------------------------------------------------
# In order to illustrate how adding the sufficent comaprison data can affect the Baeysian results
# plot mcmc results for the two different mcmc runs.

# Format the results that only used temperature as a data constraint.
size <- 10000
readRDS(here::here('analysis', 'AGU_2019', 'mcruns-conc-16.rds')) %>%
    as.data.frame() %>%
    mutate(method = 'insufficently constrained')->
    mcmc_Temp

mcmc_Temp <- mcmc_Temp[sample(1:nrow(mcmc_Temp), size), ]

# Format the rsults from the temperature - heat flux mcmc run.
size <- 10000
readRDS(here::here('analysis', 'AGU_2019', 'mcruns-conc-48.rds')) %>%
    as.data.frame() %>%
    mutate(method = 'sufficently constrained')->
    mcmc_TempHF
mcmc_TempHF <- mcmc_TempHF[sample(1:nrow(mcmc_Temp), size), ]

# Save the plot.
mcmc_Temp %>%
    bind_rows(mcmc_TempHF) %>%
    ggplot() +
    geom_point(aes(S, diff, color = method)) +
    scale_color_manual(values = c('insufficently constrained' = COLORS[[2]],
                                  'sufficently constrained' = COLORS[[3]])) +
    facet_wrap(~method) +
    THEME +
    theme(legend.position = 'none') +
    labs(x = expression(italic(S)),
         y = expression(kappa)) +
    ggsave(filename = file.path(OUTDIR, 'scatter_plots.pdf'), width = 12, height = 6)


# 5. Marginal PDf Mean vs Range -----------------------------------------------------------------------------
# In order to illsutrate how calibrating to the mean vs the range can impact the joint pdf returned by the
# mcmc plot the marignal pdfs for all of the climate paramters from two different mcmc chains.

readRDS(here::here('analysis', 'AGU_2019', 'mcruns-conc-48.rds')) %>%
    as.data.frame() %>%
    mutate(method = 'multi-model range')->
    baey_range

readRDS(here::here('analysis', 'AGU_2019', 'mcruns-conc-112.rds')) %>%
    as.data.frame() %>%
    mutate(method = 'multi-model mean')->
    baey_mean

labels <- c(expression(alpha[a]), expression(kappa), 'S', expression(alpha[v]))

baey_mean %>%
    bind_rows(baey_range) %>%
    gather(param, value, S, diff, alpha, volscl) %>%
    ggplot() +
    geom_density(aes(value, fill = method, color = method), alpha = 0.5) +
    facet_wrap(~param, nrow = 1, scales = 'free') +
    theme_bw(base_size = 28) +
    theme(legend.title = element_blank(),
          legend.position = 'bottom') +
    scale_color_manual(values = c('multi-model mean' = COLORS[6],
                                  'multi-model range' = COLORS[3])) +
    scale_fill_manual(values = c('multi-model mean' = COLORS[6],
                                 'multi-model range' = COLORS[3])) +
    labs(x = NULL,
         y = 'PDF') +
    ggsave(filename = file.path(OUTDIR, 'meav_v_range.pdf'), width = 17, height = 5.5)
