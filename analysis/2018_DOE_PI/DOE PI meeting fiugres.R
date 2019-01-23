# Produce the figures for the 2018 DOE PI Meeting poster
# Bayesian Calibration of a Simple Cliamte Model: Hector
#
# This script uses results from interactive mcmc runs, while do not have the code that produces the
# mcmc results.
#
# See section 0 for notes on script set up.

# 0. Set Up -----

# The base directory should be set to analysis/2018_DOE_PI
BASE_DIR <- './analysis/2018_DOE_PI'

# load libs
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(purrr)

# Create a list to store the plots in
output_list <- list()

# The names of the parameters
param_names <- c('ECS', 'AERO', 'DIFF', 'BETA', 'Q10', 'CO2_0', 'Sigma tgav', 'sigma co2')


# 1. Import Data ----

# Import the results
mean_rlst <- readRDS(file.path(BASE_DIR, "mean-cal-fixsigs.rds"))
con_rlst  <- readRDS(file.path(BASE_DIR, "cons-cal-sig5.rds"))

hector_data_mean_rcp85 <- read.csv(file.path(BASE_DIR, 'mean_param_hector.csv'), stringsAsFactors = FALSE)
hector_data_cons_rcp85 <- read.csv(file.path(BASE_DIR, 'cons_param_hector.csv'), stringsAsFactors = FALSE)

hector_data_mean_rcp26 <- read.csv(file.path(BASE_DIR, 'mean_param_hector_rcp26.csv'), stringsAsFactors = FALSE)
hector_data_cons_rcp26 <- read.csv(file.path(BASE_DIR, 'con_param_hector_rcp26.csv'), stringsAsFactors = FALSE)

# 2. Check the MCMC sampeler ----
# Randomly select N parameters, this is going to be used to make the mean
# paramter plots that will be used to check that the mcmc sampler is not getting
# stuck somewhere.
N <- 500
mean_param_ids <- floor(seq(1, nrow(mean_rlst$samples), length.out = N))
mean_params    <- mean_rlst$samples[mean_param_ids, ]

# Tracer plots to make sure the sampler is working
for(i in 1:ncol(mean_params)){
    plot(mean_params[, i], ylab = param_names[i], xlab = 'Sampel Index', main = 'Tracer Plot')
}

# Alright so these plots show that the sampler is not getting suck anywhere which is awesome.


# 3. Find the MAP -----
# Pull out the MAP (maxiumum a posterior value)
mean_mpa_param <- mean_rlst$samples[which.max(mean_rlst$proplp), ]
mean_params_id <- cbind(mean_params, as.character(1:nrow(mean_params)))

# Combine into a single data frame
mean_map_id        <- c(mean_mpa_param, 'MAP')
mean_hector_params <- rbind(mean_params_id, mean_map_id)
# write.csv(mean_hector_params, 'mean_hector_params.csv', row.names = FALSE)


# 3. Plots -----
## Marginal probability density -----
## A mulipaneled plot comparing the marginal probability density for each of the
## hector parameters comparing the posterior parameter space of the mean calibrated
## vs the full range calibrted mcmc restults.

# First format the results
# Save the sample results as a tibble, add paramter names and mcmc method information.
mean_data <- tibble::as_tibble(mean_rlst$samples)
cons_data <- tibble::as_tibble(con_rlst$samples)

colnames(mean_data) <- c('ECS', 'AERO', 'DIFF', 'BETA', 'Q10', 'CO2_0', 'Sigma tgav', 'sigma co2')
colnames(cons_data) <- c('ECS', 'AERO', 'DIFF', 'BETA', 'Q10', 'CO2_0')

mean_data$`Calibration Method` <- 'Mean'
cons_data$`Calibration Method` <- 'Full range'

# Combine the results into a single long data frame, this will make plotting easier.
cons_data %>%
    bind_rows(mean_data) %>%
    gather(param, value, -`Calibration Method`) %>%
    # Discard the sigma values since the sigma paramter is not applicable to the mean
    # mcmc calibration method.
    filter(!grepl('sigma', param)) ->
    to_plot

# Add factor information to to calibration method to control the plotting order.
to_plot$`Calibration Method` <- factor(to_plot$`Calibration Method` , levels = c('Mean', 'Full range'), ordered = TRUE)

# For each of the parameters plot the probability density of the mean calibrated vs the full
# range calibrated paramter space.
#
# Args
#   input: a data frame of value / param / `Calibration Method`
# Returns: a density plot that compare the posterior param values
plot_posterior_density <- function(input){

    ggplot(data = input) +
        geom_density(aes(value, fill = `Calibration Method`, color = `Calibration Method`),
                     alpha = .5) +
        theme_bw() +
        labs(x = unique(input$param), y = NULL, title = NULL) +
        theme_minimal(base_size = 44) +
        theme(strip.placement='outside') +
        # The DOE PI meeting color scheme
        scale_fill_manual(values = c(`Full range` = '#268bd2', Mean = '#dc322f')) +
        scale_color_manual(values = c(`Full range` = '#268bd2', Mean = '#dc322f')) +
        NULL
}

# Split up the paramer data to plot into a list and plot using map.
posterior_param_list <- split(to_plot, to_plot$param)
param_plots          <- map(posterior_param_list, plot_posterior_density)

# Modify the denisty plots to have the corret labels
# and y axis.
param_plots$AERO +
    guides(fill = FALSE, color = FALSE) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(x = expression(alpha)) +
    theme(axis.title.x =element_text(face="italic"))  ->
    paero

param_plots$BETA +
    guides(fill = FALSE, color = FALSE) +
    coord_cartesian(ylim = c(0, 20)) +
    labs(x = expression(beta))->
    pbeta

param_plots$CO2_0 +
    guides(fill = FALSE, color = FALSE) +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(x = expression(C[0]) ) +
    scale_x_continuous(breaks = c(275, 285, 295)) ->
    pC0

param_plots$DIFF +
    guides(fill = FALSE, color = FALSE) +
    coord_cartesian(ylim = c(0, 3)) +
    labs(x = expression(kappa))->
    pdiff

param_plots$ECS +
    guides(fill = FALSE, color = FALSE) +
    coord_cartesian(ylim = c(0, 4)) +
    labs(x = 'S')->
    pecs

param_plots$Q10 +
    guides(fill = FALSE, color = FALSE) +
    coord_cartesian(ylim = c(0, 4)) +
    labs(x = expression(Q[10]))->
    pq10

# Format all of the plots as a multi pannled plot
# First save the legend from one plot
legend <- get_legend(out$Q10)

# Arrange the plots as a grid
param_plot_grid <- grid.arrange(grobs = list(pecs, paero, pdiff, pbeta, pq10, pC0),
                      ncol = 3)
# Add the legend
posterior_param_density_plot <- plot_grid(param_plot_grid, legend, rel_widths = c(3.5, 1.5))

# Save in the output list
output_list[['posterior_param_density_plot']] <- posterior_param_density_plot


# plots
# 2048

## Tgav rcp 85 spaghetti  -----

# Format the data
hector_data_mean_rcp85 %>%
    mutate(`Calibration Method` = 'Mean') %>%
    filter(variable == 'Tgav') ->
    tgav_mean_rcp85

hector_data_cons_rcp85 %>%
    mutate(`Calibration Method` = 'Full Range') %>%
    filter(variable == 'Tgav') ->
    tgav_cons_rcp85

# Combine into a single data frame
tgav_mean_rcp85 %>%
    bind_rows(tgav_cons_rcp85) %>%
    filter(year > 1950) %>%
    mutate(rcp = 'rcp 8.5') ->
    tgav_time_rcp85

MAP_runs <- filter(tgav_time_rcp85, run_id == 'MAP')

ggplot(data = tgav_time_rcp85) +
  geom_line(aes(year, value, group = interaction(run_id, `Calibration Method`), color = `Calibration Method`),
            alpha = 0.3) +
  geom_line(data = MAP_runs, aes(year, value, color = `Calibration Method`), size = 1.5) +
  theme_minimal(base_size = 44) +
  scale_color_manual(values = c(`Full Range` = '#268bd2', Mean = '#dc322f')) +
  labs(x = 'Year',
       y = expression(degree~'C'),
       title = 'RCP 8.5') ->
  output_list[['rcp85_spaghetti']]


## Tgav rcp 26 spaghetti  -----

# Format the data
hector_data_mean_rcp26 %>%
    mutate(`Calibration Method` = 'Mean') %>%
    filter(variable == 'Tgav') ->
    tgav_mean_rcp26

hector_data_cons_rcp26 %>%
    mutate(`Calibration Method` = 'Full Range') %>%
    filter(variable == 'Tgav') ->
    tgav_cons_rcp26

# Combine into a single data frame
tgav_mean_rcp26 %>%
    bind_rows(tgav_cons_rcp26) %>%
    filter(year > 1950) %>%
    mutate(rcp = 'rcp 8.5') ->
    tgav_time_rcp26

MAP_runs <- filter(tgav_time_rcp26, run_id == 'MAP')

ggplot(data = tgav_time_rcp26) +
    geom_line(aes(year, value, group = interaction(run_id, `Calibration Method`), color = `Calibration Method`),
              alpha = 0.3) +
    geom_line(data = MAP_runs, aes(year, value, color = `Calibration Method`), size = 1.5) +
    theme_minimal(base_size = 44) +
    scale_color_manual(values = c(`Full Range` = '#268bd2', Mean = '#dc322f')) +
    labs(x = 'Year',
         y = expression(degree~'C'),
         title = 'RCP 2.6') ->
    output_list[['rcp26_spaghetti']]

# End ---
