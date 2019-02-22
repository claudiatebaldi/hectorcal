#' CMIP5 ensemble comparison data
#'
#' This dataset contains summary statistics for the CMIP5 ensemble results by year,
#' variable, and experiment.
#'
#' @section Notes:
#'
#' No effort was made to account for the fact that some models have more
#' ensemble members in an experiment than others, nor for the fact that some
#' models have several (presumably closely related) variants in the ensemble
#' (e.g., CESM has CESM1-BGC, CESM1-CAM5, CESM1-FASTCHEM, and CESM1-WACCM
#' represented).  For extreame statistics like min, max, and probably the decile
#' points, this doesn't make a big difference, but central statistics like mean and
#' mdian, the difference in representation could give some models undue
#' influence on how the statistics turn out.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{year}{Year being described (1851--2300)}
#' \item{variable}{Variable being reported ("tas" or "co2").  Temperatures are
#' given as anomalies relative to the first decade reported by the model.}
#' \item{experiment}{Experiment being reported.  One of "historical", "rcp26", "rcp45",
#' "rcp60", "rcp85", "esmHistorical", "esmrcp85"}
#' \item{mina}{Minimum value reported in the ensemble}
#' \item{maxb}{Maximum value reported in the ensemble}
#' \item{a10}{First decile reported in the ensemble}
#' \item{b90}{Ninth decile reported in the ensemble}
#' \item{cmean}{Mean value reported in the ensemble}
#' \item{cmedian}{Median value reported in the ensemble}
#' }
#' @family comparison data
'esm_comparison'

#' CMIP5 individual model comparison data
#'
#' This dataset contains outputs for the individual models in the ensemble.
#'
#' @format Data frame with 8 columns
#' \describe{
#' \item{year}{Year being described (1851--2300)}
#' \item{model}{Name of the model}
#' \item{ensemble}{Identifier for the ensemble member}
#' \item{variable}{Variable being reported ("tas" or "co2").  Temperatures are
#' anomalies relative to the first decade reported by the model.}
#' \item{experiment}{Experiment being reported. One of "historical", "rcp26", "rcp45",
#' "rcp60", "rcp85", "esmHistorical", "esmrcp85"}
#' \item{value}{Value of the variable being reported}
#' \item{esmbaseline}{Mean temperature (in Kelvin) for this model's esmHistorical run during
#' its first decade.  Temperatures for esmHistorical and esmrcp85 experiments
#' are reported as anomalies relative to this temperature.  This column will be
#' \code{NA} for mdoels that did not run the esm experiments.}
#' \item{concbaseline}{Mean temperature (in Kelvin) for this model's historical run during
#' its first decade.  Temperatures for concentration-driven runs are reported as
#' anomalies relative to this temperatrure.}
#' }
#' @family comparison data
'cmip_individual'

#' Ensemble of Hector concentration driven runs
#'
#' This dataset contains outputs for the 1000 Hector runs driven with RCP 2.6, RCP 4.5, RCP 6.0
#' and RCP 8.5 concentration. Each Hector "case" or row corresponds to results from
#' a different combination of climate sensitivity, ocean heat diffusivity, and aerosol scaling factor
#' sampled from the parameter prior distribution.
#'
#' @format Data frame with 1000 rows and 1004 columns
#' \describe{
#' \item{Rows}{Represent results from 1000 Hector runs driven with concentrations
#' and varying climate sensitivity, ocean heat diffusivity, and aerosol scaling factor}
#' \item{Columns}{Hector temperature output values named with the following
#' pattern XYYYY_variable_scenario}
#' }
#' @family hector ensemble
'concen-RCP26_concen-RCP45_concen-RCP60_concen-RCP85'

#' Ensemble of Hector emission driven runs
#'
#' This dataset contains outputs for the 1000 Hector runs driven with RCP 2.6, RCP 4.5, RCP 6.0
#' and RCP 8.5 concentration. Each Hector "case" or row corresponds to results from
#' a different combination of climate sensitivity, ocean heat diffusivity, aerosol scaling factor,
#' beta, q10, and preindustrial CO2 sampled from the parameter prior distribution.
#'
#' @format Data frame with 1000 rows and 2008 columns
#' \describe{
#' \item{Rows}{Represent results from 1000 Hector runs driven with concentrations
#' and varying climate sensitivity, ocean heat diffusivity, and aerosol scaling factor}
#' \item{Columns}{Hector temperature and atmospheric output values named with the
#'  following pattern XYYYY_variable_scenario}
#' }
#' @family hector ensemble
'emiss-CC-RCP26_emiss-CC-RCP45_emiss-CC-RCP60_emiss-CC-RCP85'

#' Ensemble of Hector emission driven runs with constant Carbon cycle parameters
#'
#' This dataset contains outputs for the 1000 Hector runs driven with RCP 2.6, RCP 4.5, RCP 6.0
#' and RCP 8.5 concentration. Each Hector "case" or row corresponds to results from
#' a different combination of climate sensitivity, ocean heat diffusivity, and
#' aerosol scaling factor sampled from the parameter prior distribution.
#'
#' @format Data frame with 1000 rows and 1004 columns
#' \describe{
#' \item{Rows}{Represent results from 1000 Hector runs driven with concentrations
#' and varying climate sensitivity, ocean heat diffusivity, and aerosol scaling factor}
#' \item{Columns}{Hector temperature output values named with the
#'  following pattern XYYYY_variable_scenario}
#' }
#' @family hector ensemble
'emiss-consatntC-RCP26_emiss-constantC-RCP45_emiss-constantC-RCP60_emiss-constantC-RCP85'


