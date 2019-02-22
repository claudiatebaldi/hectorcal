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

#' Ensemble of Hector rcp26 concentration driven runs
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 2.6 concentration, this
#' data set only contains results from runs that also solved in for the Hector concentration driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-concen-rcp26'

#' Ensemble of Hector rcp45 concentration driven runs
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 4.5 concentration, this
#' data set only contains results from runs that also solved in for the Hector concentration driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-concen-rcp45'

#' Ensemble of Hector rcp60 concentration driven runs
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 6.0 concentration, this
#' data set only contains results from runs that also solved in for the Hector concentration driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-concen-rcp60'

#' Ensemble of Hector rcp85 concentration driven runs
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 8.5 concentration, this
#' data set only contains results from runs that also solved in for the Hector concentration driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiment. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-concen-rcp85'

#' Ensemble of Hector rcp26 emission driven runs with constant carbon cycle parameters
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 2.6 emissions, this
#' data set only contains results from runs that also solved in for the Hector emissions driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs. The carbon cycle parameters were held constant here.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-emissConstantC-esmrcp26'

#' Ensemble of Hector rcp45 emission driven runs with constant carbon cycle parameters
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 4.5 emissions, this
#' data set only contains results from runs that also solved in for the Hector emissions driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs. The carbon cycle parameters were held constant here.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-emissConstantC-esmrcp45'

#' Ensemble of Hector rcp60 emission driven runs with constant carbon cycle parameters
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 6.0 emissions, this
#' data set only contains results from runs that also solved in for the Hector emissions driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs. The carbon cycle parameters were held constant here.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-emissConstantC-esmrcp60'

#' Ensemble of Hector rcp85 emission driven runs with constant carbon cycle parameters
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 8.5 emissions, this
#' data set only contains results from runs that also solved in for the Hector emissions driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs. The carbon cycle parameters were held constant here.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-emissConstantC-esmrcp85'

#' Ensemble of Hector rcp26 emission driven runs with varying carbon cycle parameters
#'
#' This dataset contains outputs for the 1000 Hector runs driven with rcp 2.6 emissions, this
#' data set only contains results from runs that also solved in for the Hector emissions driven
#' rcp 4.5, rcp 6.0, and rcp 8.5 runs. The carbon cycle parameters were varied.
#'
#' @format Data frame with 9 columns
#' \describe{
#' \item{runid}{Contains integers to indicate that results are from the same Hector run}
#' \item{variable}{Variable being reported "tas", temperatures are anomalies.}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' \item{S}{Value of the Hector run climate sensitivity.}
#' \item{alpha}{Value of the Hector run aerosol scaling factor.}
#' \item{diff}{Value of the Hector run ocean heat diffusivity.}
#' }
#' @family PCA hector ensemble
'PCA_hector_ensemble-emissConstantC-esmrcp26'
