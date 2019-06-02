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
#' \item{variable}{Variable being reported ("tas", "co2", or "heatflux").  Temperatures are
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

#' Ensembles of Hector concentration driven runs
#'
#' This dataset contains outputs for the 1000 Hector runs driven with each of
#' the four concentration pathways, along with the parameters used to drive the
#' runs.
#'
#' Temperatures (\code{tas}) temperature anomalies reported in degC.  CO2
#' concentrations (\code{co2}) are global concentrations reported in ppm.  Ocean
#' heat flux (\code{heatflux}) is reported as downward flux in W/m2
#'
#' @format List of 6 data frames: historical, rcp26, rcp45, rcp60, rcp85, and
#' params.  The four output data frames have 6 columns:
#' \describe{
#' \item{runid}{Unique identifier for each run.}
#' \item{variable}{Variable being reported: "tas" or "heatflux".}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' }
#' The params data frame has 4 columns:
#' \describe{
#' \item{runid}{Unique identifier for each run.}
#' \item{S}{Equilibrium climate sensitivity.}
#' \item{alpha}{Aerosol scaling factor.}
#' \item{volscl}{Volcanic forcing scaling factor.}
#' \item{diff}{Diffusivity parameter.}
#' }
#' @family PCA hector ensemble
'hector_conc_ensemble'


#' Ensembles of Hector emission driven runs
#'
#' This dataset contains outputs for the 1000 Hector runs driven with each of
#' the four concentration pathways, along with the parameters used to drive the
#' runs.
#'
#' Temperatures (\code{tas}) temperature anomalies reported in degC.  CO2
#' concentrations (\code{co2}) are global concentrations reported in ppm.  Ocean
#' heat flux (\code{heatflux}) is reported as downward flux in W/m2
#'
#' @format List of 6 data frames: esmHistorical, esmrcp26, esmrcp45, esmrcp60,
#' esmrcp85, and params.  The four output data frames have 6 columns:
#' \describe{
#' \item{runid}{Unique identifier for each run.}
#' \item{variable}{Variable being reported: "tas", "co2", or "heatflux".}
#' \item{scenario}{Hector scenario name.}
#' \item{experiment}{Experiment being reported, should match the cmip experiments. One of  "rcp26", "rcp45", "rcp60", "rcp85", "esmHistorical", "esmrcp85".}
#' \item{year}{Year being described (1850 -- 2100).}
#' \item{value}{Value of the variable being reported.}
#' }
#' The params data frame has 7 columns
#' \describe{
#' \item{runid}{Unique identifier for each run.}
#' \item{S}{Equilibrium climate sensitivity.}
#' \item{alpha}{Aerosol scaling factor.}
#' \item{volscl}{Volcanic forcing scaling factor.}
#' \item{diff}{Diffusivity parameter.}
#' \item{C0}{Preindustrial CO2 concentration.}
#' \item{q10_rh}{Heterotrophic respiration Q10 parameter.}
#' \item{beta}{CO2 fertilization parameter}
#' }
#' @family PCA hector ensemble
'hector_emiss_ensemble'


#' Hector concentration driven principal compoents (PCs)
#'
#' This oject contains the sd, rotation, center, and scale for the PCs calculated with
#' \code{prcomp} on an enemsble of Hector concentration driven RCP 2.6, RCP 4.5, RCP 6.0, and RCP 8.5 runs
#' along with metadata about the Hector data that went into the runs.
#'
#' @format A list with 5 elements
#' \describe{
#' \item{sdev}{A vector of the standard deviations of the principal components.}
#' \item{rotation}{A matrix of the variable loadings.}
#' \item{center}{The centering used in the PCA.}
#' \item{scale}{The scale used in the PCA.}
#' \item{meta_data}{a list of meta data information about the data that went into the PCA}
#' }
#' The metadata element is itself a list with four elements:
#' \describe{
#' \item{year}{Year used in the PCA (1850 -- 2100).}
#' \item{variable}{The variables that describe the values used in the PCA ("tas").}
#' \item{experiment}{The CMIP5 experiment names describing the data that was used in the PCA.}
#' \item{scenario}{The Hector scenario name describing the data that was used in the PCA.}
#' }
#' @family PC hector
'pc_conc'


#' Hector concentration driven principal compoents (PCs) with heat flux
#'
#' This oject contains the sd, rotation, center, and scale for the PCs calculated with
#' \code{prcomp} on an enemsble of Hector concentration driven RCP 2.6, RCP 4.5, RCP 6.0, and RCP 8.5 runs
#' along with metadata about the Hector data that went into the runs.
#'
#' @format A list with 5 elements
#' \describe{
#' \item{sdev}{A vector of the standard deviations of the principal components.}
#' \item{rotation}{A matrix of the variable loadings.}
#' \item{center}{The centering used in the PCA.}
#' \item{scale}{The scale used in the PCA.}
#' \item{meta_data}{a list of meta data information about the data that went into the PCA}
#' }
#' The metadata element is itself a list with four elements:
#' \describe{
#' \item{year}{Year used in the PCA (1850 -- 2100).}
#' \item{variable}{The variables that describe the values used in the PCA ("tas"
#' and "heatflux").}
#' \item{experiment}{The CMIP5 experiment names describing the data that was used in the PCA.}
#' \item{scenario}{The Hector scenario name describing the data that was used in the PCA.}
#' }
#' @family PC hector
'pc_conc_hflux'


#' Hector emission driven principal compoents (PCs)
#'
#' This oject contains the sd, rotation, center, and scale for the PCs calculated with
#' \code{prcomp} on an enemsble of Hector runs driven with RCP 8.5 emissions,
#' along with some meta data information about the Hector data that went into the runs.
#'
#' @format A list with 5 elements
#' \describe{
#' \item{sdev}{A vector of the standard deviations of the principal components.}
#' \item{rotation}{A matrix of the variable loadings.}
#' \item{center}{The centering used in the PCA.}
#' \item{scale}{The scale used in the PCA.}
#' \item{meta_data}{a list of meta data information about the data that went into the PCA}
#' }
#' The metadata element is itself a list with four elements:
#' \describe{
#' \item{year}{Year used in the PCA (1850 -- 2100).}
#' \item{variable}{The variables that describe the values used in the PCA
#' ("tas" and "co2")}
#' \item{experiment}{The CMIP5 experiment names describing the data that was used in the PCA.}
#' \item{scenario}{The Hector scenario name describing the data that was used in the PCA.}
#' }
#' @family PC hector
'pc_emiss'


#' Hector emission driven principal compoents (PCs) with heat flux included
#'
#' This oject contains the sd, rotation, center, and scale for the PCs calculated with
#' \code{prcomp} on an enemsble of Hector runs driven with RCP 8.5 emissions,
#' along with some meta data information about the Hector data that went into the runs.
#'
#' @format A list with 5 elements
#' \describe{
#' \item{sdev}{A vector of the standard deviations of the principal components.}
#' \item{rotation}{A matrix of the variable loadings.}
#' \item{center}{The centering used in the PCA.}
#' \item{scale}{The scale used in the PCA.}
#' \item{meta_data}{a list of meta data information about the data that went into the PCA}
#' }
#' The metadata element is itself a list with four elements:
#' \describe{
#' \item{year}{Year used in the PCA (1850 -- 2100).}
#' \item{variable}{The variables that describe the values used in the PCA
#' ("tas", "co2", and "heatflux")}
#' \item{experiment}{The CMIP5 experiment names describing the data that was used in the PCA.}
#' \item{scenario}{The Hector scenario name describing the data that was used in the PCA.}
#' }
#' @family PC hector
'pc_emiss_hflux'


#' Principal component projection coefficients for ESM outputs
#'
#' These tables give the projection coefficients onto the hector principal
#' components for each of the CMIP5 models that provided sufficient data to make
#' the projection.  The concentration-driven runs are projected onto
#' \code{\link{pc_conc}}, and the emissions-driven runs are projected onto
#' \code{\link{pc_emiss}}.
#'
#' @format A data frame with 3 columns
#' \describe{
#' \item{model}{Name of the model.}
#' \item{PC}{Index of the principal component.}
#' \item{value}{Projection coefficient.}
#' }
#' @family comparison data
#' @name esm_pcdata
NULL

#' cmip_conc_pcproj: Projection coefficients for concentration-driven ESMs
#' @rdname esm_pcdata
'cmip_conc_pcproj'

#' cmip_emiss_pcproj: Projection coefficients for emissions-driven ESMs
#' @rdname esm_pcdata
'cmip_emiss_pcproj'

#' ESM comparison data in principal components projection.
#'
#' This table gives the min-mean-max for the principal components projections of
#' the ESM ensemble.  Optionally, it \emph{may} have a single heat flux value
#' that can be used to supplement the principal components constraints.  If
#' present, this value will always be the year 2100 value for the appropriate
#' RCP 8.5 experiment (either rcp85 or esmrcp85)  Because the principal
#' components span all of the years and experiments in the comparison data,
#' there is no \code{year} or \code{experiment} column in this dataset;
#' otherwise, it is organized much like \code{\link{esm_comparison}}.
#'
#' #' @format Data frame with 7 columns
#' \describe{
#' \item{variable}{Variable being reported (PC number or "heatflux").}
#' \item{mina}{Minimum value reported in the ensemble}
#' \item{maxb}{Maximum value reported in the ensemble}
#' \item{a10}{First decile reported in the ensemble}
#' \item{b90}{Ninth decile reported in the ensemble}
#' \item{cmean}{Mean value reported in the ensemble}
#' \item{cmedian}{Median value reported in the ensemble}
#' }
#' @family comparison data
#' @name pccd
NULL

#' conc_pc_comparison: Principal components comparison data for concentration-driven runs
#' @rdname pccd
'conc_pc_comparison'

#' emiss_pc_comparison: Principal components comparison data for emissions-driven runs
#' @rdname pccd
'emiss_pc_comparison'
