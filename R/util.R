#### Misc useful functions

#' Calculate the fraction of the variance of PCA restlts
#'
#' @param pca object returned by \code{\link[stats]{prcomp}}
#' @return Vector of of the portion of variance expalained by the PCs.
#' @seealso \code{\link{plot_varfrac}}
#' @export
calc_variance <- function(pca){
    cumsum(pca$sdev^2)/sum(pca$sdev^2)
}


#' Convert a list of metrosamp structures to a \code{\link[coda]{mcmc.list}}
#' structure.
#'
#' @param mslist A list mof metrosamp structures
#' @param size Size of the mcmc objects in output coda structure. The results will be
#' thinned as necessary to get to this size.
#' @export
metrosamp2coda <- function(mslist, size=NA) {

    if(is.na(size)) {
        thin <- 1
    }
    else {
        thin <- ceiling(nrow(mslist[[1]]$samples) / size)
    }

    coda::mcmc.list(
        lapply(mslist, function(ms) {
                   samps <- ms$samples
                   idx <- seq(1,nrow(samps))
                   samps <- samps[idx %% thin == 0, ]
                   coda::mcmc(samps, thin=thin)
               })
        )
}

hv2esm <- c(Tgav='tas', Ca='co2', heatflux='heatflux')
esm2hv <- c(tas='Tgav', co2='Ca', heatflux='heatflux')
#' Translate between Hector and ESM variable names
#'
#' Translate hector variable names ("Tgav", etc.) to ESM names and vice versa.
#'
#' @name vartranslate
NULL

#' @describeIn vartranslate Translate hector variable names ('Tgav', etc.) to ESM names.
#' @param hvarnames Vector of Hector variable names
#' @return Translated names. Unknown names will be replaced with \code{NA}.
#' @export
hvar2esmvar <- function(hvarnames) {
    hv2esm[hvarnames]
}

#' @describeIn vartranslate Translate ESM variable names ('tas', etc.) to hector names.
#' @param esmvarnames Vector of ESM variable names
#' @export
esmvar2hvar <- function(esmvarnames) {
    esm2hv[esmvarnames]
}


#' Make a list of Hector cores
#'
#' The list of Hector cores will be passed into the \code{make_fn_minimize}
#'
#' @param inifile A vector of ini file paths that will be used to set up the Hector cores
#' @param name A vector of the Hector core scenario names, it should reflect CMIP5 experiment names.
#' @return A list of hector instances
#' @export

setup_hector_cores <- function(inifile, name) {

    # Check inputs
    assertthat::assert_that(length(inifile) == length(name),
                            msg = 'inifile and name must be vectors of equal lengths')

    mapply(hector::newcore, inifile=inifile, name=name, suppresslogging=TRUE)
}



#' Set parameters into a hector core
#'
#' Set hector parameters by name.  Parameters not mentioned in the input vector
#' will be left at their current values.  After the new parameters are set, the
#' core is reset.
#'
#' @param params A named vector of Hector parameters.
#' @param core A Hector core.
#' @return The input core, with parameters reset and spinup rerun.
#' @export
parameterize_core <- function(params, core) {

    pnames <- names(params)
    assertthat::assert_that(!is.null(pnames), msg = 'params must be named.')

    punits <- hector::getunits(pnames)
    assert_that(!any(is.na(punits)), msg='Bogus parameter names')

    mapply(function(x, y, z){hector::setvar(core = core, dates = NA, var = x, values = y, unit = z)},
           x = pnames, y = params, z = punits)

    hector::reset(core = core)
}



#### encoding and decoding runids

#' Encode and decode runids, and find matching runs.
#'
#' Functions to convert runids to serial numbers and flags, and vice versa, and
#' to load data matching specified runids or flags.
#'
#' Definition of runid:
#' \describe{
#'  \item{Bits 0-3:}{Serial number (used to give each run of the same
#' configuration a unique RNG seed}
#'  \item{Bit 4:}{Use pcs flag.  If false, we use the outputs directly.}
#'  \item{Bit 5:}{Use heat flux flag.}
#'  \item{Bit 6:}{Mean calibrtion flag}
#'  \item{Bit 7:}{Debug flag}
#' }
#'
#' @param serialnumber Vector of serial numbers.  Serial numbers must be between
#' 0 and 15, inclusive
#' @param pcsflag Vector of flags indicating whether or not to use principal
#' components.
#' @param hfflag Vector of flags indicating whether or not to use heat flux
#' @param meanflag Vector of flags indicating whether or not to do mean
#' calibration (as opposed to envelope calibration)
#' @param debugflag Vector of flags indicating whether or not to include
#' debugging information.
#' @return Vector of runid values
#' @export
encode_runid <- function(serialnumber, pcsflag, hfflag, meanflag, debugflag)
{
    assertthat::assert_that(all(serialnumber>=0 & serialnumber<16))
    assertthat::assert_that(is.logical(pcsflag))
    assertthat::assert_that(is.logical(hfflag))
    assertthat::assert_that(is.logical(meanflag))
    assertthat::assert_that(is.logical(debugflag))

    serialnumber + 16*pcsflag + 32*hfflag + 64*meanflag + 128*debugflag
}

#' @describeIn encode_runid Decode a vector of runids into flags and serial numbers
#'
#' @param runid Vector of integer runid values
#' @return List containing a vector of serial numbers and logical vectors for
#' each of the four flag types.
#' @export
decode_runid <- function(runid)
{
    runid <- as.integer(runid)
    serialnumber <- bitwAnd(runid, 15)
    pcsflag <- as.logical(bitwAnd(runid, 16))
    hfflag <- as.logical(bitwAnd(runid, 32))
    meanflag <- as.logical(bitwAnd(runid, 64))
    debugflag <- as.logical(bitwAnd(runid, 128))

    list(serialnumber=serialnumber, pcsflag=pcsflag, hfflag=hfflag,
         meanflag=meanflag, debugflag=debugflag)
}


#' Load all outputs matching specified runids
#'
#' Generally it's only useful to load runs that represent different chains for
#' the same config.  Therefore, this should usually only be called from
#' \code{\link{load_mc_output}}, which knows how to construct a sensible list of
#' runids.
#'
#' @param runid Runid values to match
#' @param dir Directory to load files from
#' @param filestem Stem of the data file names
#' @param niter Number of iterations. If NA, match all iteration counts
#' @keywords internal
#' @export
load_matching_mcout <- function(runid, dir='.', filestem='hectorcal', niter=NA)
{
    ## construct a regular expression to get the files we need
    if(length(niter) > 1) {
        nregexp <- paste0('(',
                          paste(niter, collapse='|'),
                          ')')
    } else if (is.na(niter)) {
        nregexp <- '[0-9]+'
    } else {
        nregexp <- as.character(niter)
    }

    if(length(runid) == 1) {
        runregexp <- as.character(runid)
    }
    else {
        runregexp <- paste0('(',
                            paste(runid, collapse='|'),
                            ')')
    }

    fileregexp <- paste0(filestem, '-', nregexp, '-', runregexp, '-',
                         'mcrslt\\.rds')


    files <- list.files(dir, fileregexp, full.names = TRUE)
    structure(lapply(files, readRDS), names=files)
}


#' @describeIn encode_runid Load available runs for a given configuration.
#'
#' @param dir Directory to load files from
#' @param filestem Stem of the data file names
#' @param niter Number of iterations. If NA, match all iteration counts
#' @export
load_mc_output <- function(pcsflag, hfflag, meanflag, debugflag, dir='.',
                           filestem='hectorcal', niter=NA)
{
    base_runid <- encode_runid(0, pcsflag, hfflag, meanflag, debugflag)
    runids <- seq(0,15) + base_runid

    load_matching_mcout(runids, dir, filestem, niter)
}
