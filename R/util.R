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
#' @export
metrosamp2coda <- function(mslist) {

    coda::mcmc.list(
        lapply(mslist, function(ms) {
                   coda::mcmc(ms$samples)
               })
        )
}

hv2esm <- c(Tgav='tas', Ca='co2')
esm2hv <- c(tas='Tgav', co2='Ca')
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
parameterize_core <- function(params, core) {

    pnames <- names(params)
    assertthat::assert_that(!is.null(pnames), msg = 'params must be named.')

    punits <- hector::getunits(pnames)
    assert_that(!any(is.na(punits)), msg='Bogus parameter names')

    mapply(function(x, y, z){hector::setvar(core = core, dates = NA, var = x, values = y, unit = z)},
           x = pnames, y = params, z = punits)

    hector::reset(core = core)
}



