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

