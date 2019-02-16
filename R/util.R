#### Misc useful functions

#' Calculate the fraction of the variance of PCA restlts
#'
#' @param pca object returned by \code{\link[stats]{prcomp}}
#' @return Vector of of the portion of variance expalained by the PCs, should be
#' plotted with using \code{\link{plot_var}} function
#' @export
calc_variance <- function(pca){
    cumsum(pca$sdev^2)/sum(pca$sdev^2)
}
