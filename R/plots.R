

#' Plot spaghetti plots for hector samples
#'
#' Spaghetti plots show the traces of hector outputs for a set of samples from a
#' Bayesian Monte Carlo calculation.
#'
#' @param mcrslt List of metrosamp structures from a Bayesian Monte Carlo run
#' @param nplot Number of traces to plot.
#' @param hcores List of hector cores initialized with the desired emissions.
#' @param pnames Hector names of the parameters in the sample matrix
#' @param vars Variables to plot
#' @param times Vector of times to include in the plot.  Default is 1850-2100.
#' @param alpha Alpha channel (transparency) for the sample traces.
#' @export
spaghetti_plot <- function(mcrslt, nplot, hcores, pnames,
                           vars=c(hector::GLOBAL_TEMP(), hector::ATMOSPHERIC_CO2()),
                           times=1850:2100, alpha=0.3)
{
    allsamps <- do.call(rbind, lapply(mcrslt, function(x) {x$samples}))

    sampout <- run_hector_samples(allsamps, nplot, hcores, times, pnames, vars)

    ## create a mean output dataset
    meanout <-
        dplyr::group_by(sampout, variable, year) %>%
          summarise(value = mean(value))

    ggplot2::ggplot(data=sampout, ggplot2::aes(x=year, y=value)) +
      ggplot2::geom_line(ggplot2::aes(group=scenario), color='blue', alpha=alpha) +
      ggplot2::geom_line(data=meanout, color='red', size=1.25) +
      ggplot2::facet_wrap(facets=~variable, scales='free_y')

}



#' Run hector runs for sampled parameters
#'
#' This function is used as a helper function in some of our other plotting
#' function, but it is also useful in its own right.  For example, the data
#' frame returned by this function is just what you need to plot density plots
#' of the Hector outputs
#' @param allsamps Matrix of all parameter samples
#' @param nsamp Number of samples to draw.  Must not be greater than
#' \code{nrow(allsamps)}.
#' @param hcores List of hector cores initialized with the desired emissions.
#' @param times Times to pull from the hector output.
#' @param pnames Hector internal names for the parameters in the columns of
#' \code{allsamps}.
#' @param vars Variables to retrieve.
#' @export
run_hector_samples <- function(allsamps, nsamp, hcores, times, pnames, vars)
{
    ntot <- nrow(allsamps)
    samp_indices <- sample.int(ntot, nsamp)

    if(!is.list(hcores)) {
        ## If a single core was passed, convert it into a list of length one
        hcores <- list(hcores)
    }

    ncore <- length(hcores)
    hout <-
            foreach(icore=seq(1,ncore), .combine=rbind) %dopar% {
                idx <- samp_indices[seq(icore, nsamp, ncore)]
                hcore <- hcores[[icore]]
                foreach(i=seq_along(idx), .combine=rbind) %do% {
                    irow <- idx[i]
                    runid <- i*ncore + icore
                    for(iparm in seq_along(pnames)) {
                        hector::setvar(hcore, NA, pnames[iparm], allsamps[irow,iparm],
                                       hector::getunits(pnames[iparm]))
                    }
                    hector::reset(hcore)
                    hector::run(hcore, max(times))
                    hector::fetchvars(hcore, times, vars,
                                      scenario=as.character(runid))
                }
            }

}


#' Make plots of principal components
#'
#' @param pc Principal components structure
#' @param nc Number of components to plot.  Default is to plot up to
#' \code{vartot} total variance.
#' @param vartot Cumulative variance cutoff for determining how many components
#' to plot.  Ignored if \code{nc} is specified.
#' @param rscl If true, reverse the rotation and scale factors applied by the
#' PCA. This will have the effect of displaying the variables in their natural
#' units.
#' @param swap_layout Swap the layout of the facet grid from the default.
#' Ignored if \code{rscl} is true, since then we need to have variables in
#' columns to make the scales come out right.
#' @return ggplot plot structure.
#' @export
plot_pc <- function(pc, nc=NA, vartot=0.995, rscl=FALSE, swap_layout=FALSE)
{
    if(is.na(nc)) {
        pcvar <- pc$sdev^2
        cumvar <- cumsum(pcvar)/sum(pcvar)
        nc <- min(which(cumvar >= vartot))
    }

    r <- pc$rotation[ , 1:nc, drop=FALSE]

    ## swap sign of PCs so that final value is always positive.
    rt <- t(r)
    signs <- sign(rt[ , nrow(r)])       # signs of the last value of each PC
    r <- t(rt*signs)

    if(rscl) {
        r <- r * pc$scale
        r <- r + pc$center
        facet_scale='free'
    }
    else {
        facet_scale='free_x'
    }

    if(swap_layout) {
        facet_layout <- ~var+expt
    }
    else {
        facet_layout <- ~expt+var
    }


    df <- as.data.frame(r)
    df$key <- row.names(r)

    df <- tidyr::extract(df, key, c('expt', 'var', 'year'),
                         '([:alnum:]+)\\.([:alnum:]+)\\.([:alnum:]+)')
    df$year <- as.integer(df$year)

    pltdata <- tidyr::gather(df, 'PC', 'value', -expt, -var, -year)

    ggplot2::ggplot(data=pltdata, ggplot2::aes(x=year, y=value, color=PC)) +
      ggplot2::geom_line() + ggplot2::facet_wrap(facet_layout, scales=facet_scale) +
      ggthemes::theme_solarized_2(base_size=16, light=FALSE) +
      ggthemes::scale_color_solarized()
}


#' Plot the fraction of variance captured by the PCs
#'
#' @param pca_l A list of pca structures
#' @param nc Number of components to plot (default is all)
#' @param cvthresh Cumulative variance threshold to mark with a horizontal
#' line.  If omitted, no line will be drawn.
#' @param labels Labels for the principal components being plotted
#' @return  ggplot dot plot of the fraction of the variance accounted for by
#' each PC
#' @export
plot_varfrac <- function(pca_l, nc=NA, cvthresh=NA, labels=NULL)
{
    if(inherits(pca_l, 'prcomp')) {
        ## user passed a bare prcomp strucutre
        pca_l <- list(pca_l)
    }

    if(is.null(labels)) {
        labels <-
            if(is.null(names(pca_l))) {
                as.character(seq_along(pca_l))
            }
            else {
                names(pca_l)
            }
    }

    pltdata <- dplyr::bind_rows(
        lapply(seq_along(pca_l),
               function(i) {
                   pca <- pca_l[[i]]
                   pcvar <- pca$sdev^2
                   cumvar <- cumsum(pcvar)/sum(pcvar)

                   if(!is.na(nc)) {
                       cumvar <- cumvar[1:nc]
                   }

                   data.frame(npc=seq_along(cumvar), cv=cumvar,
                              pc_set=labels[i], stringsAsFactors=FALSE)
               }))

    plt <-
      ggplot2::ggplot(data=pltdata, mapping=ggplot2::aes(x=npc, y=cv, col=pc_set)) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::labs(y = 'Cumulative variance fraction',
                    x = 'Number of PCs') +
      ggthemes::theme_solarized_2(base_size=16, light=FALSE) +
      ggthemes::scale_color_solarized()

    if(!is.na(cvthresh)) {
        plt <- plt + ggplot2::geom_hline(yintercept=cvthresh, color='lightgrey',
                                         linetype=2)
    }
    plt
}
