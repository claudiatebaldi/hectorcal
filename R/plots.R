

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
          dplyr::summarise(value = mean(value))

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

#' Plot the first two PC coordinates for a table of ESM projections
#'
#' @param pctable Data frame containing the ESM projections by model
#' @importFrom dplyr %>%
#' @export
esm_pcplot <- function(pctable)
{
    pltdata <- dplyr::filter(pctable, PC %in% c(1,2)) %>%
        tidyr::spread('PC','value')
    names(pltdata) <- c('model','PC1','PC2')
    xmid <- 0.5*(min(pltdata$PC1) + max(pltdata$PC1))
    xrng <- 1.1*(max(pltdata$PC1) - min(pltdata$PC1))
    xlo <- xmid - 0.5*xrng
    xhi <- xmid + 0.5*xrng

    ggplot2::ggplot(data=pltdata, ggplot2::aes(x=PC1, y=PC2, label=model)) +
        ggplot2::geom_label() +
        ggplot2::xlim(c(xlo,xhi)) +
        ggthemes::theme_solarized_2(light=FALSE, base_size=14)
}

#' Plot the rotation of the principal components into the output space.
#'
#' Show what a gate defined in the principal components space looks like in the
#' output space.
#'
#' Six principal components can be used in the plot.  Two of them will be used
#' to draw the boxes in the output space.  The other four will be assigned to
#' the color, linetype, facet row, and facet column aesthetics.  These
#' assignments are referred to as x, y, col, lty, facetrow, and facetcol,
#' respectively.  By default they are assigned to PCs 1-6 in that order.  This
#' assignment can be changed by passing a named vector to the pccords
#' parameter.  The names give the roles to rebind, and the (integer) values
#' give the PC to bind to them.  You only need to pass the ones you want to change
#' from their defaults, so for example, passing \code{c(x=3, lty=1)} would swap
#' the roles of PC1 and PC3, while leaving the other roles unchanged.
#'
#' Note that if any of the PCs used in pccords are not in the pclimits table, or
#' if any of them occur more than once, you will get an error.
#'
#' @param pclimits Data frame giving the high, medium, and low values for each
#' PC.
#' @param pcstruct Principal components structure.
#' @param yrx Year to plot as the x-value in the output space
#' @param yry Year to plot as tye y-value in the output space
#' @param var Variable to plot in the output space
#' @param expt Experiment to use for the output variable assignment.
#' @param pccords Assignment of coordinates to principal components (see details)
#' @param HML Names of the columns in pclimits containing the High, Medium, and
#' Low values for the gates.  (Also, my grandmother's initials.  Miss you.)
#' @importFrom foreach foreach %do%
#' @importFrom dplyr %>%
#' @export
pc_rotplot <- function(pclimits, pcstruct, yrx=2006, yry=2100, var='tas', expt='rcp85',
                       pccords=NULL, HML=c('max', 'mean', 'min'))
{
    year <- variable <- NULL

    npc <- 6
    pccords_default <- seq(1,npc)
    names(pccords_default) <- c('x','y','col','lty','facetrow', 'facetcol')
    if(is.null(pccords)) {
        pccords <- pccords_default
    }
    for(n in names(pccords_default)) {
        if(!(n %in% names(pccords))) {
            pccords[[n]] <- pccords_default[[n]] # Use double braket for all
                                        # this stuff so that if someone passes
                                        # in a list instead of a vector it still
                                        # works.
        }
    }

    ## combinations of x and y roles to use to draw a box.  The indices refer to
    ## the entries in the HML vector
    boxx <- c(3,3,3, 2, 1,1,1, 2)
    boxy <- c(3,2,1, 1, 1,2,3, 3)

    ## Display names of the levels for facets
    levelnames <- c('zhi', 'med', 'low')
    ## Rejigger the names for linetype so that "medium" gets the solid line.  I
    ## should really figure out how to do this the right way.
    ltylevelnames <- c('zhi', 'amed', 'low')

    ## Discard any tibble attributes that pclimits might have
    pclimits <- as.data.frame(pclimits)

    plotcolnames <- paste0('PC', pccords[names(pccords_default)])

    pltdata <-
        foreach (col3 = 1:3, .combine=dplyr::bind_rows) %do% {
            x3 <- pclimits[pclimits$PC == pccords[['col']], HML[[col3]]]
            foreach (col4 = 1:3, .combine=dplyr::bind_rows) %do% {
                x4 <- pclimits[pclimits$PC == pccords[['lty']], HML[[col4]]]
                foreach (col5 = 1:3, .combine=dplyr::bind_rows) %do% {
                    x5 <- pclimits[pclimits$PC == pccords[['facetrow']], HML[[col5]]]
                    foreach(col6 = 1:3, .combine=dplyr::bind_rows) %do% {
                        x6 <- pclimits[pclimits$PC == pccords[['facetcol']], HML[[col6]]]

                        foreach(i = seq_along(boxx), .combine=dplyr::bind_rows) %do% {
                            x1 <- pclimits[pclimits$PC == pccords[['x']],
                                           HML[[boxx[i]]]]
                            x2 <- pclimits[pclimits$PC == pccords[['y']],
                                           HML[[boxy[i]]]]
                            pcvec <- rep(0, npc)
                            ## The names of pccords_default are guaranteed to be
                            ## in the right order.  pccords, not so much.
                            pcvec[pccords[names(pccords_default)]] <-
                                c(x1,x2,x3,x4,x5,x6)

                            cd <- reconstruct_climate(pcvec, pcstruct, npc) %>%
                              dplyr::filter(year %in% c(yrx, yry),
                                            variable==var, experiment==expt) %>%
                              as.data.frame()

                            rslt <-
                                data.frame(v1 = cd[cd$year==yrx, 'value'],
                                           v2 = cd[cd$year==yry, 'value'],
                                           v3 = levelnames[col3],
                                           v4 = ltylevelnames[col4],
                                           v5 = levelnames[col5],
                                           v6 = levelnames[col6], stringsAsFactors = FALSE)
                            names(rslt) <- plotcolnames

                            rslt
                        }
                    }
                }
            }
        }

    ggplot2::ggplot(data=pltdata,
                    mapping=
                      ggplot2::aes_string(plotcolnames[1], plotcolnames[2],
                                          color=plotcolnames[3],
                                          linetype=plotcolnames[4])) +
      ggplot2::geom_path() +
      ggplot2::facet_grid(paste(plotcolnames[5], '~', plotcolnames[6])) +
      xlab(paste(var,'(', yrx, ')')) + ylab(paste(var,'(', yry, ')'))
}
