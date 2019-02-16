#' Plot the fraction of variance captured by the PCs
#'
#' @param frac a vector of cumulative variance values
#' @param title optional string for plot title
#' @return  ggplot dot plot of the fraction of the variance
#' @export
plot_frac <- function(frac, title = NULL){

    ggplot2::ggplot(mapping=ggplot2::aes(x=seq_along(frac), y=frac)) +
      ggplot2::geom_point(size = 1.5, col=ggthemes::solarized_pal()(1)) +
      ggplot2::labs(y = 'Cumulative variance fraction',
                    x = 'Number of PCs',
                    title = title) +
      ggthemes::theme_solarized_2(base_size=16, light=FALSE)
}


#' Plot N PCs from the object returend from prcomp
#'
#' @param pca the PCA results, an object returned by \code{\link[stats]{prcomp}}
#' @param pc_indx vector of the PC indices to plot
#' @param title optional title for the plot
#' @return ggplot line plot of the PCs
#' @importFrom foreach foreach %do%
#' @export
plot_pcs <- function(pca, pc_indx, title = NULL) {

    # For each of the pcs indicated in the pc_indx vector transform the rotation
    # value into something to plot (mulitply by the sign of the last pc value and then
    # the scale value if applicable). Save as a long data frame with information about the
    # scenario name, year, and variable.
    cfunc <- dplyr::bind_rows
    data <- foreach(i = pc_indx, .combine = cfunc) %do% {

        # Pull out the pieces of information from the prcomp object
        tibble::tibble(rotation = pca$rotation[ , i],
                       rowname = row.names(pca$rotation),
                       sign = sign(pca$rotation[ , i]),
                       scale = pca$scale) %>%
            # Parse out information from the row name
            tidyr::separate(rowname, into = c('scenario', 'variable', 'year'), sep = '_') %>%
            dplyr::mutate(variable = dplyr::if_else(variable == 1, 'Tgav', 'CO2')) %>%
            # Make sure that the year is an integer
            dplyr::mutate(year = as.integer(year)) ->
            inter_med

        # Figure out which sign to mulitply
        inter_med %>%
            dplyr::filter(year == max(year)) %>%
            dplyr::select(scenario, variable, mulitply_by = sign) ->
            sign_tibble


        # Join the unscaled_data with the sign tibble and muliply the rotation by the appropriate sign
        # and scale if need be.
        if(!isTRUE(unique(inter_med$scale))) {

            # The scale is FALSE so there is no need to multiply the rotation by the scale
            inter_med %>%
                dplyr::left_join(sign_tibble, by = c("scenario", "variable")) %>%
                dplyr::mutate(value = rotation * mulitply_by) %>%
                dplyr::mutate(pc = as.character(i))

        } else {

            # The data was scale so the rotation must be muliplied by the scale
            inter_med %>%
                dplyr::left_join(sign_tibble, by = c("scenario", "variable")) %>%
                dplyr::mutate(value = rotation * scale * mulitply_by) %>%
                dplyr::mutate(pc = as.character(i))
        }


    }


    # Basic ggplot line
    ggplot2::ggplot(data = data) +
        ggplot2::geom_line(ggplot2::aes(year, value, color = scenario, linetype = pc),
                  size = 1) +
        ggthemes::theme_solarized_2(base_size = 16, light=FALSE) +
        ggthemes::scale_color_solarized() +
        ggplot2::labs(y = 'PC', title = title)
}
