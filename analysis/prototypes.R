library('hectorcal')
library('ggplot2')
library('dplyr')

#### More functions that will come in handy later.

## detrend and standardize a sequence
detstd <- function(x) {
    indx <- 1:length(x)
    lin <- lm(x~indx)
    xdet <- x - predict(lin)

    xctr <- xdet - mean(xdet)
    xctr / sd(xctr)   # standardized values
}

## This thing we're calling "erf" is missing the sqrt(2) factor inside the pnorm
## call. Our calls *should* be erf((b-a)/(sqrt(2)*sig)), so dropping the sqrt(2)
## here allows us to drop it to from the args too, saving us a lot of sqrt(2)
## factors in our code that will only end up canceling out anyway.  It does,
## however, mean that this "erf" is not technically the real erf.
erf <- function(x) {2*pnorm(x)-1}
mesa <- function(x, a, b, sig) {(erf((b-x)/sig) - erf((a-x)/sig))/(2*(b-a))}

## Consensus plot, illustrating the style we will be using.  Colors are from the solarized palette
ggplot(data=filter(esm_comparison, experiment=='esmrcp85', year<2100), aes(x=year)) +
    geom_ribbon(aes(ymin=mina, ymax=maxb), color='#268bd2', fill='#268bd2', alpha=0.5) +
    geom_line(aes(y=cmean), color='#dc322f', size=1.1) +
    facet_wrap(facets=~variable, scales='free', strip.position='left',
               labeller = as_labeller(c(tas = "Temperature anomaly (C)",
                                        co2 = "CO2 concentration (ppm)"))) +
    theme_minimal(base_size=44) + ylab('') + theme(strip.placement='outside')

a <- 3.363861
b <- 7.34919
x <- seq(a-1, b+1, length.out=1000)
y <- mesa(x, a, b, 1)
pltdata <- data.frame(x=x, y=y)
ggplot(data=pltdata, aes(x=x, y=y)) +
    geom_line(size=1.1, color='#268bd2') +
    theme_minimal(base_size=44) + xlab('T') + ylab('p(T)')
