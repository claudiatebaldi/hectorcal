library(ggplot2)
library(ggthemes)

## Function for generating a log-mesa function for use in plotting functions
logmesaplt <- function(a,b,sig) {
    function(x) {
        -log(mesa(x, a, b, sig))
    }
}

mesaplt <- ggplot(data=data.frame(x=c(-1.5, 1.5)), aes(x)) + stat_function(fun=logmesaplt(-1,1,0.1), size=1.1) +
    ylab('L(x)') +
    theme_bw(base_size=14)
