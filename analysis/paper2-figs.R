library(ggplot2)
library(hectorcal)

## Run from repository top level
datadir <- file.path('analysis','mcmc-output','conc')
mcruns <- proc_mc_rslts(datadir)

## Four panels in the pairs plot
pairplot_A <- pairplot(mcruns$mcobjs$`48`) + theme_bw(base_size=10) + labs(tag='A')
pairplot_B <- pairplot(mcruns$mcobjs$`16`) + theme_bw(base_size=10) + labs(tag='B')
pairplot_C <- pairplot(mcruns$mcobjs$`32`) + theme_bw(base_size=10) + labs(tag='C')
pairplot_D <- pairplot(mcruns$mcobjs$`0`) + theme_bw(base_size=10) + labs(tag='D')

figdir <- file.path('analysis', 'figs-paper2')
ggsave(file.path(figdir,'pairplot_A.pdf'), plot=pairplot_A, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_B.pdf'), plot=pairplot_B, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_C.pdf'), plot=pairplot_C, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_D.pdf'), plot=pairplot_D, device='pdf', width=4, height=4, units='in')
