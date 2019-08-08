library(ggplot2)
library(hectorcal)
library(metrosamp)
library(tibble)
library(tidyr)
library(dplyr)
library(latex2exp)
library(ggthemes)

## Run from repository top level
datadir_conc <- file.path('analysis','mcmc','conc', 'calibration')   # XXX Temporary, replace with final runs when available
mcruns_conc <- proc_mc_rslts(datadir_conc, 'hectorcal-conc')
datadir_emiss <- file.path('analysis','mcmc','emiss', 'calibration')   # XXX Temporary, replace with final runs when available
mcruns_emiss <- proc_mc_rslts(datadir_emiss, 'hectorcal-emiss')
figdir <- file.path('analysis', 'figs-paper2')

#### Four panels in the pairs plot
pairplot_conc_env_A <- pairplot(mcruns_conc$mcobjs$`48`) + theme_bw(base_size=8) + labs(tag='A', title='Protocol 48')
pairplot_conc_env_B <- pairplot(mcruns_conc$mcobjs$`16`) + theme_bw(base_size=8) + labs(tag='B', title='Protocol 16')
pairplot_conc_env_C <- pairplot(mcruns_conc$mcobjs$`32`) + theme_bw(base_size=8) + labs(tag='C', title='Protocol 32')
pairplot_conc_env_D <- pairplot(mcruns_conc$mcobjs$`0`) + theme_bw(base_size=8) + labs(tag='D', title='Protocol 0')

ggsave(file.path(figdir,'pairplot_conc_env_A.pdf'), plot=pairplot_conc_env_A, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_env_B.pdf'), plot=pairplot_conc_env_B, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_env_C.pdf'), plot=pairplot_conc_env_C, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_env_D.pdf'), plot=pairplot_conc_env_D, device='pdf', width=4, height=4, units='in')

#### Pairs plot for mean calibration runs
pairplot_conc_mean_A <- pairplot(mcruns_conc$mcobjs$`112`) + theme_bw(base_size=8) + labs(tag='A', title='Protocol 112')
pairplot_conc_mean_B <- pairplot(mcruns_conc$mcobjs$`80`) + theme_bw(base_size=8) + labs(tag='B', title='Protocol 80')
pairplot_conc_mean_C <- pairplot(mcruns_conc$mcobjs$`96`) + theme_bw(base_size=8) + labs(tag='C', title='Protocol 96')
pairplot_conc_mean_D <- pairplot(mcruns_conc$mcobjs$`64`) + theme_bw(base_size=8) + labs(tag='D', title='Protocol 64')

ggsave(file.path(figdir,'pairplot_conc_mean_A.pdf'), plot=pairplot_conc_mean_A, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_mean_B.pdf'), plot=pairplot_conc_mean_B, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_mean_C.pdf'), plot=pairplot_conc_mean_C, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_mean_D.pdf'), plot=pairplot_conc_mean_D, device='pdf', width=4, height=4, units='in')

#### Table of expectation values and 95% CIs for all of these protocols
protostats <- function(protoid, mcruns, round=FALSE) {
    cproto <- as.character(protoid)
    config <- decode_runid(protoid)
    ev <- EV(mcruns$mcobjs[[cproto]])
    map <- MAP(mcruns$mcobjs[[cproto]])
    ci <- CI(mcruns$mcobjs[[cproto]])
    ci25 <- ci[,1]
    ci975 <- ci[,2]
    if(round) {
        ev <- signif(ev, 2)
        map <- signif(map, 2)
        ci25 <- signif(ci25, 2)
        ci975 <- signif(ci975, 2)
    }
    tibble(protocol=protoid, PCA=config$pcsflag, meancal=config$meanflag,
           parm=names(ev), MAP=map, EV=ev, CI025=ci25, CI975=ci975)
}
protos <- seq(0,7) * 16
ptbl <- bind_rows(lapply(protos, protostats, mcruns=mcruns_conc))
ptbl_round <- bind_rows(lapply(protos, protostats, mcruns=mcruns_conc, round=TRUE))
map_conc <- spread(select(ptbl_round, parm, protocol, MAP), parm, MAP)
ev_conc <- spread(select(ptbl_round, parm, protocol, EV), parm, EV)
ci_conc <- mutate(ptbl_round, CI=paste(CI025, CI975, sep='--')) %>% select(parm, protocol, CI) %>%
    spread(protocol, CI)

#### Make a plot of expectation values and credibile intervals
Stbl <- filter(ptbl, parm=='S') %>% rename(S=EV, Smin=CI025, Smax=CI975)
Ktbl <- filter(ptbl, parm=='diff') %>% rename(diff=EV, dmin=CI025, dmax=CI975)
SKtbl <- bind_cols(Stbl, Ktbl)
plt_sk <- ggplot(data=SKtbl, aes(x=S, y=diff, color=PCA, shape=meancal, group=protocol)) + geom_point(size=2.5) +
    #geom_label(aes(label=protocol), size=4, nudge_x= 0.025, nudge_y=0.025) +
    geom_errorbar(mapping=aes(ymin=dmin, ymax=dmax), alpha=0.5) +
    geom_errorbarh(mapping=aes(xmin=Smin, xmax=Smax), alpha=0.5) +
    xlab(TeX('$S$')) + ylab(TeX('$\\kappa$')) +
    scale_color_colorblind() +
    theme_bw(base_size = 12)

alphatbl <- filter(ptbl, parm=='alpha') %>% rename(alpha=EV, amin=CI025, amax=CI975)
voltbl <- filter(ptbl, parm=='volscl') %>% rename(volscl=EV, vmin=CI025, vmax=CI975)
vatbl <- bind_cols(alphatbl, voltbl)
plt_va <- ggplot(data=vatbl, aes(x=alpha, y=volscl, color=PCA, shape=meancal)) + geom_point(size=2.5) +
    geom_errorbar(mapping=aes(ymin=vmin, ymax=vmax), alpha=0.5) +
    geom_errorbarh(mapping=aes(xmin=amin, xmax=amax), alpha=0.5) +
    xlab(TeX('$\\alpha_a$')) + ylab(TeX('$\\alpha_v$')) +
    scale_color_colorblind() +
    theme_bw(base_size = 12)
ggsave(file.path(figdir,'ev_sk_conc.pdf'), plot=plt_sk, device='pdf', width=4, height=3, units='in')
ggsave(file.path(figdir,'ev_va_conc.pdf'), plot=plt_va, device='pdf', width=4, height=3, units='in')

#### Pairs plots for emissions driven
pairplot_emiss_env_A <- pairplot(mcruns_emiss$mcobjs$`48`) + theme_bw(base_size=8) + labs(tag='A', title='Protocol 48')
pairplot_emiss_env_B <- pairplot(mcruns_emiss$mcobjs$`16`) + theme_bw(base_size=8) + labs(tag='B', title='Protocol 16')
pairplot_emiss_env_C <- pairplot(mcruns_emiss$mcobjs$`32`) + theme_bw(base_size=8) + labs(tag='C', title='Protocol 32')
pairplot_emiss_env_D <- pairplot(mcruns_emiss$mcobjs$`0`) + theme_bw(base_size=8) + labs(tag='D', title='Protocol 0')

ggsave(file.path(figdir,'pairplot_conc_env_A.pdf'), plot=pairplot_conc_env_A, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_env_B.pdf'), plot=pairplot_conc_env_B, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_env_C.pdf'), plot=pairplot_conc_env_C, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_env_D.pdf'), plot=pairplot_conc_env_D, device='pdf', width=4, height=4, units='in')


#### Pairs plot for mean calibration runs for emissions driven
pairplot_emiss_mean_A <- pairplot(mcruns_emiss$mcobjs$`112`) + theme_bw(base_size=8) + labs(tag='A', title='Protocol 112')
pairplot_emiss_mean_B <- pairplot(mcruns_emiss$mcobjs$`80`) + theme_bw(base_size=8) + labs(tag='B', title='Protocol 80')
pairplot_emiss_mean_C <- pairplot(mcruns_emiss$mcobjs$`96`) + theme_bw(base_size=8) + labs(tag='C', title='Protocol 96')
pairplot_emiss_mean_D <- pairplot(mcruns_emiss$mcobjs$`64`) + theme_bw(base_size=8) + labs(tag='D', title='Protocol 64')

ggsave(file.path(figdir,'pairplot_conc_mean_A.pdf'), plot=pairplot_conc_mean_A, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_mean_B.pdf'), plot=pairplot_conc_mean_B, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_mean_C.pdf'), plot=pairplot_conc_mean_C, device='pdf', width=4, height=4, units='in')
ggsave(file.path(figdir,'pairplot_conc_mean_D.pdf'), plot=pairplot_conc_mean_D, device='pdf', width=4, height=4, units='in')
