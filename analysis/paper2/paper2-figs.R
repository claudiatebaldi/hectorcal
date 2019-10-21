library(ggplot2)
library(hectorcal)
library(metrosamp)
library(tibble)
library(tidyr)
library(dplyr)
library(latex2exp)
library(ggthemes)

## Run from repository top level
datadir_conc <- file.path('analysis','mcmc','conc', 'primary')   # XXX Temporary, replace with final runs when available
mcruns_conc <- proc_mc_rslts(datadir_conc, 'hectorcal-conc')
datadir_emiss <- file.path('analysis','mcmc','emiss', 'primary')   # XXX Temporary, replace with final runs when available
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

#### Expectation values and CIs for emissions driven runs
ptbl_emiss <- bind_rows(lapply(protos, protostats, mcruns=mcruns_emiss))
ptbl_emiss_round <- bind_rows(lapply(protos, protostats, mcruns=mcruns_emiss, round=TRUE))
map_emiss <- spread(select(ptbl_emiss_round, parm, protocol, MAP), parm, MAP)
ev_emiss <- spread(select(ptbl_emiss_round, parm, protocol, EV), parm, EV)
ci_emiss <- mutate(ptbl_emiss_round, CI=paste(CI025, CI975, sep='--')) %>% select(parm, protocol, CI) %>%
    spread(protocol, CI)

Stbl_emiss <- filter(ptbl_emiss, parm=='S') %>% rename(S=EV, Smin=CI025, Smax=CI975)
Ktbl_emiss <- filter(ptbl_emiss, parm=='diff') %>% rename(diff=EV, dmin=CI025, dmax=CI975)
SKtbl_emiss <- bind_cols(Stbl_emiss, Ktbl_emiss)
plt_sk_emiss <- ggplot(data=SKtbl_emiss, aes(x=S, y=diff, color=PCA, shape=meancal, group=protocol)) + geom_point(size=2.5) +
    geom_errorbar(mapping=aes(ymin=dmin, ymax=dmax), alpha=0.5) +
    geom_errorbarh(mapping=aes(xmin=Smin, xmax=Smax), alpha=0.5) +
    xlab(TeX('$S$')) + ylab(TeX('$\\kappa$')) +
    scale_color_colorblind() +
    theme_bw(base_size = 12)

alphatbl_emiss <- filter(ptbl_emiss, parm=='alpha') %>% rename(alpha=EV, amin=CI025, amax=CI975)
voltbl_emiss <- filter(ptbl_emiss, parm=='volscl') %>% rename(volscl=EV, vmin=CI025, vmax=CI975)
vatbl_emiss <- bind_cols(alphatbl_emiss, voltbl_emiss)
plt_va_emiss <- ggplot(data=vatbl_emiss, aes(x=alpha, y=volscl, color=PCA, shape=meancal)) + geom_point(size=2.5) +
    geom_errorbar(mapping=aes(ymin=vmin, ymax=vmax), alpha=0.5) +
    geom_errorbarh(mapping=aes(xmin=amin, xmax=amax), alpha=0.5) +
    xlab(TeX('$\\alpha_a$')) + ylab(TeX('$\\alpha_v$')) +
    scale_color_colorblind() +
    theme_bw(base_size = 12)

## Geez, we really should have made a function to do these plots.
betatbl_emiss <- filter(ptbl_emiss, parm=='beta') %>% rename(beta=EV, bmin=CI025, bmax=CI975)
q10tbl_emiss <- filter(ptbl_emiss, parm=='q10_rh') %>% rename(q10=EV, qmin=CI025, qmax=CI975)
bqtbl_emiss <- bind_cols(betatbl_emiss, q10tbl_emiss)
plt_bq_emiss <- ggplot(data=bqtbl_emiss, aes(x=beta, y=q10, color=PCA, shape=meancal)) + geom_point(size=2.5) +
    geom_errorbar(mapping=aes(ymin=qmin, ymax=qmax), alpha=0.5, width=0) +
    geom_errorbarh(mapping=aes(xmin=bmin, xmax=bmax), alpha=0.5) +
    xlab(TeX('$\\beta$')) + ylab(TeX('$\\Q_{10}$')) +
    scale_color_colorblind() +
    theme_bw(base_size = 12)

## We have one odd parameter left.
c0tbl_emiss <- filter(ptbl_emiss, parm=='C0') %>% rename(C0=EV, c0min=CI025, c0max=CI975)
plt_c0_emiss <- ggplot(data=c0tbl_emiss, aes(x=factor(protocol), y=C0, color=PCA, shape=meancal)) + geom_point(size=2.5)+
    geom_errorbar(mapping=aes(ymin=c0min, ymax=c0max), alpha=0.5, width=0) +
    xlab('Protocol') + ylab(TeX('$C_0')) +
    scale_color_colorblind() +
    theme_bw(base_size = 12)

ggsave(file.path(figdir,'ev_sk_emiss.pdf'), plot=plt_sk_emiss, device='pdf', width=4, height=3, units='in')
ggsave(file.path(figdir,'ev_va_emiss.pdf'), plot=plt_va_emiss, device='pdf', width=4, height=3, units='in')
ggsave(file.path(figdir,'ev_bq_emiss.pdf'), plot=plt_bq_emiss, device='pdf', width=4, height=3, units='in')
ggsave(file.path(figdir,'ev_c0_emiss.pdf'), plot=plt_c0_emiss, device='pdf', width=4, height=3, units='in')


#### Pairs plots for emissions driven
pairplot_emiss_env_A <- pairplot(mcruns_emiss$mcobjs$`48`) + theme_bw(base_size=8) + labs(tag='A', title='Protocol 48')
pairplot_emiss_env_B <- pairplot(mcruns_emiss$mcobjs$`16`) + theme_bw(base_size=8) + labs(tag='B', title='Protocol 16')
pairplot_emiss_env_C <- pairplot(mcruns_emiss$mcobjs$`32`) + theme_bw(base_size=8) + labs(tag='C', title='Protocol 32')
pairplot_emiss_env_D <- pairplot(mcruns_emiss$mcobjs$`0`) + theme_bw(base_size=8) + labs(tag='D', title='Protocol 0')

ggsave(file.path(figdir,'pairplot_emiss_env_A.pdf'), plot=pairplot_emiss_env_A, device='pdf', width=6, height=6, units='in')
ggsave(file.path(figdir,'pairplot_emiss_env_B.pdf'), plot=pairplot_emiss_env_B, device='pdf', width=6, height=6, units='in')
ggsave(file.path(figdir,'pairplot_emiss_env_C.pdf'), plot=pairplot_emiss_env_C, device='pdf', width=6, height=6, units='in')
ggsave(file.path(figdir,'pairplot_emiss_env_D.pdf'), plot=pairplot_emiss_env_D, device='pdf', width=6, height=6, units='in')


#### Pairs plot for mean calibration runs for emissions driven
pairplot_emiss_mean_A <- pairplot(mcruns_emiss$mcobjs$`112`) + theme_bw(base_size=8) + labs(tag='A', title='Protocol 112')
pairplot_emiss_mean_B <- pairplot(mcruns_emiss$mcobjs$`80`) + theme_bw(base_size=8) + labs(tag='B', title='Protocol 80')
pairplot_emiss_mean_C <- pairplot(mcruns_emiss$mcobjs$`96`) + theme_bw(base_size=8) + labs(tag='C', title='Protocol 96')
pairplot_emiss_mean_D <- pairplot(mcruns_emiss$mcobjs$`64`) + theme_bw(base_size=8) + labs(tag='D', title='Protocol 64')

ggsave(file.path(figdir,'pairplot_emiss_mean_A.pdf'), plot=pairplot_emiss_mean_A, device='pdf', width=6, height=6, units='in')
ggsave(file.path(figdir,'pairplot_emiss_mean_B.pdf'), plot=pairplot_emiss_mean_B, device='pdf', width=6, height=6, units='in')
ggsave(file.path(figdir,'pairplot_emiss_mean_C.pdf'), plot=pairplot_emiss_mean_C, device='pdf', width=6, height=6, units='in')
ggsave(file.path(figdir,'pairplot_emiss_mean_D.pdf'), plot=pairplot_emiss_mean_D, device='pdf', width=6, height=6, units='in')

#### Plot the mesa function
fmesa <- function(x) {mesa(x, -1, 1, 0.1)}
plt_mesa <- ggplot(data.frame(x=c(-2,2)), aes(x)) + stat_function(fun=fmesa, color='MidnightBlue', size=1) +
    ylab('M(x)') +
    theme_bw(base_size = 12)
ggsave(file.path(figdir, 'mesaplot.pdf'), plot=plt_mesa, device='pdf', width=6, height=4, units='in')

#### Make the plot of the esm comparison data
colors <- c('#424242', solarized_pal()(4))
scens <- c('historical','rcp26','rcp45', 'rcp60', 'rcp85')
names(colors) <- scens
esmcmp <- filter(esm_comparison, experiment %in% scens,
                 variable=='tas',
                 (experiment=='historical' & year<2006) | year > 2005)
plt_cmpdata <-
    ggplot(data=esmcmp, aes(x=year, color=experiment, fill=experiment)) +
    geom_line(mapping=aes(y=cmean), size=1.25) +
    geom_ribbon(mapping=aes(ymin=mina, ymax=maxb), alpha=0.5, size=0.1) +
    xlab('Temperature (\u00B0C)') +
    theme_bw(base_size=12) +
    scale_color_manual(values=colors, aesthetics=c('color','fill'))
ggsave(file.path(figdir, 'cmpdata.pdf'), plot=plt_cmpdata, device='pdf', width=6, height=4, units='in')
