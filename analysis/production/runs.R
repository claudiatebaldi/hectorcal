funcfile <- system.file('scripts/production-runs-conc.R', package='hectorcal')
source(funcfile)
args <- commandArgs(trailingOnly=TRUE)
strt <- as.integer(args[1])
runids <- 128 + seq(strt, 127, 4)

for(runid in runids) {
    message('runid: ', runid)
    production_run_conc(runid, 10000)
}
