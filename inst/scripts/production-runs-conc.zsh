#!/bin/zsh

#SBATCH -N 1
#SBATCH -t 60
#SBATCH -A GCAM

date


program=`Rscript -e 'cat(system.file("scripts", "production-runs-conc.R", package="hectorcal"))'`
nsamp=500
outfile="hectorcal-$nsamp"

runid=$SLURM_ARRAY_TASK_ID

echo "Run command:"
echo "source('$program'); production_run_conc($runid, $nsamp, '$outfile')"

Rscript -e "source('$program'); production_run_conc($runid, $nsamp, '$outfile')"

date


