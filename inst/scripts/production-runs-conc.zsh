#!/bin/zsh

#SBATCH -N 1
#SBATCH -t 60
#SBATCH -A GCAM

date

module load R/3.4.3
module load gcc/8.1.0

hostname

program=`Rscript -e 'cat(system.file("scripts", "production-runs-conc.R", package="hectorcal"))'`
nsamp=100
npc=10
outfile="hectorcal-$nsamp"

runid=$SLURM_ARRAY_TASK_ID

echo "Run command:"
echo "source('$program'); production_run_conc($runid, $nsamp, '$outfile', npc=$npc)"

Rscript -e "source('$program'); production_run_conc($runid, $nsamp, '$outfile', npc=$npc)"

date


