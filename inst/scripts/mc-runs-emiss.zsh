#!/bin/zsh

#SBATCH -N 1
#SBATCH -t 60
#SBATCH -A GCAM

echo "start:  " `date`

module load R/3.4.3
module load gcc/8.1.0

hostname

## Check to see if the variable hectorcalNPC is a number
if [[ $hectorcalNPC =~ [0-9]+ ]]; then
    npc=$hectorcalNPC
else
    npc=10
fi

## Check to see if hectorcalNSAMP is a number
if [[ $hectorcalNSAMP =~ [0-9]+ ]]; then
    nsamp=$hectorcalNSAMP
else
    nsamp=100
fi


## Check to see if hectorcalNAME is not empty
if [[ -n $hectorcalNAME ]]; then
    ofname=$hectorcalNAME
else
    ofname="hectorcal-emiss"
fi

program=`Rscript -e 'cat(system.file("scripts", "mc-runs-emiss.R", package="hectorcal"))'`
outfile="$ofname-$nsamp"

runid=$SLURM_ARRAY_TASK_ID

if [[ -n $hectorcalRESTART ]]; then
    echo "Run command:"
    echo "source('$program'); mc_run_emiss($runid, $nsamp, '$outfile', npc=$npc, restart='$hectorcalRESTART')"

    Rscript -e "source('$program'); mc_run_emiss($runid, $nsamp, '$outfile', npc=$npc, restart='$hectorcalRESTART')"
else
    echo "Run command:"
    echo "source('$program'); mc_run_emiss($runid, $nsamp, '$outfile', npc=$npc)"

    Rscript -e "source('$program'); mc_run_emiss($runid, $nsamp, '$outfile', npc=$npc)"
fi

echo "end:  "  `date`


