# Using the MCMC batch scripts

These scripts are set up to allow you to easily launch a collection of
MCMC runs on a cluster using the Slurm resource manager.  Most of the
scripts are generic and should work on any cluster using Slurm.  The
exceptions are the allocation account name (which will be unique to
each cluster) and the module load statements, which ensure that the
correct versions of R and g++ are used on our system.

## Usage

First, install the `hectorcal` package, if you haven't already.  You
can do this by starting R and running
```
devtools::install('/path/to/package/dir')
```
This should also install any dependencies, such as the `hector`
and `metrosamp` packages.  You only have to do this once.

Once you have installed the package, you are ready to start running
Monte Carlo jobs.  To submit a job with the default settings

```
jobscript=`Rscript -e "cat(system.file('scripts/production-runs-conc.zsh', package='hectorcal'))"`
sbatch -a 0 $jobscript
```

Your job will be submitted to your cluster's job queue, and it will
run when resources are available.  In addition to a file of console
output, the job should produce a `.png` file and a `.rds` file.  The
former is a diagnostic plot that should allow you to tell at a glance
whether the job ran correctly.  The latter is a saved copy of the R
object returned by the Monte Carlo run.  You can load this into a
future R session to do additional analysis with it.

## Settings

The default settings are meant for testing whether the job
configuration works correctly; they don't produce a useful Monte Carlo
run.  To produce useful output you need to change the settings.

### The "runid"

The number after `-a` in the command above is the "runid".  The
meaning of the runid is given in the documentation for the
`decode_runid` function.  Briefly, the first 4 bits are a serial
number, and the remaining bits set or clear various flags used in
setting up the log-posterior for the Bayesian Monte Carlo.  The `-a`
option actually activates Slurm's batch array function, so you can
specify a list of runids, and Slurm will generate a job for each of
them.  The runids can be specified as a list (e.g., `-a 1,2,13,23`) or
as a range (e.g., `-a 0-127`).

### Environment variables

Several environment variables can be used to change settings in the
job:

* `hectorcalNAME` : Change the stem used to generate the filename for
  the output files.
  
* `hectorcalNSAMP` : Change the number of samples that the Monte Carlo
  will be run for.  This number will be encoded in the output
  filenames.

* `hectorcalNPC` : Number of principal components to keep in runs that
  use principal components (has no effect in runs that compare to the
  raw outputs).  This number is _not_ automatically encoded in the
  output file names, so consider setting `hectorcalNAME` accordingly.
  
To set an environment variable in bash (or bash-derived shells), use
(for example) `export hectorcalNSAMP=10000`.  If you're using a
C-shell variant, use (again, for example) `setenv hectorcalNAME
myruns`.

### Arguments to sbatch

There are a plethora of arguments to sbatch that you can use to
customize the job.  The most commonly useful ones are

* `-t` Set the time limit for the job.  The default time limit
  configured in the job script is short so as to allow the test jobs
  to be run on the test and debug queue on our system.  You will
  probably need something longer in order to do a production run.
  
* `-A` Set the account name.  Unless you happen to be a member of an
  account called "GCAM" on your system, you will need to reset this.
  
