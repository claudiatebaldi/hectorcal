This repository contains the code and data supporting the two papers in the
Hector calibration series:  

1. "Calibrating simple climate models to individual Earth system models: Lessons
   learned from calibrating Hector" 
2. "Calibrating the Hector simple climate model II:  Bayesian calibration of
   Hector parameters"
   
The repository comprises the following subdirectories:

* `R`:  Code for the `hectorcal` package.  This comprises all of the utility
  subroutines used by the analysis and calibration code.  
* `analysis`:  Code for doing analysis on the results of the calibration runs.
  This includes various experiments, plus the code for producing figures and
  tables from the two papers.
* `data-raw`:  Raw data and code for processing it.  
* `data`:  R package data.  All of this data is created as the output of code in
  the `data-raw` subdirectory.  
* `inst/scripts`:  Code for running the Bayesian Monte Carlo from Paper 2  
* `man`:  R help files for the functions in the `R` subdir and the datasets in
  the `data` subdir.  These are all accessible through the R help system when
  the `hectorcal` package is installed.
* `tests`: Tests for the subroutines in the `hectorcal` package.  These are run
  with the `test_package` function from the `testthat` package; trying to
  `source` them in an R session won't work.  
  
To reproduce the results in the papers:
1. Clone this repository.
2. From the top directory, install the package by running `R CMD INSTALL .`
3. Run the code from `data-raw`, `analysis`, or `inst/scripts` as desired.  (The
   latter has additional instructions in a README located therein.)
   
