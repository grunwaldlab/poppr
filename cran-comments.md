# Poppr version 2.8.7

This fixes the help links to other packages on CRAN. There is a warning from windows devel, but I believe it is a false positive because the warning refers to not being able to find the current package.

## Test environments

* local macOS install, R 4.0.2
* local ubuntu 20.04 install, R 4.0.2
* windows R Under development (unstable) 

## R CMD check results

Ubuntu Linux: OK
macOS:        OK
Windows:      WARNING

From https://win-builder.r-project.org/C2pr8e03v0mW/00check.log

* checking data for non-ASCII characters ... WARNING
  Error loading dataset 'Pinf':
   Error in .requirePackage(package) : 
    unable to find required package 'poppr'
  
  Error loading dataset 'Pram':
   Error in .requirePackage(package) : 
    unable to find required package 'poppr'
  
  Error loading dataset 'monpop':
   Error in .requirePackage(package) : 
    unable to find required package 'poppr'
  
  Error loading dataset 'old_Pinf':
   Error in .requirePackage(package) : 
    unable to find required package 'poppr'

## Downstream dependencies

- popprxl:  OK
- vcfR:     OK
- adegenet: OK

