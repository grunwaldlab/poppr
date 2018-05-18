# Poppr version 2.7.1

This version corrects a major documentation error and moves polysat to imports

## Test environments

* local macOS install, R 3.5.0
* ubuntu 12.04 (on travis-ci), devel (2018-03-15 r74411) and R 3.5.0
* winbuilder devel (2018-03-12 r74405) 

## R CMD check results

Ubuntu Linux: OK
macOS:        OK
Windows:      NOTE -- Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

I believe this is a false positive as it did not show up on the other systems.

## Downstream dependencies

- popprxl:  OK
- vcfR:     OK
- adegenet: OK
