# Poppr version 2.7.1

This version corrects a major documentation error and moves polysat to imports

## Test environments

* local macOS install, R 3.5.0
* ubuntu 12.04 (on travis-ci), devel (2018-03-15 r74411) and R 3.5.0
* windows R Under development (unstable) (2018-05-15 r74727)

## R CMD check results

Ubuntu Linux: OK
macOS:        OK
Windows:      OK 

I believe this is a false positive as it did not show up on the other systems.

## Downstream dependencies

- popprxl:  OK
- vcfR:     OK
- adegenet: OK
