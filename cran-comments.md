# Poppr version 2.9.3

This update rearranges the C files to declare all non-R headers before R headers to address the openmp failures in clang13 as pointed out by Prof. Ripley. 

## Test environments

* local macOS install, R 4.1.1
* local ubuntu 20.04 install, R 4.1.1, R devel
* windows R Under development (unstable) 

## R CMD check results

Ubuntu Linux: OK
macOS:        OK
Windows:      OK


## Downstream dependencies

- popprxl:  OK
- vcfR:     OK
- adegenet: OK

