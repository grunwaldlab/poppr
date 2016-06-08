## Test environments

* local OS X install, R 3.3.0
* ubuntu 12.04 (on travis-ci), devel [2016-06-08 r70732] and R 3.3.0
* win-builder (devel [2016-06-06 r70718] and release [3.2.4])

## R CMD check results

### NOTE on win-builder devel:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Zhian N. Kamvar <kamvarz@science.oregonstate.edu>'

Found the following (possibly) invalid URLs:
  URL: https://developer.apple.com/xcode
    From: README.md
    Status: 401
    Message: Unauthorized
```

### WARNING on win-builder devel:

```
* checking re-building of vignette outputs ... WARNING
Error in re-building vignettes:
  ...
Loading required package: adegenet
Loading required package: ade4

   /// adegenet 2.0.1 is loaded ////////////

   > overview: '?adegenet'
   > tutorials/doc/questions: 'adegenetWeb()' 
   > bug reports/feature requests: adegenetIssues()


This is poppr version 2.1.1.99.67. To get started, type package?poppr
OMP parallel support: available

This version of poppr is under development.
If you find any bugs, please report them at https://github.com/grunwaldlab/poppr/issues
Loading required package: ape
Quitting from lines 287-295 (mlg.Rmd) 
Error: processing vignette 'mlg.Rmd' failed with diagnostics:
invalid gray level, must be in [0,1].
Execution halted
```

This is indeed a valid URL

## Downstream dependencies

- popprxl: OK
- vcfR: OK
