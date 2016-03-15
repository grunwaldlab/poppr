## Resubmission

This resubmission is to address a bizarre accumulation of user time from two
examples (spotted by Kurt):

```
* checking examples ... [69s/23s] OK
Examples with CPU or elapsed time > 5s
                  user system elapsed
cutoff_predictor 33.812  0.012   2.500
genclone-class   14.164  0.180   1.316
```

Behavior of these examples is normally well under a second, but I have set them
to not run. 

## Test environments

* local OS X install, R 3.2.4
* ubuntu 12.04 (on travis-ci), devel [2016-03-15 r70332] and R 3.2.3
* win-builder (devel [2016-03-14 r70331] and release [3.2.4])

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

This is indeed a valid URL

## Downstream dependencies

- popprxl: OK
- vcfR: OK