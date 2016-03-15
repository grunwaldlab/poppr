## Resubmission

This resubmission is to address a failure building the package documentation,
The location of the error was idenitified by Kurt Hornik, and I believe it was 
caused by an errant '\#' in a URL. 

This has been fixed and re-tested.

## Test environments

* local OS X install, R 3.2.4
* ubuntu 12.04 (on travis-ci), devel [2016-03-14 r70328] and R 3.2.3
* win-builder (devel [2016-03-14 r70328] and release [3.2.4])

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