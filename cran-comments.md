## Test environments

* local OS X install, R 3.3.0
* ubuntu 12.04 (on travis-ci), devel [2016-06-08 r70732] and R 3.3.0
* win-builder (devel [2016-06-08 r70728])

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
