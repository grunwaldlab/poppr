## Test environments
* local OS X install, R 3.2.4
* ubuntu 12.04 (on travis-ci), R 3.2.4
* win-builder (devel [2015-11-28 r69714] and release [3.2.4])
* Linux version 4.0.9-boot2docker (root@b20c898782d0), 
  R-devel [2015-11-21 r69678] with address sanitizer: 
  https://hub.docker.com/r/zkamvar/poppr-devel-docker/builds/bm5c3vv9uidmf3gui8rizcm/

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


### WARNING with version 3.2.4

```
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! pdfTeX error (ext4): \pdfendlink ended up in different nesting level than \pd
fstartlink.
\AtBegShi@Output ...ipout \box \AtBeginShipoutBox 
                                                  \fi \fi 
l.197 \end{Section}
!  ==> Fatal error occurred, no output PDF file produced!
```


## Downstream dependencies
popprxl: pending
vcfR: pending