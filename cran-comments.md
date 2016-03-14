## Test environments

* local OS X install, R 3.2.4
* ubuntu 12.04 (on travis-ci), devel [2016-03-14 r70328] and R 3.2.3
* win-builder (devel [2015-11-28 r69714] and release [3.2.4])

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

This is a warning that I find both on Windows and OSX with 3.2.3. 
It should be noted that I do not find this warning in either version 3.2.3
or the devel version (2016-03-14 r70328). 

## Downstream dependencies

- popprxl: OK
- vcfR: OK