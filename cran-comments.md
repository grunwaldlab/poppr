## Test environments

* local OS X install, R 3.3.1
* ubuntu 12.04 (on travis-ci), devel [2016-08-28 r71162] and R 3.3.1
* win-builder (devel [2016-08-27 r71160])

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

## valgrind

There is currently a valgrind note on CRAN. It was brought up on 2016-06-08.
I responded on 2016-06-09:

> It was my understanding that this was an artifact due to the behavior of
> OpenMP and valgrind (possibly summarized in a similar situation here:
> https://gcc.gnu.org/bugzilla/show_bug.cgi?id=36298), which is why the leak
> summary does not show any memory that is definitely or indirectly lost.
> 
> FWIW, the summary points to line 255 in the code, which is in the middle of an
> OMP pragma statement:
> https://github.com/grunwaldlab/poppr/blob/fda8606ab9c9965
> 42cabc8c9092e024ceb86e4cd/src/bitwise_distance.c#L254-L259

Kurt's response on 2016-06-13:

> Right.  So there's perhaps nothing you can do about this.
> 
> On its way to CRAN ...
> 
> Best
> -k

## Downstream dependencies

- popprxl: OK
- vcfR: OK
