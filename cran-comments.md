## Resubmission

This is a resubmission because I previously submitted a day ago and subsequently
found a silent user-facing bug that impacts users negatively. This update fixes
that bug.

## Test environments

* local OS X install, R 3.3.3
* ubuntu 12.04 (on travis-ci), devel (2017-04-07 r72495) and R 3.3.3
* winbuilder R version 3.4.0 beta (2017-04-07 r72496)
* Rhub fedora-gcc-devel (2017-04-08 r72499)

## R CMD check results

All good on the above platforms.

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
