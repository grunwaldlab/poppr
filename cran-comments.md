## Test environments

* local OS X install, R 3.3.2
* ubuntu 12.04 (on travis-ci), devel [2016-11-14 r71659] and R 3.3.1
* R-hub windows, R 3.3.2
* R-hub with valgrind, R 3.3.1
* R-hub with ASAN/UBSAN, R devel [2016-09-18 r71304]
* R-hub fedora, R devel [2016-10-30 r71610]


## R CMD check results

Windows gave threw a warning due to MiKTeX errors.

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
