# Poppr version 2.6.1

This version fixes a heap-buffer-overflow in clang-ASAN. Details of the fix are documented at https://github.com/grunwaldlab/poppr/issues/169

## Test environments

* local OS X install, R 3.4.3
* ubuntu 12.04 (on travis-ci), devel (2018-01-13 r74115) and R 3.4.3
* winbuilder R-devel (2018-01-13 r74113) 

## R CMD check results

## clang-ASAN

The not currently on cran is fixed in this version.

## valgrind

There is currently a valgrind note on CRAN. I believe this fix addresses all the locations that were originally flagged in this comment.

It was brought up on 2016-06-08. I responded on 2016-06-09:

> It was my understanding that this was an artifact due to the behavior of OpenMP and valgrind (possibly summarized in a similar situation here: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=36298), which is why the leak summary does not show any memory that is definitely or indirectly lost.
> 
> FWIW, the summary points to line 255 in the code, which is in the middle of an OMP pragma statement: https://github.com/grunwaldlab/poppr/blob/fda8606ab9c996542cabc8c9092e024ceb86e4cd/src/bitwise_distance.c#L254-L259

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
