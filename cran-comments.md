## Resubmission

This is a resubmission. 

This submission fixes two bugs:

 - Parallel code that threw an error on Fedora systems was serialized
   (https://github.com/grunwaldlab/poppr/issues/138)
 - A significant and silent user-facing bug was fixed
   (https://github.com/grunwaldlab/poppr/issues/136)

## Test environments

* local OS X install, R 3.3.3
* ubuntu 12.04 (on travis-ci), devel (2017-04-12 r72509) and R 3.3.3
* winbuilder R version 3.4.0 beta (2017-04-12 r72509)
* Rhub fedora-gcc-devel (2017-04-08 r72499)

## R CMD check results

On winbuilder:
> 
> * checking CRAN incoming feasibility ... NOTE
> Maintainer: 'Zhian N. Kamvar <zkamvar@gmail.com>'
> 
> Days since last update: 3

I have been in communication with Professor Ligges about this update.

## Solaris SPARC

There is currently an ERROR on poppr's results, but this originates in another
package. I have submitted code to fix this: 
https://github.com/thibautjombart/adegenet/pull/176

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
