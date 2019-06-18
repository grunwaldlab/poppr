# Poppr version 2.8.3

This is an attempt at fixing a LaTeX issue with a straddling link [0]

[0]: https://hypatia.math.ethz.ch/pipermail/r-package-devel/2019q2/004093.html


Thank you for implementing the extra checks for poppr. I really
appreciate the hard work that goes into maintenance of CRAN. Because
of the extra checks, I was able to find old code that was poorly
written and subsequently was fixed.

This resubmission addresses the following concerns:

> If there are references describing the methods in your package, please
> add these in the Description field of your DESCRIPTION file in the form
> authors (year) doi:...
> authors (year) arXiv:...
> authors (year, ISBN:...)
> with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.

Done.

> Please do not modify the global environment (e.g. by using <<-) in your
> functions.

Done. There is currently one function that still contains this, but it is
properly scoped and does not modify the global environment.

> Please always reset the user's options, working directory, ... in your
> examples (e.g. man/getfile.Rd).

Done.

> Please ensure that your functions do not write by default or in your
> examples/vignettes/tests in the user's home filespace. That is not allow
> by CRAN policies. Please only write/save files if the user has specified
> a directory. In your examples/vignettes/tests you can write to tempdir()

Done.

> Please replace \dontrun{} by \donttest{} or if(interactive()){} in your
> Rd-file.

I have not done this and am hesitant to do so as several of the examples that
use \dontrun will have long-running code to show more in-depth examples or have
code that will fail to show more realistic scenarios for user file interactions.
Because of the long-running examples, I fear that wrapping them in \donttest
would cause my package to fail because of their runtime.

I do realize I could wrap these in if (interactive()), but then I would not be
able to illustrate the longer-running examples on my documentation website 
(i.e. https://grunwaldlab.github.io/poppr/reference/bitwise.dist.html#examples)

## Test environments

* local macOS install, R 3.6.0
* local ubuntu 16.04 install, R 3.6.0
* ubuntu 14.04 (on travis-ci), devel (2019-06-17 r76705) and R 3.6.0
* windows R Under development (unstable) 

## R CMD check results

Ubuntu Linux: OK
macOS:        OK
Windows:      OK 

## Downstream dependencies

- popprxl:  OK
- vcfR:     OK
- adegenet: OK

