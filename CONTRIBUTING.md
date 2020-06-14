Contributing to *poppr*
=======================

Thank you for taking interest in contributing to *poppr*. All are welcome to participate in any way they feel is appropriate as long as the [community code of conduct](CONDUCT.md) is respected :+1:.

There are four ways you can contribute to the development of poppr: reporting a bug, fixing documentation errors, contributing new code, or commenting on issues/pull requests. Note that the majority of these activities do not require you to be proficient in R.

All new code or documentation must be approved by @zkamvar. 

Commenting on Issues / Pull Requests
------------------------------------

All are welcome to comment on issues or pull requests regardless if they are directly involved in creating them or maintaining *poppr*. Of course, please be respectful of others and follow the [code of conduct](CONDUCT.md).

Raising an Issue / Reporting a Bug
----------------------------------

The issue tracker is for both Bug Reports and Feature Requests. If you think you have found a bug, you can raise an issue by going to <https://github.com/grunwaldlab/poppr/issues> and submitting a reproducible example. 

Here are some guidelines to help make the process go smoothly:

#### Ensure your R version and poppr version are up-to-date.

To check your versions of R and poppr, you can run the following code

```r
R.version.string
packageVersion("poppr")
```

Your output will look something like this:

```r
R.version.string
#> [1] "R version 3.4.1 (2017-06-30)"
packageVersion("poppr")
#> [1] '2.5.0'
```

You can check this against the versions reported at <https://r-project.org> and <https://cran.r-project.org/package=poppr>.

#### Make sure your problem is reproducible and minimal

A minimal reproducible example (reprex) ensures that the maintainers will be able to address the issue in a timely manner. 

The easiest way of doing this is to use the [reprex package](https://cran.r-project.org/package=reprex), which has excellent documentation on how to create a reprex: <https://github.com/jennybc/reprex#what-is-a-reprex>.

Here's an example of how to do this for a poppr-specific example:

```r
library("reprex")
reprex({
  suppressPackageStartupMessages(library("poppr"))
  data(Pinf)
  Pinf.bd <- bruvo.dist(Pinf)
})
```

This will copy the rendered code to your clipboard so you can paste it directly
to the github issue like this:

``` r
suppressPackageStartupMessages(library("poppr"))
data(Pinf)
Pinf.bd <- bruvo.dist(Pinf)
#> Warning in bruvo.dist(Pinf): 
#>  -------------------------------------------------------------- 
#>                         !!! ALERT !!! 
#>  
#>  This warning will become an ERROR in future versions of poppr. 
#>  Please define your repeat lengths to avoid this error. 
#>  -------------------------------------------------------------- 
#> 
#> Repeat length vector for loci is not equal to the number of loci represented.
#> Estimating repeat lengths from data:
#>  c(2, 2, 6, 2, 2, 2, 2, 2, 3, 3, 2)
```

Please make sure that your examples run in a reasonable amount of time. Usually maintainers will have to run these examples several times in order to fix the issue. If an example takes more than 1 minute to run, please find a way to reduce the size of the example. 

Updating Documentation
----------------------

This is one of the best ways to contribute to *poppr*. If you see a spelling error or have a better way of phrasing an error message, vignette entry, or manual page, feel free to contribute!

All functions are documented using [roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html) comments above the functions themselves. If you want to edit documentation, please make the changes there and make sure you include the roxygen comment (`#'`) in front of your additions. An example of a documentation commit can be found here: <https://github.com/grunwaldlab/poppr/commit/d721c5b>. 

#### If you are not experienced with git/github

You can make changes directly on github by clicking the little pencil icon in the top right corner of the file when you view it on github. When you do that, github will create a copy of the *poppr* repository in your account so you can edit the file. When you are done editing the file, you may click on "Propose File Change" at the bottom of the page (You can also add a short commit message or leave the entries blank). From there, you will be taken the form to create a Pull Request, where you can reveiw your changes and click "Create Pull Request". Details can be found here: <https://help.github.com/articles/editing-files-in-another-user-s-repository/>


#### If you are experienced with git/github

1. fork poppr to your github account
2. clone poppr to your desktop 
3. create a new branch from main and push/sync to your repository
4. make and commit your changes
5. update the documentation with roxygen


Contributing New Code
---------------------

If you want to contribute code, you must ensure that it is accurate and follows
the tidyverse style guidelines <http://style.tidyverse.org/>. 

You can use the following checklist for your contributions:

 - [ ] I have [written tests](http://r-pkgs.had.co.nz/tests.html) to ensure my code works (here's an [example of a code change with a test](https://github.com/grunwaldlab/poppr/commit/43be276))
 - [ ] I have followed the tidyverse style guidelines
 - [ ] I have documented my changes in the [NEWS](NEWS) file
 - [ ] I have added myself as a contributor to the [DESCRIPTION](DESCRIPTION) file
 - [ ] All the checks on <https://travis-ci.org/grunwaldlab/poppr> pass
 
