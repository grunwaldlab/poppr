## Welcome

Poppr is an R package designed for analysis of populations with mixed modes of 
sexual and clonal reproduction. It is built around the framework of [adegenet's](http://adegenet.r-forge.r-project.org/)
genind object and offers the following implementations:

- clone censoring of populations at any of multiple levels of a hierarchy
- convenient counting of multilocus genotypes and sub-setting of populations with multiple levels of hierarchy
- define multilocus genotypes
- calculation of indices of genotypic diversity, evenness, richness, and rarefaction
- drawing of dendrograms with bootstrap support for Bruvo's distance
- drawing of minimum spanning networks for genetic distances
- calculation of the index of association (![equation](http://latex.codecogs.com/gif.latex?I_A)) or ![equation](http://latex.codecogs.com/gif.latex?\\bar{r}_d)
- batch processing on any server that has R ( &ge; 2.15.1) installed
- calculation of Bruvo's distance for microsatellite (SSR) markers (implemented in C for speed)
- import of data from and export to [GenAlEx](http://biology.anu.edu.au/GenAlEx/Welcome.html "GenAlEx Homepage")


[![Build Status](https://travis-ci.org/poppr/poppr.png?branch=master)](https://travis-ci.org/poppr/poppr?branch=master)

If the image above says "Passing", then that means it should be safe to install with the latest version of R. If it does not say "Passing", I am probably trying to fix whatever problem is causing it as fast as I can.

## Installation

### From CRAN

Binary versions for mac and windows are available for R &ge; 2.15.1 [**here**](http://cran.r-project.org/web/packages/poppr/index.html).

To install, make sure R is at least version 2.15.1 (the authors recommend 3.0), and in your console, type:

```s
install.packages("poppr")
```

If you want the absolute latest version of *poppr*, see about installing from github below.

***

### Stable and Development versions

To install this package from github, make sure you have the following:

- [Xcode](https://developer.apple.com/xcode/) (OSX)
    OR [Rtools](http://cran.r-project.org/bin/windows/Rtools/) (Windows)
- [devtools](https://github.com/hadley/devtools) (to install, use: `install.packages("devtools")`)

First, install the required dependencies:

```s
install.packages(c("adegenet", "pegas", "vegan", "ggplot2", "phangorn", "ape", "igraph"))
```

Now you can use the `install_github()` function:

#### For the latest stable release:    

```s
library(devtools)
install_github(repo = "poppr", user = "grunwaldlab")
library(poppr)
```

#### For the bleeding edge development version:

```s
library(devtools)
install_github(repo = "poppr", user = "grunwaldlab", ref = "devel")
library(poppr)
```

You can view the manual by typing: `vignette("poppr_manual")` after installation or by clicking [**here**](http://grunwaldlab.cgrb.oregonstate.edu/sites/default/files/u5/poppr_manual.pdf)

Be sure to visit our [Webpage](http://grunwaldlab.cgrb.oregonstate.edu/poppr-r-package-population-genetics) with FAQs and tutorials coming soon.
	
Enjoy!
