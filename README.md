## Welcome

Poppr is an R package designed for analysis of populations with mixed modes of 
sexual and clonal reproduction. It is built around the framework of [adegenet's](http://adegenet.r-forge.r-project.org/)
genind object and offers the following implementations:

- clone censoring of populations at any of multiple levels of a hierarchy
- convenient counting of multilocus genotypes and sub-setting of populations with multiple levels of hierarchy
- define multilocus genotypes
- calculation of indices of genotypic diversity, evenness, richness, and rarefaction
- drawing of dendrograms with bootstrap support for genetic distances
- drawing of minimum spanning networks for genetic distances
- calculation of the index of association 
(![equation](http://latex.codecogs.com/gif.latex?I_A))
or (![equation](http://latex.codecogs.com/gif.latex?%5Cbar%7Br%7D_d))
- batch processing on any server that has R ( &ge; 2.15.1) installed
- calculation of Bruvo's distance for microsatellite (SSR) markers (implemented in C for speed)
- import of data from and export to [GenAlEx](http://biology.anu.edu.au/GenAlEx/Welcome.html "GenAlEx Homepage")

## Citation

Please cite poppr as:

> Kamvar ZN, Tabima JF, Grünwald NJ. (2014) Poppr: an R package for genetic analysis of populations with clonal, partially clonal, and/or sexual reproduction. PeerJ 2:e281 [http://dx.doi.org/10.7717/peerj.281](http://dx.doi.org/10.7717/peerj.281)
  
You can obtain citation information in R by typing:

```s
citation(package = "poppr")
```

## Installation

### From CRAN

Binary versions for mac and windows are available for R &ge; 2.15.1 [**here**](http://cran.r-project.org/web/packages/poppr/index.html).

To install, make sure R is at least version 2.15.1 (the authors recommend &ge; 3.0), and in your console, type:

```s
install.packages("poppr")
```

If you want the absolute latest version of *poppr*, see about installing from github below.

***

### Stable and Development versions

[![Build Status](https://travis-ci.org/grunwaldlab/poppr.png?branch=master)](https://travis-ci.org/grunwaldlab/poppr?branch=master)

If the image above says "Passing", then that means it should be safe to install with the latest version of R. If it does not say "Passing", I am probably trying to fix whatever problem is causing it as fast as I can.

To install this package from github, make sure you have the following:

- [Xcode](https://developer.apple.com/xcode/) (OSX)
    OR [Rtools](http://cran.r-project.org/bin/windows/Rtools/) (Windows)
- [devtools](https://github.com/hadley/devtools) (to install, use: `install.packages("devtools")`)

For Linux users, make sure that the function `getOption("unzip")` returns `"unzip"` or `"internal"`. If it doesn't, then run `options(unzip = "internal")`.

Now you can use the `install_github()` function:

#### For the latest stable release:    

```s
library(devtools)
install_github(repo = "grunwaldlab/poppr")
library(poppr)
```

#### For the bleeding edge (development) version:

```s
library(devtools)
install_github(repo = "grunwaldlab/poppr", ref = "devel")
library(poppr)
```

## Help / Documentation

### User Group

Users who have any questions/comments/suggestions regarding any version of poppr (stable or development) should direct their comments to the [Poppr google group](http://groups.google.com/group/poppr)

### Vignettes

Two vignettes have been written for poppr:

1. Data Import and Manipulation (`vignette("poppr_manual", package = "poppr")`)
2. Algorithms and Equation Utilized (`vignette("algo", package = "poppr")`)

### Book/Primer

In Spring of 2014, Dr. Niklaus J. Grünwald, Dr. Sydney E. Everhart and Zhian N. Kamvar wrote a primer for population genetic analysis in R located at http://grunwaldlab.github.io/Population_Genetics_in_R.

