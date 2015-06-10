[![Build Status](https://travis-ci.org/grunwaldlab/poppr.svg?branch=diversity_table)](https://travis-ci.org/grunwaldlab/poppr)
[![Coverage Status](http://codecov.io/github/grunwaldlab/poppr/coverage.svg?branch=diversity_table)](http://codecov.io/github/grunwaldlab/poppr?branch=diversity_table)

# Poppr version 2.0

This version of poppr is version 2.0 and is only compatible with 
[adegenet version 2.0](https://github.com/thibautjombart/adegenet). If you want 
to try out poppr version 2.0, you must install adegenet 2.0 first. See the 
Installation section below.

```R
devtools::install_github("thibautjombart/adegenet")
devtools::install_github("grunwaldlab/poppr@diversity_table")
```

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
(<img src="http://latex.codecogs.com/gif.latex?I_A" alt = "Index of association">)
or (<img src="http://latex.codecogs.com/gif.latex?%5Cbar%7Br%7D_d" alt = "Standardized index of association">)
- batch processing on any server that has R ( &ge; 2.15.1) installed
- calculation of Bruvo's distance for microsatellite (SSR) markers (implemented in C for speed)
- import of data from and export to [GenAlEx](http://biology.anu.edu.au/GenAlEx/Welcome.html "GenAlEx Homepage")

### New in version 2.0:

- handling of genomic SNP data
- custom multilocus genotype definitions
- collapse multilocus lineages by genetic distance
- calculate reticulate minimum spanning networks
- calculate index of association in a sliding window across snps
- bootstrapping of MLG diversity statistics
- and more!

For full details, see the NEWS file or type in your R console:

```R
news(Version == "1.1.5.99", package = "poppr")
```

## Citation

Please cite poppr as:

> Kamvar ZN, Tabima JF, Grünwald NJ. (2014) Poppr: an R package for genetic analysis of populations with clonal, partially clonal, and/or sexual reproduction. PeerJ 2:e281 [http://dx.doi.org/10.7717/peerj.281](http://dx.doi.org/10.7717/peerj.281)

and version 2.0 with new features as

> Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of genome-wide population genetic data with emphasis on clonality. Front. Genet. 6:208. doi: [10.3389/fgene.2015.00208](http://journal.frontiersin.org/article/10.3389/fgene.2015.00208/abstract)
  
You can obtain citation information in R by typing:

```R
citation(package = "poppr")
```

## Installation

### From CRAN

[![Downloads from Rstudio mirror per month](http://cranlogs.r-pkg.org/badges/poppr)](http://cran.r-project.org/web/packages/poppr/index.html)
[![Downloads from Rstudio mirror](http://cranlogs.r-pkg.org/badges/grand-total/poppr)](http://cran.r-project.org/web/packages/poppr/index.html)
[![CRAN version](http://www.r-pkg.org/badges/version/poppr)](http://cran.r-project.org/web/packages/poppr/index.html)

Binary versions for mac and windows are available for R &ge; 2.15.1 [**here**](http://cran.r-project.org/web/packages/poppr/index.html).

To install, make sure R is at least version 2.15.1 (the authors recommend &ge; 3.0), and in your console, type:

```R
install.packages("poppr")
```

If you want the absolute latest version of *poppr*, see about installing from github below.

***

### Stable and Development versions

To install this package from github, make sure you have the following:

- [Xcode](https://developer.apple.com/xcode/) (OSX)
    OR [Rtools](http://cran.r-project.org/bin/windows/Rtools/) (Windows)
- [devtools](https://github.com/hadley/devtools) (to install, use: `install.packages("devtools")`)

For Linux users, make sure that the function `getOption("unzip")` returns `"unzip"` or `"internal"`. If it doesn't, then run `options(unzip = "internal")`.

Now you can use the `install_github()` function:

#### For the latest stable release:    

```R
devtools::install_github(repo = "grunwaldlab/poppr@current")
library("poppr")
```

#### For the bleeding edge (development) version:

```R
devtools::install_github(repo = "grunwaldlab/poppr")
library("poppr")
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

