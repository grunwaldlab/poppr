# Welcome

Poppr is an R package designed for analysis of populations with mixed modes of 
sexual and clonal reproduction. It is built around the framework of [adegenet's](http://adegenet.r-forge.r-project.org/)
genind object and offers the following implementations:

- clone censoring of populations at any of multiple levels of a hierarchy
- convenient counting of multilocus genotypes and sub-setting of populations with multiple levels of hierarchy
- define multilocus genotypes
- calculation of indices of genotypic diversity, evenness, richness, and rarefaction
- drawing of dendrograms with bootstrap support for Bruvo's distance
- drawing of minimum spanning networks for Bruvo's distance
- calculation of the index of association (![equation](http://latex.codecogs.com/gif.latex?I_A)) or ![equation](http://latex.codecogs.com/gif.latex?\\bar{r}_D)
- batch processing on any server that has R (>2.15.1) installed
- calculation of Bruvo's distance for microsatellite (SSR) markers (implemented in C for speed)
- import of data from and export to [GenAlEx](http://biology.anu.edu.au/GenAlEx/Welcome.html "GenAlEx Homepage")

# Installation

To install this package from github, make sure you have the following:

- [devtools](https://github.com/hadley/devtools) (to install, use: `install.packages("devtools")`)
- [Xcode](https://developer.apple.com/xcode/) (OSX)
- [Rtools](http://cran.r-project.org/bin/windows/Rtools/) (Windows)

Now you can use the `install_github()` function:

    library(devtools)
    install_github("poppr", "poppr")
    library(poppr)

You can view the manual by typing: `vignette("poppr_manual")` after installation or by clicking [here](http://grunwaldlab.cgrb.oregonstate.edu/sites/default/files/u5/poppr_manual.pdf)

Be sure to visit our [Webpage](http://grunwaldlab.cgrb.oregonstate.edu/poppr-r-package-population-genetics) with FAQs and tutorials coming soon.
	
Enjoy!
