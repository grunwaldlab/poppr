% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bruvo.r
\name{fix_replen}
\alias{fix_replen}
\title{Find and fix inconsistent repeat lengths}
\usage{
fix_replen(gid, replen, e = 1e-05, fix_some = TRUE)
}
\arguments{
\item{gid}{a \code{\link[adegenet]{genind}} or \code{\link{genclone}} object}

\item{replen}{a numeric vector of repeat motif lengths.}

\item{e}{a number to be subtracted or added to inconsistent repeat lengths to
allow for proper rounding.}

\item{fix_some}{if \code{TRUE} (default), when there are inconsistent repeat 
lengths that cannot be fixed by subtracting or adding e, those than can be 
fixed will. If \code{FALSE}, the original repeat lengths will not be fixed.}
}
\value{
a numeric vector of corrected repeat motif lengths.
}
\description{
Attempts to fix inconsistent repeat lengths found by \code{test_replen}
}
\details{
This function is modified from the version used in 
  \doi{10.5281/zenodo.13007}.\cr Before being fed into the
  algorithm to calculate Bruvo's distance, the amplicon length is divided by
  the repeat unit length. Because of the amplified primer sequence attached
  to sequence repeat, this division does not always result in an integer and
  so the resulting numbers are rounded. The rounding also protects against
  slight mis-calls of alleles. Because we know that \deqn{\frac{(A - e) - (B
  - e)}{r}}{((A - e) - (B - e))/r} is equivalent to \deqn{\frac{A - B}{r}}{(A
  - B)/r}, we know that the primer sequence will not alter the relationships
  between the alleles. Unfortunately for nucleotide repeats that have powers
  of 2, rounding in R is based off of the IEC 60559 standard (see
  \code{\link{round}}), that means that any number ending in 5 is rounded to
  the nearest \emph{even} digit. This function will attempt to alleviate this
  problem by adding a very small amount to the repeat length so that division
  will not result in a 0.5. If this fails, the same amount will be
  subtracted. If neither of these work, a warning will be issued and it is up
  to the user to determine if the fault is in the allele calls or the repeat
  lengths.
}
\examples{

data(Pram)
(Pram_replen <- setNames(c(3, 2, 4, 4, 4), locNames(Pram)))
fix_replen(Pram, Pram_replen)
# Let's start with an example of a tetranucleotide repeat motif and imagine
# that there are twenty alleles all 1 step apart:
(x <- 1:20L * 4L)
# These are the true lengths of the different alleles. Now, let's add the
# primer sequence to them. 
(PxP <- x + 21 + 21)
# Now we make sure that x / 4 is equal to 1:20, which we know each have
# 1 difference.
x/4
# Now, we divide the sequence with the primers by 4 and see what happens.
(PxPc <- PxP/4)
(PxPcr <- round(PxPc))
diff(PxPcr) # we expect all 1s

# Let's try that again by subtracting a tiny amount from 4
(PxPc <- PxP/(4 - 1e-5))
(PxPcr <- round(PxPc))
diff(PxPcr)
}
\references{
Zhian N. Kamvar, Meg M. Larsen, Alan M. Kanaskie, Everett M. 
  Hansen, & Niklaus J. Grünwald. Sudden_Oak_Death_in_Oregon_Forests: Spatial
  and temporal population dynamics of the sudden oak death epidemic in Oregon
  Forests. ZENODO, \doi{10.5281/zenodo.13007}, 2014.
  
  Kamvar, Z. N., Larsen, M. M., Kanaskie, A. M., Hansen, E. M., & Grünwald,
  N. J. (2015). Spatial and temporal analysis of populations of the sudden
  oak death pathogen in Oregon forests. Phytopathology 105:982-989.
  doi: \doi{10.1094/PHYTO-12-14-0350-FI}
  
  Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza, and Hinrich 
  Schulenburg. A simple method for the calculation of microsatellite genotype
  distances irrespective of ploidy level. Molecular Ecology, 13(7):2101-2106,
  2004.
}
\seealso{
\code{\link{test_replen}} \code{\link{bruvo.dist}}
  \code{\link{bruvo.msn}} \code{\link{bruvo.boot}}
}
\author{
Zhian N. Kamvar
}
