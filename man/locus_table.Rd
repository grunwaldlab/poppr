\name{locus_table}
\alias{locus_table}
\title{Create a table of summary statistics per locus.}
\usage{
locus_table(x, index = "simpson", lev = "allele", population = "ALL",
  information = TRUE)
}
\arguments{
  \item{x}{a \code{\linkS4class{genind}} or
  \code{\linkS4class{genclone}} object.}

  \item{index}{Which diversity index to use. Choices are
  \itemize{ \item \code{"simpson"} (Default) to give
  Simpson's index \item \code{"shannon"} to give the
  Shannon-Wiener index \item \code{"invsimpson"} to give
  the Inverse Simpson's index aka the Stoddard and Tayor
  index.}}

  \item{lev}{At what level do you want to analyze
  diversity? Choices are \code{"allele"} (Default) or
  \code{"genotype"}.}

  \item{population}{Select the populations to be analyzed.
  This is the parameter \code{sublist} passed on to the
  function \code{\link{popsub}}.  Defaults to
  \code{"ALL"}.}

  \item{information}{When \code{TRUE} (Default), this will
  print out a header of information to the R console.}
}
\value{
a table with 4 columns indicating the Number of
alleles/genotypes observed, Diversity index chosen, Nei's
1978 expected heterozygosity, and Evenness.
}
\description{
Create a table of summary statistics per locus.
}
\note{
This will calculate statistics for polyploids as well by
only counting observed allelic states.
}
\examples{
# Analyze locus statistics for the North American population of P. infestans.
data(Pinf)
locus_table(Pinf, population = "North America")
}
\author{
Zhian N. Kamvar
}
\references{
Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre
Legendre, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson,
Peter Solymos, M. Henry H.  Stevens, and Helene Wagner.
vegan: Community Ecology Package, 2012. R package version
2.0-5.

Niklaus J. Gr\"unwald, Stephen B. Goodwin, Michael G.
Milgroom, and William E. Fry. Analysis of genotypic
diversity data for populations of microorganisms.
Phytopathology, 93(6):738-46, 2003

J.A. Ludwig and J.F. Reynolds. Statistical Ecology. A
Primer on Methods and Computing. New York USA: John Wiley
and Sons, 1988.

E.C. Pielou. Ecological Diversity. Wiley, 1975.

J.A. Stoddart and J.F. Taylor. Genotypic diversity:
estimation and prediction in samples. Genetics,
118(4):705-11, 1988.

Masatoshi Nei. Estimation of average heterozygosity and
genetic distance from a small number of individuals.
Genetics, 89(3):583-590, 1978.

Claude Elwood Shannon. A mathematical theory of
communication. Bell Systems Technical Journal,
27:379-423,623-656, 1948
}
\seealso{
\code{\link{vegan::diversity}}, \code{\link{poppr}}
}
