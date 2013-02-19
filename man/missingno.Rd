\name{missingno}
\alias{missingno}
\title{How to deal with missing data in a genind object.}
\usage{
  missingno(pop, missing, quiet = FALSE)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object.}

  \item{missing}{a character string: can be "zero", "mean",
  "loci", or "geno" (see \code{Details} for definitions).]}

  \item{quiet}{if \code{TRUE}, it will print to the screen
  the action performed.}
}
\value{
  a \code{\link{genind}} object.
}
\description{
  missingno gives the user four options to deal with
  missing data.
}
\note{
  \emph{"wild missingno appeared!"}
}
\section{Details}{
  The default way that functions in \code{poppr} deal with
  missing data is to simply ignore it. These methods
  provide a way to deal with systematic missing data and to
  give a wrapper for \code{adegenet}'s \code{
  \link{na.replace}} function. ALL OF THESE ARE TO BE USED
  WITH CAUTION.

  \code{"loci"} - removes all loci containing missing data
  in the entire data set.

  \code{"geno"} - removes any
  genotypes/isolates/individuals with missing data.

  \code{"mean"} - replaces all NA's with the mean of the
  alleles for the entire data set.

  \code{"zero"} or \code{"0"} - replaces all NA's with "0".
  Introduces more diversity.
}
\examples{
data(nancycats)

# Removing 3 loci with missing data.
nancy.locina <- missingno(nancycats, "loci")

# Removing 38 individuals/isolates/genotypes with missing data.
nancy.genona <- missingno(nancycats, "geno")

# Replacing all NA with "0" (see na.replace in the adegenet package).
nancy.0 <- missingno(nancycats, "0")

# Replacing all NA with the mean of each column (see na.replace in the
# adegenet package).
nancy.mean <- missingno(nancycats, "mean")
}
\seealso{
  \code{\link{na.replace}}, \code{\link{poppr}}
}

