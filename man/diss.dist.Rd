\name{diss.dist}
\alias{diss.dist}
\title{Calculate a distance matrix based on relative dissimilarity}
\usage{
diss.dist(x, diff = TRUE, frac = TRUE, mat = FALSE)
}
\arguments{
  \item{x}{a \code{\link{genind}} object.}

  \item{diff}{\code{logical}. If \code{TRUE}, this will
  count the number of differences between individuals.
  \code{FALSE} will count the number of similarities.
  Default set to \code{TRUE}}

  \item{frac}{\code{logical}. Should the distance be
  represented from 0 to 1? Default set to \code{TRUE}.
  \code{FALSE} will return the distance represented as
  integers from 1 to n where n is the number of
  loci*ploidy.}

  \item{mat}{\code{logical}. Return a matrix object.
  Default set to \code{FALSE}, returning a dist object.
  \code{TRUE} returns a matrix object.}
}
\value{
Pairwise distances between individuals present in the
genind object.
}
\description{
diss.dist uses the same discreet dissimilarity matrix
utilized by the index of association (see \code{\link{ia}}
for details). It returns a distance reflecting a ratio of
the number of observed differences by the number of
possible differences. Eg. two individuals who share half of
the same alleles will have a distance of 0.5. This function
can analyze distances for any marker system.
}
\note{
This is exactly the same as \code{\link{provesti.dist}},
except that it performs better for large numbers of
individuals (n > 125).
}
\examples{
# A simple example. Let's analyze the mean distance among populations of A.
# euteiches.

data(Aeut)
mean(diss.dist(popsub(Aeut, 1)))
\dontrun{
mean(diss.dist(popsub(Aeut, 2)))
mean(diss.dist(Aeut))
}
}
\author{
Zhian N. Kamvar
}

