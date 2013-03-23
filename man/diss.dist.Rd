\name{diss.dist}
\alias{diss.dist}
\title{Calculate a distance matrix based on relative dissimilarity}
\usage{
  diss.dist(pop)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object.}
}
\value{
  A distance object.
}
\description{
  diss.dist uses the same discreet dissimilarity matrix
  utilized by the index of association (see
  \code{\link{ia}} for details). It returns a distance
  reflecting a ratio of the number of observed differences
  by the number of possible differences. Eg. two
  individuals who share half of the same alleles will have
  a distance of 0.5. This function can analyze distances
  for any marker system.
}
\examples{
# A simple example. Let's analyze the mean distance among populations of A.
# euteiches.

data(Aeut)
mean(diss.dist(popsub(Aeut, 1)))
mean(diss.dist(popsub(Aeut, 2)))
mean(diss.dist(Aeut))
}

