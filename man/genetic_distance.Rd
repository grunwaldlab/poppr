\name{nei.dist}
\alias{edward.dist}
\alias{nei.dist}
\alias{provesti.dist}
\alias{reynold.dist}
\alias{roger.dist}
\title{Calculate Genetic Distance for a genind or genclone object.}
\usage{
nei.dist(x, warning = TRUE)

edward.dist(x)

roger.dist(x)

reynold.dist(x)

provesti.dist(x)
}
\arguments{
  \item{x}{a \linkS4class{genind}, \linkS4class{genclone},
  or matrix object.}

  \item{warning}{If \code{TRUE}, a warning will be printed
  if any infinite values are detected and replaced. If
  \code{FALSE}, these values will be replaced without
  warning. See Details below.}
}
\value{
an object of class dist with the same number of
observations as the number of individuals in your data.
}
\description{
These functions are modified from the function
\link[adegenet]{dist.genpop} to be applicable for distances
between individuals.
}
\details{
It is important to be careful with the interpretation of
these distances as they were originally intended for
calculation of between-population distance. As Nei's
distance is the negative log of 0:1, this means that it is
very possible to obtain distances of infinity. When this
happens, infinite values are corrected to be 10 * max(D)
where D is the distance matrix without infinite values.
}
\note{
Provesti's distance is identical to
\code{\link{diss.dist}}, except that
\code{\link{diss.dist}} is optimized for a larger number of
individuals (n > 125) at the cost of required memory.
}
\examples{
data(nancycats)
nan9 <- popsub(nancycats, 9)
neinan <- nei.dist(nan9)
ednan <- edward.dist(nan9)
rodnan <- roger.dist(nan9)
reynan <- reynold.dist(nan9)
pronan <- provesti.dist(nan9)
}

