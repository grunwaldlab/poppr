\name{nei.dist}
\alias{edwards.dist}
\alias{nei.dist}
\alias{provesti.dist}
\alias{reynolds.dist}
\alias{rogers.dist}
\title{Calculate Genetic Distance for a genind or genclone object.}
\usage{
nei.dist(x, warning = TRUE)

edwards.dist(x)

rogers.dist(x)

reynolds.dist(x)

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
ednan <- edwards.dist(nan9)
rodnan <- rogers.dist(nan9)
reynan <- reynolds.dist(nan9)
pronan <- provesti.dist(nan9)
}
\keyword{angular}
\keyword{coancestry}
\keyword{edwards}
\keyword{nei}
\keyword{provesti}
\keyword{reynolds}
\keyword{rodgers}
\keyword{rogers}

