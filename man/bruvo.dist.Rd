\name{bruvo.dist}
\alias{bruvo.dist}
\title{Calculate the average Bruvo's Distance over all loci in a population.}
\usage{
  bruvo.dist(pop, replen = 1, add = TRUE, loss = TRUE)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object}

  \item{replen}{a \code{vector} of \code{integers}
  indicating the length of the nucleotide repeats for each
  microsatellite locus.}

  \item{add}{if \code{TRUE}, genotypes with zero values
  will be treated under the genome addition model presented
  in Bruvo et al. 2004.}

  \item{loss}{if \code{TRUE}, genotypes with zero values
  will be treated under the genome loss model presented in
  Bruvo et al. 2004.}
}
\value{
  a \code{distance matrix}
}
\description{
  Calculate the average Bruvo's Distance over all loci in a
  population.
}
\note{
  The result of both \code{add = TRUE} and \code{loss =
  TRUE} is that the distance is averaged over both values.
  If both are set to \code{FALSE}, then the infinite
  alleles model is used. For genotypes with all missing
  values, the result will be NA.

  If the user does not provide a vector of appropriate
  length for \code{replen} , it will be estimated by taking
  the minimum difference among represented alleles at each
  locus. It is not recommended to rely on this estimation.
}
\examples{
# Please note that the data presented is assuming that the nancycat dataset
# contains all dinucleotide repeats, it most likely is not an accurate
# representation of the data.

# Load the nancycats dataset and construct the repeat vector.
data(nancycats)
ssr <- rep(1,9)

# Analyze the 1st population in nancycats

bruvo.dist(popsub(nancycats, 1), replen=ssr)

# View each population as a heatmap.
\dontrun{
sapply(nancycats$pop.names, function(x)
heatmap(as.matrix(bruvo.dist(popsub(nancycats, x))), symm=TRUE))
}
}
\author{
  Zhian N. Kamvar
}
\references{
  Bruvo, R., Michiels, N. K., D'Sousa, T. G., and
  Schulenberg, H. (2004) A simple method for calculation of
  microsatellite genotypes irrespective of ploidy level.
  \emph{Molecular Ecology} \strong{13}, 2101–2106.
}
\seealso{
  \code{\link{nancycats}}
}

