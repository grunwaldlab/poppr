\name{bruvo.dist}
\alias{bruvo.dist}
\title{Calculate the average Bruvo's Distance over all loci in a population.}
\usage{
  bruvo.dist(pop, replen = c(2))
}
\arguments{
  \item{pop}{a \code{\link{genind}} object}

  \item{replen}{a \code{vector} of \code{integers}
  indicating the length of the nucleotide repeats for each
  microsatellite locus.}
}
\value{
  a \code{distance matrix}
}
\description{
  Calculate the average Bruvo's Distance over all loci in a
  population.
}
\note{
  This function calculates bruvo's distance for non-special
  cases (ie. the ploidy and all alleles are known).
  Currently there is no way to import polyploid partial
  heterozygote data into adegenet. For Bruvo's Distance
  concerning special cases, see the package \code{polysat}.

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
\seealso{
  \code{\link{nancycats}}
}

