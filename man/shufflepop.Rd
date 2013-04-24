\name{shufflepop}
\alias{shufflepop}
\title{Shuffle individuals in a \code{\link{genind}} object independently over each
locus.}
\usage{
  shufflepop(pop, method = 1)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object}

  \item{method}{an integer between 1 and 4. See details
  below.}
}
\value{
  a \code{\link{genind}} object shuffled by a specified
  method
}
\description{
  Shuffle individuals in a \code{\link{genind}} object
  independently over each locus.
}
\section{Details}{
  This function will shuffle each locus in the data set
  independently of one another, rendering them essentially
  unlinked. The following methods are available to shuffle
  your data: \enumerate{ \item \strong{Multilocus Style
  Permutation} This will shuffle the genotypes at each
  locus, maintaining the heterozygosity and allelic
  structure.  \item \strong{Permute Alleles} This will
  redistribute all alleles in the sample throughout the
  locus. Missing data is fixed in place. This maintains
  allelic structure, but heterozygosity is variable.  \item
  \strong{Parametric Bootstrap} This will redistribute
  available alleles within the locus based on their allelic
  frequencies. This means that both the allelic state and
  heterozygosity will vary. The resulting data set will not
  have missing data.  \item \strong{Non-Parametric
  Bootstrap} This will shuffle the allelic state for each
  individual. Missing data is fixed in place. }
}
\examples{
# load the microbov dataset
data(microbov)
# Let's look at a single population for now. Howsabout Zebu
Zebu <- popsub(microbov, "Zebu")
summary(Zebu)

# Take note of the Number of alleles per population and the Observed
# heterozygosity as we go through each method.

# Multilocus Style: maintain allelic state and heterozygosity.
summary(shufflepop(Zebu, method=1))

# Permute Alleles: maintain allelic state; heterozygosity varies.
summary(shufflepop(Zebu, method=2))

# Parametric Bootstrap: do not maintain allelic state or heterozygosity
summary(shufflepop(Zebu, method=3))

# Non-Parametric Bootstrap: do not maintain allelic state or heterozygosity.
summary(shufflepop(Zebu, method=4))
}
\author{
  Zhian N. Kamvar
}

