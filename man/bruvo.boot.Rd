\name{bruvo.boot}
\alias{bruvo.boot}
\title{Create a tree using Bruvo's Distance with non-parametric bootstrapping.}
\usage{
  bruvo.boot(pop, replen = 1, add = TRUE, loss = TRUE,
    sample = 100, tree = "upgma", showtree = TRUE,
    cutoff = NULL, quiet = FALSE, ...)
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

  \item{sample}{an \code{integer} indicated the number of
  bootstrap replicates desired.}

  \item{tree}{choose between "nj" for neighbor-joining and
  "upgma" for a upgma tree to be produced.}

  \item{showtree}{\code{logical} if \code{TRUE}, a tree
  will be plotted with nodelabels.}

  \item{cutoff}{\code{integer} the cutoff value for
  bootstrap node label values (between 0 and 100).}

  \item{quiet}{\code{logical} defaults to \code{FALSE}. If
  \code{TRUE}, a progress bar and messages will be
  supressed.}

  \item{...}{any argument to be passed on to
  \code{\link{boot.phylo}}. eg. \code{quiet = TRUE}.}
}
\value{
  a tree with nodelables
}
\description{
  Create a tree using Bruvo's Distance with non-parametric
  bootstrapping.
}
\note{
  This function calculates bruvo's distance for non-special
  cases (ie. the ploidy and all alleles are known).
  Currently there is no way to import polyploid partial
  heterozygote data into adegenet. For Bruvo's Distance
  concerning special cases, see the package \code{polysat}.
  Missing data is ignored, but be sure that missing data is
  NOT set to 0 in the genind object. This is not easy to
  detect and will result in an error. Please use any other
  method in \code{\link{na.replace}} or
  \code{\link{missingno}}.

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

bruvo.boot(popsub(nancycats, 1), replen=ssr)
}
\author{
  Javier F. Tabima, Zhian N. Kamvar
}
\references{
  Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza,
  and Hinrich Schulenburg. A simple method for the
  calculation of microsatellite genotype distances
  irrespective of ploidy level. Molecular Ecology,
  13(7):2101-2106, 2004.
}
\seealso{
  \code{\link{nancycats}}, \code{\link{upgma}},
  \code{\link{nj}}, \code{\link{boot.phylo}},
  \code{\link{nodelabels}}, \code{\link{na.replace}},
  \code{\link{missingno}}.
}

