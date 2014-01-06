\name{ia}
\alias{ia}
\title{Index of Association}
\usage{
ia(pop, sample = 0, method = 1, quiet = FALSE, missing = "ignore",
  hist = TRUE)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object OR any fstat,
  structure, genetix, genpop, or genalex formatted files.}

  \item{sample}{an integer indicating the number of
  permutations desired (eg 999).}

  \item{method}{an integer from 1 to 4 indicating the
  sampling method desired. see \code{\link{shufflepop}} for
  details.}

  \item{quiet}{Should the function print anything to the
  screen while it is performing calculations?

  \code{TRUE} prints nothing.

  \code{FALSE} (defualt) will print the population name and
  progress bar.}

  \item{missing}{a character string. see
  \code{\link{missingno}} for details.}

  \item{hist}{\code{logical} if \code{TRUE}, a histogram
  will be printed for each population if there is
  sampling.}
}
\value{
\emph{If no sampling has occured:}

A named number vector of length 2 giving the Index of
Association, "Ia"; and the Standardized Index of
Association, "rbarD"

\emph{If there is sampling:}

A a named number vector of length 4 with the following
values: \item{Ia}{numeric. The index of association.}
\item{p.Ia}{A number indicating the p-value resulting from
a one-sided permutation test based on the number of samples
indicated in the original call.} \item{rbarD}{numeric. The
standardized index of association.} \item{p.rD}{A factor
indicating the p-value resutling from a one-sided
permutation test based on the number of samples indicated
in the original call.}
}
\description{
Calculate the Index of Association and Standardized Index
of Association. Obtain p-values from one-sided permutation
tests.
}
\examples{
data(nancycats)
ia(nancycats)

\dontrun{
# Get the index for each population.
lapply(seppop(nancycats), ia)
# With sampling
lapply(seppop(nancycats), ia, sample=999)
}
}
\author{
Zhian N. Kamvar
}
\references{
Paul-Michael Agapow and Austin Burt. Indices of multilocus
linkage disequilibrium. \emph{Molecular Ecology Notes},
1(1-2):101-102, 2001

A.H.D. Brown, M.W. Feldman, and E. Nevo. Multilocus
structure of natural populations of hordeum spontaneum.
\emph{Genetics}, 96(2):523-536, 1980.

J M Smith, N H Smith, M O'Rourke, and B G Spratt. How
clonal are bacteria? Proceedings of the National Academy of
Sciences, 90(10):4384-4388, 1993.
}
\seealso{
\code{\link{poppr}}, \code{\link{missingno}},
\code{\link{import2genind}}, \code{\link{read.genalex}},
\code{\link{clonecorrect}}
}

