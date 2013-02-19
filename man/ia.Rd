\name{ia}
\alias{ia}
\title{Index of Association}
\usage{
  ia(pop, sample = 0, method = 1, quiet = "minimal",
    missing = "ignore", hist = TRUE)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object OR any fstat,
  structure, gentix, genpop, or genalex formatted files.}

  \item{sample}{an integer indicating the number of
  permutations desired (eg 999).}

  \item{method}{an integer from 1 to 4 indicating the
  sampling method desired. see \code{popsample} for
  details.}

  \item{quiet}{Should the function print anything to the
  screen while it is performing calculations?

  \code{TRUE} prints nothing.

  \code{FALSE} same as "minimal".

  \code{"minimal"} (defualt) will print the population name
  and dots indicating permutation progress.

  \code{"noisy"} will print out the individual indecies as
  they are produced.}

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
  \item{p.Ia}{A number indicating the p-value resulting
  from a one-sided permutation test based on the number of
  samples indicated in the original call.}
  \item{rbarD}{numeric. The standardized index of
  association.} \item{p.rD}{A factor indicating the p-value
  resutling from a one-sided permutation test based on the
  number of samples indicated in the original call.}
}
\description{
  Calculate the Index of Association and Standardized Index
  of Association. Obtain p-values from one-sided
  permutation tests.
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
\seealso{
  \code{\link{poppr}}, \code{\link{missingno}},
  \code{\link{import2genind}}, \code{\link{read.genalex}},
  \code{\link{clonecorrect}}
}

