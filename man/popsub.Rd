% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_subset.r
\name{popsub}
\alias{popsub}
\title{Subset data by population}
\usage{
popsub(
  gid,
  sublist = "ALL",
  exclude = NULL,
  blacklist = NULL,
  mat = NULL,
  drop = TRUE
)
}
\arguments{
\item{gid}{a \code{\link[adegenet:genind-class]{genind}}, \code{\linkS4class{genclone}},
\code{\link[adegenet:genlight-class]{genlight}}, or \code{\linkS4class{snpclone}} object.}

\item{sublist}{a \code{vector} of population names or indexes that the user
wishes to keep. Default to \code{"ALL"}.}

\item{exclude}{a \code{vector} of population names or indexes that the user
wishes to discard. Default to \code{NULL}.}

\item{blacklist}{DEPRECATED, use exclude.}

\item{mat}{a \code{matrix} object produced by \code{\link{mlg.table}} to be
subsetted. If this is present, the subsetted matrix will be returned instead
of the genind object}

\item{drop}{\code{logical}. If \code{TRUE}, unvarying alleles will be dropped
from the population.}
}
\value{
A \code{genind} object or a matrix.
}
\description{
Create a new dataset with specified populations or exclude specified
populations from the dataset.
}
\examples{
# Load the dataset microbov.
data(microbov)

# List the population names.
popNames(microbov)

# Analyze only the populations with exactly 50 individuals
mic.50 <- popsub(microbov, sublist=c(1:6, 11:15), exclude=c(3,4,13,14))

\dontrun{
# Analyze the first 10 populations, except for "Bazadais"
mic.10 <- popsub(microbov, sublist=1:10, exclude="Bazadais")

# Take out the two smallest populations
micbig <- popsub(microbov, exclude=c("NDama", "Montbeliard"))

# Analyze the two largest populations
miclrg <- popsub(microbov, sublist=c("BlondeAquitaine", "Charolais"))
}
}
\author{
Zhian N. Kamvar
}
