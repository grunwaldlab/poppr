\name{genind2genalex}
\alias{genind2genalex}
\title{Exporting data from genind objects to genalex formatted *.csv files.}
\usage{
  genind2genalex(pop, filename = "genalex.csv",
    quiet = FALSE)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object.}

  \item{filename}{a string indicating the name and/or path
  of the file you wish to create.}

  \item{quiet}{\code{logical} If \code{FALSE} a message
  will be printed to the screen.}
}
\description{
  genind2genalex will export a genind object to a *.csv
  file formatted for use in genalex.
}
\note{
  If you enter a file name that exists, that file will be
  overwritten. If your data set lacks a population
  structure, it will be coded in the new file as a single
  population lableled "Pop". Likewise, if you don't have
  any labels for your individuals, they will be labeled as
  "ind1" through "ind\emph{N}", with \emph{N} being the
  size of your population.
}
\examples{
\dontrun{
data(nancycats)
genind2genalex(nancycats, "~/Documents/nancycats.csv")
}
}
\seealso{
  \code{\link{clonecorrect}}, \code{\link{genind}}
}

