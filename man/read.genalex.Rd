\name{read.genalex}
\alias{read.genalex}
\title{Importing data from genalex formatted *.csv files.}
\usage{
  read.genalex(genalex, ploidy = 2, geo = FALSE,
    region = FALSE)
}
\arguments{
  \item{genalex}{a *.csv file exported from genalex}

  \item{ploidy}{indicate the ploidy of the dataset}

  \item{geo}{indicates the presence of geographic data in
  the file.}

  \item{region}{indicates the presence of regional data in
  the file.}
}
\value{
  A \code{\link{genind}} object.
}
\description{
  read.genalex will read in a genalex-formatted file that
  has been exported in a comma separated format and will
  parse most types of genalex data. The output is a
  \code{\link{genind}} object.
}
\note{
  this function cannot handle raw allele frequency data.

  This function cannot handle raw allele frequency data. In
  the case that there are duplicated names within the file,
  this function will assume separate individuals and rename
  each one to a sequence of integers from 1 to the number
  of individuals. A vector of the original names will be
  saved in the \code{other} slot under
  \code{original.names}. A data frame called
  \code{population_hierarchy} will also be present in the
  \code{other} slot and it is the place you should insert
  population structure for clone correction analysis.
}
\examples{
Aeut <- read.genalex(system.file("files/rootrot.csv", package="poppr"))

\dontrun{
genalex2 <- read.genalex("genalex2.csv", geo=TRUE)
# A genalex file with geographic coordinate data.

genalex3 <- read.genalex("genalex3.csv", region=TRUE)
# A genalex file with regional information.
}
}
\seealso{
  \code{\link{clonecorrect}}, \code{\link{genind}}
}

