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
  This function cannot handle raw allele frequency data.

  The resulting genind object will have a data frame in the
  \code{other} slot called population_hierarchy. This will
  contain a column for your population data and a column
  for your Regional data if you have set the flag.

  If there is geographic data, it will be included in a
  data frame called xy in the \code{other} slot.

  In the case that there are duplicated names within the
  file, this function will assume separate individuals and
  rename each one to a sequence of integers from 1 to the
  number of individuals. A vector of the original names
  will be saved in the \code{other} slot under
  \code{original_names}.
}
\examples{
Aeut <- read.genalex(system.file("files/rootrot.csv", package="poppr"))

\dontrun{
genalex2 <- read.genalex("genalex2.csv", geo=TRUE)
# A genalex file with geographic coordinate data.

genalex3 <- read.genalex("genalex3.csv", region=TRUE)
# A genalex file with regional information.

genalex4 <- read.genalex("genalex4.csv", region=TRUE, geo=TRUE)
# A genalex file with both regional and geographic information.
}
}
\author{
  Zhian N. Kamvar
}
\seealso{
  \code{\link{clonecorrect}}, \code{\link{genind}}
}

