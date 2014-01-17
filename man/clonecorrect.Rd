\name{clonecorrect}
\alias{clonecorrect}
\title{Remove potential bias caused by cloned genotypes in genind object.}
\usage{
clonecorrect(pop, hier = 1, dfname = "population_hierarchy",
  combine = FALSE, keep = 1)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object}

  \item{hier}{a \code{numeric or character list}. This is
  the list of vectors within a data frame (specified in
  \code{dfname}) in the 'other' slot of the
  \code{\link{genind}} object. The list should indicate the
  population hierarchy to be used for clone correction.}

  \item{dfname}{a \code{character string}. This is the name
  of the data frame or list containing the vectors of the
  population hierarchy within the \code{other} slot of the
  \code{\link{genind}} object.}

  \item{combine}{\code{logical}. When set to TRUE, the
  heirarchy will be combined to create a new population for
  the genind object.}

  \item{keep}{\code{integer}. When \code{combine} is set to
  \code{FALSE}, you can use this flag to choose the levels
  of your population hierarchy. For example: if your clone
  correction hierarchy is set to "Pop", "Subpop", and
  "Year", and you want to analyze your populations with
  respect to year, you can set \code{keep = c(1,3)}.}
}
\value{
a clone corrected \code{\link{genind}} object.
}
\description{
This function removes any duplicated multi locus genotypes
from any specified population hierarchy.
}
\note{
This function will clone correct to the population level
indicated in the \code{pop} slot of the
\code{\link{genind}} object if there is no data frame
specified in dfname. If there is no population structure
and there is no specified data frame, it will clone correct
the entire \code{\link{genind}} object.
}
\examples{
# LOAD A. euteiches data set
data(Aeut)

# Check the number of multilocus genotypes
mlg(Aeut)
Aeut$pop.names

# Clone correct at the population level.
Aeut.pop <- clonecorrect(Aeut, hier="Pop")
mlg(Aeut.pop)
Aeut.pop$pop.names

\dontrun{
# Clone correct at the subpopulation level with respect to population and
# combine.
Aeut.subpop <- clonecorrect(Aeut, hier=c("Pop", "Subpop"), combine=TRUE)
mlg(Aeut.subpop)
Aeut.subpop$pop.names

# Do the same, but set to the population level.
Aeut.subpop2 <- clonecorrect(Aeut, hier=c("Pop", "Subpop"), keep=1)
mlg(Aeut.subpop2)
Aeut.subpop2$pop.names

# LOAD H3N2 dataset
data(H3N2)

# Extract only the individuals located in China
country <- clonecorrect(H3N2, hier=c("country"), dfname="x")

# How many isolates did we have from China before clone correction?
length(which(other(H3N2)$x$country=="China")) # 155

# How many unique isolates from China after clone correction?
length(which(other(country)$x$country=="China")) # 79

# Something a little more complicated. (This could take a few minutes on
# slower computers)

# setting the hierarchy to be Country > Year > Month
c.y.m <- clonecorrect(H3N2, hier=c("year","month","country"), dfname="x")

# How many isolates in the original data set?
length(other(H3N2)$x$country) # 1903

# How many after we clone corrected for country, year, and month?
length(other(c.y.m)$x$country) # 1190
}
}
\author{
Zhian N. Kamvar
}

