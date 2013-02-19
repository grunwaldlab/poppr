\name{clonecorrect}
\alias{clonecorrect}
\title{Remove potential bias caused by cloned genotypes in genind object.}
\usage{
  clonecorrect(pop, hier = c(1),
    dfname = "population_hierarchy", combine = FALSE)
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
}
\value{
  a clone corrected \code{\link{genind}} object.
}
\description{
  This function removes any duplicated multi locus
  genotypes from any specified population hierarchy.
}
\note{
  This function will clone correct to the population level
  indicated in the \code{pop} slot of the
  \code{\link{genind}} object if there is no data frame
  specified in dfname. If there is no population structure
  and there is no specified data frame, it will clone
  correct the entire \code{\link{genind}} object.
}
\examples{
# LOAD H3N2 dataset
data(H3N2)

# Extract only the individuals located in China
country <- clonecorrect(H3N2, hier=c("country"), dfname="x")

# How many isolates did we have from China before clone correction?
length(which(other(H3N2)$x$country=="China")) # 155

# How many unique isolates from China after clone correction?
length(which(other(country)$x$country=="China")) # 79

\dontrun{
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

