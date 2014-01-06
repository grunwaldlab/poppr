\name{splitcombine}
\alias{splitcombine}
\title{Split a or combine items within a data frame in \code{\link{genind}} objects.}
\usage{
  splitcombine(pop, method = 1,
    dfname = "population_hierarchy", sep = "_",
    hier = c(1), setpopulation = TRUE, fixed = TRUE)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object.}

  \item{method}{an \code{integer}, 1 for splitting, 2 for
  combining.}

  \item{dfname}{the name of the data frame containing the
  population structure. for the splitting method, the
  combined population structure must be in the first
  column.}

  \item{sep}{The separator used for separating or combining
  the data. See note.}

  \item{hier}{a \code{vector} containing the population
  hierarchy you wish to split or combine.}

  \item{setpopulation}{\code{logical}. if \code{TRUE}, the
  population of the resulting genind object will be that of
  the highest population structure (split method) or the
  combined populations (combine method).}

  \item{fixed}{\code{logical}. An argument to be passed
  onto \code{\link{strsplit}}. If \code{TRUE}, \code{sep}
  must match exactly to the populations for the split
  method.}
}
\value{
  a \code{\link{genind}} object with a modified data frame
  in the \code{\link{other}} slot.
}
\description{
  Often, one way a lot of file formats fail is that they do
  not allow multiple population hierarchies. This can be
  circumvented, however, by coding all of the hierarchies
  in one string in the input file with a common separator
  (eg. "_"). \code{splitcombine} will be able to recognise
  those separators and create a data frame of all the
  population structures for whatever subsetting you might
  need.
}
\note{
  The separator field is sensitive to regular expressions.
  If you do not know what those are, please use the default
  underscore to separate your populations. Use \code{fixed
  = TRUE} to ignore regular expressions. If you do not set
  the \code{hier} flag for the split method, your new data
  frame will have the names "comb", "h1", "h2" and so on;
  for the combine method, your data frame will return the
  first column of your data frame.
}
\examples{
\dontrun{
# Method 1: Splitting.
Aeut <- read.genalex(system.file("files/rootrot.csv", package="poppr"))

# We have 19 different "populations", but really, there is a hierarchy.
Aeut$pop.names

# Let's split them up. The default data frame from read.genalex is the same
# as the default for this function.
Aeut <- splitcombine(Aeut, hier=c("Pop", "Subpop"))

# Much better!
Aeut$pop.names
}

# Method 2: Combining.

data(H3N2)
# Create a new data set combining the population factors of year and country
H.comb <- splitcombine(H3N2, method=2, dfname="x", hier=c("year", "country"))

# Checking to make sure they were actually combined.
head(H.comb$other$x$year_country)
\dontrun{
# Creating new data frame in the object to mess around with.
H.comb$other$year_country <- data.frame(H.comb$other$x$year_country)

# Splitting those factors into their original components and setting the
# population to year.
H.comb <- splitcombine(H.comb, method=1, dfname="year_country", hier=c("year", "country"))
}
}
\author{
  Zhian N. Kamvar
}

