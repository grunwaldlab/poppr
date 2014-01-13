\docType{methods}
\name{as.genclone}
\alias{as.genclone}
\alias{as.genclone,as.genind,genclone-method}
\alias{as.genind}
\alias{as.genind,as.genclone,genclone-method}
\title{Create a genclone object from a genind object.}
\usage{
as.genclone(gen, hierarchy = NULL)

as.genind(gcl)
}
\arguments{
  \item{gen}{a \code{\linkS4class{genind}} object}

  \item{hierarchy}{a data frame representing the population
  hierarchy.}

  \item{gcl}{a \code{\linkS4class{genclone}}}
}
\description{
Wrapper for genclone initializer. The hierarchy must have
the same number of rows as the number of observations in
the genind object. If no hierarchy is defined, the
population will be set as the hierarchy under the label
"Pop".
}
\examples{
data(Aeut)
Aeut
Aeut.gc <- as.genclone(Aeut)
Aeut.gc
Aeut.gc <- as.genclone(Aeut, other(Aeut)$population_hierarchy[-1])
Aeut.gc
Aeut <- as.genind(Aeut.gc)
}

