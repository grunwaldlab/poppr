\name{poppr.amova}
\alias{poppr.amova}
\title{Perform Analysis of Molecular Variance (AMOVA) on genind or genclone objects.}
\usage{
poppr.amova(x, hier = NULL, clonecorrect = FALSE, within = TRUE,
  dist = NULL, squared = TRUE, correction = "quasieuclid",
  dfname = "population_hierarchy", sep = "_", missing = "loci",
  cutoff = 0.05, quiet = FALSE)
}
\arguments{
  \item{x}{a \code{\linkS4class{genind}} or
  \code{\linkS4class{genclone}} object}

  \item{hier}{a hierarchical \code{\link{formula}} that
  defines your population hieararchy. (e.g.:
  ~Population/Subpopulation). \strong{See Details below.}}

  \item{clonecorrect}{\code{logical} if \code{TRUE}, the
  data set will be clone corrected with respect to the
  lowest level of the hierarchy. The default is set to
  \code{FALSE}. See \code{\link{clonecorrect}} for
  details.}

  \item{within}{\code{logical}. When this is set to
  \code{TRUE} (Default), variance within individuals are
  calculated as well. If this is set to \code{FALSE}, The
  lowest level of the hierarchy will be the sample level.
  See Details below.}

  \item{dist}{an optional distance matrix calculated on
  your data.}

  \item{squared}{if a distance matrix is supplied, this
  indicates whether or not it represents squared
  distances.}

  \item{correction}{a \code{character} defining the
  correction method for non-euclidean distances. Options
  are \code{\link[ade4]{quasieuclid}} (Default),
  \code{\link[ade4]{lingoes}}, and
  \code{\link[ade4]{cailliez}}. See Details below.}

  \item{dfname}{if the input data set is a
  \code{\linkS4class{genind}} object, specify the name of
  the data frame in the \code{\link[adegenet]{other}} slot
  defining the population hierarchy. Defaults to
  \code{"population_hierarchy"}}

  \item{sep}{What character is used to separate the
  popualtion hierarchy?}

  \item{missing}{specify method of correcting for missing
  data utilizing options given in the function
  \code{\link{missingno}}. Default is \code{"loci"}.}

  \item{cutoff}{specify the level at which missing data
  should be removed/modified. See \code{\link{missingno}}
  for details.}

  \item{quiet}{\code{logical} If \code{FALSE} (Default),
  messages regarding any corrections will be printed to the
  screen. If \code{TRUE}, no messages will be printed.}
}
\value{
a list of class \code{amova} from the ade4 package. See
\code{\link[ade4]{amova}} for details.
}
\description{
This function utilizes the ade4 implementation of AMOVA.
See \code{\link[ade4]{amova}} for details on the specific
implementation.
}
\examples{
data(Aeut)
agc <- as.genclone(Aeut)
agc
amova.result <- poppr.amova(agc, ~Pop/Subpop)
amova.result[1:4]
\dontrun{
amova.cc.result <- poppr.amova(agc, ~Pop/Subpop, clonecorrect = TRUE)
amova.cc.result[1:4]
}
}

