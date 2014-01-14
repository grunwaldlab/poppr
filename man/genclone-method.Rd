\docType{methods}
\name{[,genclone,ANY,ANY,ANY-method}
\alias{[,genclone,ANY,ANY,ANY-method}
\alias{initialize,genclone-method}
\alias{show,genclone-method}
\title{Methods used for the genclone object}
\usage{
\S4method{[}{genclone,ANY,ANY,ANY}(x, i, j, ..., drop = FALSE)

\S4method{initialize}{genclone}(.Object, gen, hierarchy, mlg)

\S4method{show}{genclone}(object)
}
\arguments{
  \item{x}{a genclone object}

  \item{i}{vector of numerics indicating number of
  individuals desired}

  \item{j}{a vector of numerics corresponding to the loci
  desired.}

  \item{...}{passed on to the \code{\linkS4class{genind}}
  object.}

  \item{drop}{set to \code{FALSE}}

  \item{.Object}{a character, "genclone"}

  \item{gen}{\code{"\linkS4class{genind}"} object}

  \item{hierarchy}{a data frame where each row i represents
  the different population assignments of individual i in
  the data set. If this is empty, the hierarchy will be
  created from the population factor.}

  \item{mlg}{a vector where each element assigns the
  multilocus genotype of that individual in the data set.}

  \item{object}{a genclone object}
}
\description{
words.
}

