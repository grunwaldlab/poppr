\docType{methods}
\name{[,bootgen,ANY,ANY,ANY-method}
\alias{[,bootgen,ANY,ANY,ANY-method}
\alias{dim,bootgen-method}
\alias{initialize,bootgen-method}
\title{Methods used for the bootgen object.}
\usage{
\S4method{[}{bootgen,ANY,ANY,ANY}(x, i, j, ..., drop = FALSE)

\S4method{dim}{bootgen}(x)

\S4method{initialize}{bootgen}(.Object, gen, replen)
}
\arguments{
  \item{x}{a \code{"\linkS4class{bootgen}"} object}

  \item{i}{vector of numerics indicating number of
  individuals desired}

  \item{j}{a vector of numerics corresponding to the loci
  desired.}

  \item{...}{unused.}

  \item{drop}{set to \code{FALSE}}

  \item{.Object}{a character, "bootgen"}

  \item{gen}{\code{"\linkS4class{genind}"} object}

  \item{replen}{a vector of numbers indicating the repeat
  length for each microsatellite locus.}
}
\description{
This is not designed for user interaction.
}

