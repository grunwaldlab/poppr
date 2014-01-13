\docType{methods}
\name{gethierarchy}
\alias{gethierarchy}
\alias{gethierarchy,genclone-method}
\alias{sethierarchy}
\alias{sethierarchy,genclone-method}
\alias{sethierarchy<-}
\alias{sethierarchy<-,genclone-method}
\title{Access and manipulate the population hierarchy for genclone objects.}
\usage{
gethierarchy(x, formula = NULL, combine = TRUE)

sethierarchy(x) <- value

sethierarchy(x, value)
}
\arguments{
  \item{x}{a genclone object}

  \item{formula}{a nested formula indicating the order of
  the population hierarchy.}

  \item{combine}{if \code{TRUE}, the levels will be
  combined according to the formula argument. If it is
  \code{FALSE}, the levels will not be combined.}

  \item{value}{a data frame giving the population hierarchy
  of each individual.}
}
\description{
The following methods allow the user to quickly change the
hierarchy or population of a genclone object.
}

