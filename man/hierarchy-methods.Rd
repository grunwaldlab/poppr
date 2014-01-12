\docType{methods}
\name{gethierarchy}
\alias{gethierarchy}
\alias{gethierarchy,genclone-method}
\alias{sethierarchy<-}
\alias{sethierarchy<-,genclone-method}
\alias{setpop}
\alias{setpop,genclone-method}
\alias{setpop<-}
\alias{setpop<-,genclone-method}
\title{Access and organize the population hierarchy for genclone objects.}
\usage{
gethierarchy(x, formula = NULL, combine = TRUE)

\S4method{gethierarchy}{genclone}(x, formula = NULL, combine = TRUE)

sethierarchy(x) <- value

\S4method{sethierarchy}{genclone}(x) <- value

setpop(x, formula = NULL)

\S4method{setpop}{genclone}(x, formula)

setpop(x) <- value

\S4method{setpop}{genclone}(x) <- value
}
\arguments{
  \item{x}{a genclone object}

  \item{formula}{a nested formula indicating the order of
  the population hierarchy. Defaults to \code{NULL}.}

  \item{combine}{if \code{TRUE}, the levels will be
  combined according to the formula argument. If it is
  \code{FALSE}, the levels will not be combined.}

  \item{value}{a data frame to replace the population
  hierarchy with.}

  \item{formula}{a formula defining the levels of the
  population hierarchy to set the population factor with.}

  \item{value}{a formula to replace the population
  hierarchy with.}
}
\description{
word
}

