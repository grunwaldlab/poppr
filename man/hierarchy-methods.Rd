\docType{methods}
\name{gethierarchy}
\alias{gethierarchy}
\alias{gethierarchy,genclone-method}
\alias{sethierarchy}
\alias{sethierarchy,genclone-method}
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

sethierarchy(x, df) <- value

\S4method{sethierarchy}{genclone}(x, df) <- value

sethierarchy(x, df)

\S4method{sethierarchy}{genclone}(x, df)

setpop(x, formula = NULL)

\S4method{setpop}{genclone}(x, formula = NULL)

setpop(x, formula = NULL) <- value

\S4method{setpop}{genclone}(x, formula) <- value
}
\arguments{
  \item{x}{a genclone object}

  \item{formula}{a nested formula indicating the order of
  the population hierarchy.}

  \item{combine}{if \code{TRUE}, the levels will be
  combined according to the formula argument. If it is
  \code{FALSE}, the levels will not be combined.}

  \item{df}{a data frame giving the population hierarchy of
  each individual.}
}
\description{
word
}

