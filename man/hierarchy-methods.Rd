\docType{methods}
\name{gethierarchy}
\alias{gethierarchy}
\alias{gethierarchy,genclone-method}
\alias{namehierarchy}
\alias{namehierarchy,genclone-method}
\alias{namehierarchy<-}
\alias{namehierarchy<-,genclone-method}
\alias{sethierarchy}
\alias{sethierarchy,genclone-method}
\alias{sethierarchy<-}
\alias{sethierarchy<-,genclone-method}
\title{Access and manipulate the population hierarchy for genclone objects.}
\usage{
gethierarchy(x, formula = NULL, combine = TRUE)

sethierarchy(x, value)

sethierarchy(x) <- value

namehierarchy(x, value)

namehierarchy(x) <- value
}
\arguments{
  \item{x}{a genclone object}

  \item{formula}{a nested formula indicating the order of
  the population hierarchy.}

  \item{combine}{if \code{TRUE}, the levels will be
  combined according to the formula argument. If it is
  \code{FALSE}, the levels will not be combined.}

  \item{value}{a data frame OR vector OR formula (see
  details).}
}
\description{
The following methods allow the user to quickly change the
hierarchy or population of a genclone object.
}
\details{
These functions allow the user to seamlessly assign the
hierarchy of their \code{\linkS4class{genclone}} object.
Note that there are two ways of performing
\code{sethierarchy()} and \code{namehierarchy()}, they
essentially do the same thing except that the assignment
method (the one with the "\code{<-}") will modify the
object in place whereas the non-assignment method will not
modify the original object. Due to convention, everything
right of the assignment is termed \code{value}. To avoid
confusion, here is a guide: \enumerate{ \item
\strong{sethierarchy()} This will be a
\code{\link{data.frame}} that defines the hierarchy for
each individual in the rows.  \item
\strong{namehierarchy()} This will be either a
\code{\link{vector}} or a \code{\link{formula}} that will
define the names. }
}
\examples{
data(Aeut)
Aeut.gc <- as.genclone(Aeut)
head(gethierarchy(Aeut.gc))
popsub <- other(Aeut)$population_hierarchy[-1]
Aeut.gc <- sethierarchy(Aeut.gc, popsub)
sethierarchy(Aeut.gc) <- popsub
head(gethierarchy(Aeut.gc))
namehierarchy(Aeut.gc) <- ~Population/Subpopulation
Aeut.gc
Aeut.old <- namehierarchy(Aeut.gc, ~Pop/Subpop)
Aeut.old
}

