\name{any.boot}
\alias{any.boot}
\title{Calculate a dendrogram with bootstrap support using any distance applicable
to genind or genclone objects.}
\usage{
\method{any}{boot}(x, tree = "upgma", distance = "nei.dist", sample = 100,
  cutoff = 0, showtree = TRUE, ...)
}
\arguments{
  \item{x}{a \linkS4class{genind}, \linkS4class{genclone},
  or matrix object.}

  \item{tree}{one of "upgma" (Default) or "nj" defining the
  type of dendrogram to be produced, UPGMA or
  Neighbor-Joining.}

  \item{distance}{a character or function defining the
  distance to be applied to x. Defaults to
  \code{\link{nei.dist}}.}

  \item{sample}{An integer representing the number of
  bootstrap replicates Default is 100.}

  \item{cutoff}{An integer from 0 to 100 setting the cutoff
  value to return the bootstrap values on the nodes.
  Default is 0.}

  \item{showtree}{If \code{TRUE} (Default), a dendrogram
  will be plotted. If \code{FALSE}, nothing will be
  plotted.}

  \item{...}{any parameters to be passed off to the
  distance method.}
}
\value{
an object of class \code{\link[ape]{phylo}}.
}
\description{
Calculate a dendrogram with bootstrap support using any
distance applicable to genind or genclone objects.
}
\details{
This function utilizes an internal class called
\code{\linkS4class{bootgen}} that allows bootstrapping of
objects that inherit the genind class. This is necessary
due to the fact that columns in the genind matrix are
defined as alleles and are thus interrelated. This function
will specifically bootstrap loci so that results are
biologically relevant. With this function, the user can
also define a custom distance to be performed on the genind
or genclone object.
}
\examples{
data(nancycats)
nan9 <- popsub(nancycats, 9)

set.seed(9999)
# Generate a tree using nei's distance
neinan <- any.boot(nan9, dist = nei.dist)

set.seed(9999)
# Generate a tree using custom distance
bindist <- function(x) dist(x, method = "binary")
binnan <- any.boot(nan9, dist = bindist)

\dontrun{
# AFLP data
data(Aeut)

# Nei's distance
anei <- any.boot(Aeut, dist = nei.dist, sample = 1000, cutoff = 50)

# Rodger's distance
arog <- any.boot(Aeut, dist = rodger.dist, sample = 1000, cutoff = 50)

}
}
\seealso{
\code{\link{nei.dist}} \code{\link{edward.dist}}
\code{\link{rodger.dist}} \code{\link{reynold.dist}}
\code{\link{provesti.dist}} \code{\link{diss.dist}}
\code{\link{bruvo.boot}} \code{\link[ape]{boot.phylo}}
\code{\link[adegenet]{dist.genpop}} \code{\link{dist}}
}

