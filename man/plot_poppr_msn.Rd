\name{plot_poppr_msn}
\alias{plot_poppr_msn}
\title{Plot minimum spanning networks produced in poppr.}
\usage{
plot_poppr_msn(x, poppr_msn, gadj = 3, glim = c(0, 0.8), gweight = 1,
  inds = "ALL", quantiles = TRUE, nodelab = 2, cutoff = NULL,
  palette = NULL, layfun = layout.auto, beforecut = FALSE, ...)
}
\arguments{
  \item{x}{a \code{\linkS4class{genind}} or
  \code{\linkS4class{genclone}} object from which
  \code{poppr_msn} was derived.}

  \item{poppr_msn}{a \code{list} produced from either
  \code{\link{poppr.msn}} or \code{\link{bruvo.msn}}. This
  list should contain a graph, a vector of population names
  and a vector of hexadecimal color definitions for each
  popualtion.}

  \item{inds}{a character vector indicating which
  individual names to label nodes with. See details.}

  \item{quantiles}{\code{logical}. When set to \code{TRUE}
  (default), the scale bar will be composed of the
  quantiles from the observed edge weights. When set to
  \code{FALSE}, the scale bar will be composed of a smooth
  gradient from the minimum edge weight to the maximum edge
  weight.}

  \item{nodelab}{an \code{integer} specifying the smallest
  size of node to label. See details.}

  \item{cutoff}{a number indicating the longest distance to
  display in your graph. This is performed by removing
  edges with weights greater than this number.}

  \item{palette}{a function or character corresponding to a
  specific palette you want to use to delimit your
  populations. The default is whatever palette was used to
  produce the original graph.}

  \item{layfun}{a function specifying the layout of nodes
  in your graph. It defaults to
  \code{\link[igraph]{layout.auto}}.}

  \item{beforecut}{if \code{TRUE}, the layout of the graph
  will be computed before any edges are removed with
  \code{cutoff}. If \code{FALSE} (Default), the layout will
  be computed after any edges are removed.}

  \item{...}{any other parameters to be passed on to
  \code{\link[igraph]{plot.igraph}}.}

  \item{gadj}{"grey adjust". a positive \code{integer}
  greater than zero that will serve as the exponent to the
  edge weight to scale the grey value to represent that
  weight.}

  \item{glim}{"grey limit". Two numbers between zero and
  one. They determine the upper and lower limits for the
  \code{\link{gray}} function. Default is 0 (black) and 0.8
  (20\% black).}

  \item{gweight}{"grey weight". an \code{integer}. If it's
  1, the grey scale will be weighted to emphasize the
  differences between closely related nodes. If it is 2,
  the grey scale will be weighted to emphasize the
  differences between more distantly related nodes.}
}
\description{
This function allows you to take the output of poppr.msn
and bruvo.msn and customize the plot by labeling groups of
individuals, size of nodes, and adjusting the palette and
scale bar.
}
\details{
The previous incarnation of msn plotting in poppr simply
plotted the minimum spanning network with the legend of
populations, but did not provide a scale bar and it did not
provide the user a simple way of manipulating the layout or
labels. This function allows the user to manipulate many
facets of graph creation, making the creation of minimum
spanning networks ever so slightly more user friendly. Note
that this function will only plot individual names, not MLG
names since the naming convention for those are arbitrary.
This will also not give the user the option to remove the
shading or widths of the edges as they are informative to
the distance.

This function must have both the source data and the output
msn to work.  The source data must also have the same
population structure as the graph.  Every other parameter
has a default setting.

Parameter details: \itemize{ \item \code{inds} This will
take in the name of a query individual in your data set and
will use that to query any other individuals that share
multilocus genotypes and label their node on the graph. The
default is to label all the nodes, but you can set it to a
name that doesn't exist to label none of the nodes. \item
\code{nodelab} If a node is not labeled by individual, this
will label the size of the nodes greater than or equal to
this value. If you don't want to label the size of the
nodes, simply set this to a very high number. \item
\code{cutoff} This is useful for when you want to
investigate groups of multilocus genotypes separated by a
specific distance or if you have two distinct populations
and you want to physically separate them in your network.
\item \code{beforecut} This is an indicator useful if you
want to maintain the same position of the nodes before and
after removing edges with the \code{cutoff} argument. This
works best if you set a seed before you run the function.}
}
\examples{
# Using a data set of the Aphanomyces eutieches root rot pathogen.
data(Aeut)
amsn <- poppr.msn(Aeut, diss.dist(Aeut), showplot = FALSE)

# Default
library(igraph) # To get all the layouts.
set.seed(500)
plot_poppr_msn(Aeut, amsn, gadj = 15, beforecut = TRUE)

# Removing link between populations and labelling no individuals
set.seed(500)
plot_poppr_msn(Aeut, amsn, inds = "none", gadj = 15, beforecut = TRUE, cutoff = 0.2)

# Labelling individual #57 because it is an MLG that crosses popualtions
# Showing clusters of MLGS with at most 5\% variation
# Notice that the Mt. Vernon population appears to be more clonal
set.seed(50)
plot_poppr_msn(Aeut, amsn, gadj = 15, cutoff = 0.05, inds = "57")


\dontrun{
data(partial_clone)
pcmsn <- bruvo.msn(partial_clone, replen = rep(1, 10))
plot_poppr_msn(partial_clone, pcmsn, palette = rainbow, inds = "sim 20")

# Something pretty
data(microbov)
micmsn <- poppr.msn(microbov, diss.dist(microbov), showplot = FALSE)

plot_poppr_msn(microbov, micmsn, palette = "terrain.colors", quantiles = FALSE)
plot_poppr_msn(microbov, micmsn, palette = "terrain.colors", cutoff = 0.3, quantiles = FALSE)
}
}
\author{
Zhian N. Kamvar
}
\seealso{
\code{\link[igraph]{layout.auto}}
\code{\link[igraph]{plot.igraph}} \code{\link{poppr.msn}}
\code{\link{bruvo.msn}} \code{\link{greycurve}}
\code{\link[igraph]{delete.edges}} \code{\link{palette}}
}

