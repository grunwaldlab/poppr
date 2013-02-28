\name{bruvo.msn}
\alias{bruvo.msn}
\title{Create minimum spanning network of selected populations using Brvuo's
distance.}
\usage{
  bruvo.msn(pop, replen = c(1), palette = topo.colors,
    sublist = "All", blacklist = NULL,
    vertex.label = "MLG", gscale = TRUE, glim = c(0, 0.8),
    wscale = TRUE, ...)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object}

  \item{replen}{a \code{vector} of \code{integers}
  indicating the length of the nucleotide repeats for each
  microsatellite locus.}

  \item{palette}{a \code{function} defining the color
  palette to be used to color the populations on the graph.
  It defaults to \code{\link{topo.colors}}, but you can
  easily create new schemes by using
  \code{\link{colorRampPalette}} (see examples for
  details)}

  \item{sublist}{a \code{vector} of population names or
  indexes that the user wishes to keep. Default to "ALL".}

  \item{blacklist}{a \code{vector} of population names or
  indexes that the user wishes to discard. Default to
  \code{NULL}}

  \item{vertex.label}{a \code{vector} of characters to
  label each vertex. There are two defaults: \code{"MLG"}
  will label the nodes with the multilocus genotype from
  the original data set and \code{"inds"} will label the
  nodes with the representative individual names.}

  \item{gscale}{"grey scale". If this is \code{TRUE}, this
  will scale the color of the edges proportional to Bruvo's
  distance, with the lines becoming darker for more related
  nodes.}

  \item{glim}{"grey limit". Two numbers between zero and
  one. They determine the upper and lower limits for the
  \code{\link{gray}} function. Default is 0 (black) and 0.8
  (20\% black).}

  \item{wscale}{"width scale". If this is \code{TRUE}, the
  edge widths will be scaled proportional to Bruvo's
  distance, with the lines becoming thicker for more
  related nodes.}

  \item{...}{any other arguments that could go into
  \code{\link{igraph:::plot.igraph}}}
}
\value{
  a minimum spanning network with nodes corresponding to
  MLGs within the data set. Colors of the nodes represent
  population membership, and length of edges represent
  Bruvo's distance.
}
\description{
  Create minimum spanning network of selected populations
  using Brvuo's distance.
}
\note{
  The edges of these graphs may cross each other if the
  graph becomes too large.
}
\examples{
# Load the data set.
data(nancycats)

# View populations 8 and 9 with default colors.
bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds",
vertex.label.cex=0.7, vertex.label.dist=0.4)
\dontrun{
# View heat colors.
bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds",
palette=heat.colors, vertex.label.cex=0.7, vertex.label.dist=0.4)

# View custom colors. Here, we use black and orange.
bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds",
palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7,
vertex.label.dist=0.4)

# View with darker shades of grey (setting the upper limit to 1/2 black 1/2 white).
bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds",
palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7,
vertex.label.dist=0.4, glim=c(0, 0.5))

# View with no grey scaling.
bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds",
palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7,
vertex.label.dist=0.4, gscale=FALSE)

# View with no line widths.
bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds",
palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7,
vertex.label.dist=0.4, wscale=FALSE)

# View with no scaling at all.
bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds",
palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7,
vertex.label.dist=0.4, vscale=FALSE, gscale=FALSE)

# View the whole population, but without labels.
bruvo.msn(nancycats, replen=rep(1, 9), vertex.label=NA)
}
}
\seealso{
  \code{\link{nancycats}}, \code{\link{upgma}},
  \code{\link{nj}}, \code{\link{boot.phylo}},
  \code{\link{nodelabels}}, \code{\link{na.replace}},
  \code{\link{missingno}}, \code{\link{bruvo.boot}}.
}

