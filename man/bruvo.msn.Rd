% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bruvo.r
\name{bruvo.msn}
\alias{bruvo.msn}
\alias{msn.bruvo}
\title{Create minimum spanning network of selected populations using Bruvo's 
distance.}
\usage{
bruvo.msn(
  gid,
  replen = 1,
  add = TRUE,
  loss = TRUE,
  mlg.compute = "original",
  palette = topo.colors,
  sublist = "All",
  exclude = NULL,
  blacklist = NULL,
  vertex.label = "MLG",
  gscale = TRUE,
  glim = c(0, 0.8),
  gadj = 3,
  gweight = 1,
  wscale = TRUE,
  showplot = TRUE,
  include.ties = FALSE,
  threshold = NULL,
  clustering.algorithm = NULL,
  ...
)
}
\arguments{
\item{gid}{a \code{\link[adegenet]{genind}} or \code{\link{genclone}} object}

\item{replen}{a \code{vector} of \code{integers} indicating the length of the
nucleotide repeats for each microsatellite locus.}

\item{add}{if \code{TRUE}, genotypes with zero values will be treated under 
the genome addition model presented in Bruvo et al. 2004.}

\item{loss}{if \code{TRUE}, genotypes with zero values will be treated under 
the genome loss model presented in Bruvo et al. 2004.}

\item{mlg.compute}{if the multilocus genotypes are set to "custom" (see 
\code{\link{mll.custom}} for details) in your genclone object, this will 
specify which mlg level to calculate the nodes from. See details.}

\item{palette}{a \code{vector} or \code{function} defining the color palette 
to be used to color the populations on the graph. It defaults to 
\code{\link{topo.colors}}. See examples for details.}

\item{sublist}{a \code{vector} of population names or indexes that the user 
wishes to keep. Default to "ALL".}

\item{exclude}{a \code{vector} of population names or indexes that the user
wishes to discard. Default to \code{NULL}.}

\item{blacklist}{DEPRECATED, use exclude.}

\item{vertex.label}{a \code{vector} of characters to label each vertex. There
are two defaults: \code{"MLG"} will label the nodes with the multilocus 
genotype from the original data set and \code{"inds"} will label the nodes 
with the representative individual names.}

\item{gscale}{"grey scale". If this is \code{TRUE}, this will scale the color
of the edges proportional to the observed distance, with the lines becoming
darker for more related nodes. See \code{\link{greycurve}} for details.}

\item{glim}{"grey limit". Two numbers between zero and one. They determine 
the upper and lower limits for the \code{\link{gray}} function. Default is 
0 (black) and 0.8 (20\% black). See \code{\link{greycurve}} for details.}

\item{gadj}{"grey adjust". a positive \code{integer} greater than zero that 
will serve as the exponent to the edge weight to scale the grey value to 
represent that weight. See \code{\link{greycurve}} for details.}

\item{gweight}{"grey weight". an \code{integer}. If it's 1, the grey scale 
will be weighted to emphasize the differences between closely related 
nodes. If it is 2, the grey scale will be weighted to emphasize the 
differences between more distantly related nodes. See 
\code{\link{greycurve}} for details.}

\item{wscale}{"width scale". If this is \code{TRUE}, the edge widths will be 
scaled proportional to the inverse of the observed distance , with the 
lines becoming thicker for more related nodes.}

\item{showplot}{logical. If \code{TRUE}, the graph will be plotted. If 
\code{FALSE}, it will simply be returned.}

\item{include.ties}{logical. If \code{TRUE}, the graph will include all edges
that were arbitrarily passed over in favor of another edge of equal weight.
If \code{FALSE}, which is the default, one edge will be arbitrarily 
selected when two or more edges are tied, resulting in a pure minimum
spanning network.}

\item{threshold}{numeric. By default, this is \code{NULL}, which will have no
effect. Any threshold value passed to this argument will be used in
\code{\link{mlg.filter}} prior to creating the MSN. If you have a data set
that contains contracted MLGs, this argument will override the threshold in
the data set. See Details.}

\item{clustering.algorithm}{string. By default, this is \code{NULL}. If
\code{threshold = NULL}, this argument will have no effect. When supplied
with either "farthest_neighbor", "average_neighbor", or "nearest_neighbor",
it will be passed to \code{\link{mlg.filter}} prior to creating the MSN. If
you have a data set that contains contracted MLGs, this argument will
override the algorithm in the data set. See Details.}

\item{...}{any other arguments that could go into plot.igraph}
}
\value{
\item{graph}{a minimum spanning network with nodes corresponding to 
  MLGs within the data set. Colors of the nodes represent population 
  membership. Width and color of the edges represent distance.} 
  \item{populations}{a vector of the population names corresponding to the 
  vertex colors} \item{colors}{a vector of the hexadecimal representations of
  the colors used in the vertex colors}
}
\description{
Create minimum spanning network of selected populations using Bruvo's 
distance.
}
\details{
The minimum spanning network generated by this function is generated
  via igraph's \code{\link[igraph:mst]{minimum.spanning.tree}}. The resultant
  graph produced can be plotted using igraph functions, or the entire object
  can be plotted using the function \code{\link{plot_poppr_msn}}, which will
  give the user a scale bar and the option to layout your data.
  \subsection{node sizes}{
  The area of the nodes are representative of the number of samples. Because
  \pkg{igraph} scales nodes by radius, the node sizes in the graph are 
  represented as the square root of the number of samples.}
  \subsection{mlg.compute}{
  Each node on the graph represents a different multilocus genotype. 
  The edges on the graph represent genetic distances that connect the
  multilocus genotypes. In genclone objects, it is possible to set the
  multilocus genotypes to a custom definition. This creates a problem for
  clone correction, however, as it is very possible to define custom lineages
  that are not monophyletic. When clone correction is performed on these
  definitions, information is lost from the graph. To circumvent this, The
  clone correction will be done via the computed multilocus genotypes, either
  "original" or "contracted". This is specified in the \code{mlg.compute}
  argument, above.}
  \subsection{contracted multilocus genotypes}{
  If your incoming data set is of the class \code{\linkS4class{genclone}},
  and it contains contracted multilocus genotypes, this function will retain
  that information for creating the minimum spanning network. You can use the
  arguments \code{threshold} and \code{clustering.algorithm} to change the
  threshold or clustering algorithm used in the network. For example, if you
  have a data set that has a threshold of 0.1 and you wish to have a minimum
  spanning network without a threshold, you can simply add 
  \code{threshold = 0.0}, and no clustering will happen. 
  
  The \code{threshold} and \code{clustering.algorithm} arguments can also be
  used to filter un-contracted data sets.
  }
}
\note{
\itemize{ \item \strong{Please see the documentation for
  \code{\link{bruvo.dist}} for details on the algorithm}. \item The edges of
  these graphs may cross each other if the graph becomes too large. \item The
  nodes in the graph represent multilocus genotypes. The colors of the nodes
  are representative of population membership. It is not uncommon to see
  different populations containing the same multilocus genotype.}
}
\examples{

# Load the data set.
data(nancycats)

# View populations 8 and 9 with default colors. 
bruvo.msn(nancycats, replen = rep(2, 9), sublist=8:9, vertex.label="inds", 
          vertex.label.cex=0.7, vertex.label.dist=0.4)
\dontrun{
# View heat colors.
bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
palette=heat.colors, vertex.label.cex=0.7, vertex.label.dist=0.4)

# View custom colors. Here, we use black and orange.
bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
vertex.label.dist=0.4)

# View with darker shades of grey (setting the upper limit to 1/2 black 1/2 white).
bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
vertex.label.dist=0.4, glim=c(0, 0.5))

# View with no grey scaling.
bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
vertex.label.dist=0.4, gscale=FALSE)

# View with no line widths.
bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
vertex.label.dist=0.4, wscale=FALSE)

# View with no scaling at all.
bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
vertex.label.dist=0.4, gscale=FALSE)

# View the whole population, but without labels.
bruvo.msn(nancycats, replen=rep(2, 9), vertex.label=NA)
}
}
\references{
Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza, and 
  Hinrich Schulenburg. A simple method for the calculation of microsatellite 
  genotype distances irrespective of ploidy level. Molecular Ecology, 
  13(7):2101-2106, 2004.
}
\seealso{
\code{\link{bruvo.dist}}, \code{\link[adegenet]{nancycats}}, 
  \code{\link{plot_poppr_msn}}, \code{\link[igraph]{mst}} 
  \code{\link{bruvo.boot}}, \code{\link{greycurve}} \code{\link{poppr.msn}}
}
\author{
Zhian N. Kamvar, Javier F. Tabima
}
