#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#
# This software was authored by Zhian N. Kamvar and Javier F. Tabima, graduate 
# students at Oregon State University; and Dr. Nik Grünwald, an employee of 
# USDA-ARS.
#
# Permission to use, copy, modify, and distribute this software and its
# documentation for educational, research and non-profit purposes, without fee, 
# and without a written agreement is hereby granted, provided that the statement
# above is incorporated into the material, giving appropriate attribution to the
# authors.
#
# Permission to incorporate this software into commercial products may be
# obtained by contacting USDA ARS and OREGON STATE UNIVERSITY Office for 
# Commercialization and Corporate Development.
#
# The software program and documentation are supplied "as is", without any
# accompanying services from the USDA or the University. USDA ARS or the 
# University do not warrant that the operation of the program will be 
# uninterrupted or error-free. The end-user understands that the program was 
# developed for research purposes and is advised not to rely exclusively on the 
# program for any reason.
#
# IN NO EVENT SHALL USDA ARS OR OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY 
# FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
# LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, 
# EVEN IF THE OREGON STATE UNIVERSITY HAS BEEN ADVISED OF THE POSSIBILITY OF 
# SUCH DAMAGE. USDA ARS OR OREGON STATE UNIVERSITY SPECIFICALLY DISCLAIMS ANY 
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND ANY STATUTORY 
# WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
# BASIS, AND USDA ARS AND OREGON STATE UNIVERSITY HAVE NO OBLIGATIONS TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#==============================================================================#
# BRUVO'S DISTANCE
#
# This distance is already applied in polysat, but that is more appropriate for
# polyploids. This is being implemented here so that the user does not have to 
# convert their data in order to perform the analysis.
#==============================================================================#
#==============================================================================#
#
#' Calculate the average Bruvo's Distance over all loci in a population.
#'
#' @param pop a \code{\link{genind}} object
#'
#' @param replen a \code{vector} of \code{integers} indicating the length of the
#' nucleotide repeats for each microsatellite locus.
#' 
#' @param add if \code{TRUE}, genotypes with zero values will be treated under
#' the genome addition model presented in Bruvo et al. 2004.
#' 
#' @param loss if \code{TRUE}, genotypes with zero values will be treated under
#' the genome loss model presented in Bruvo et al. 2004.
#'
#' @return a \code{distance matrix}
#'
#' @seealso \code{\link{nancycats}}
#'
#' @note The result of both \code{add = TRUE} and \code{loss = TRUE} is that
#' the distance is averaged over both values. If both are set to \code{FALSE},
#' then the infinite alleles model is used. For genotypes with all missing
#' values, the result will be NA. 
#'
#' If the user does not provide a vector of appropriate length for \code{replen}
#' , it will be estimated by taking the minimum difference among represented
#' alleles at each locus. It is not recommended to rely on this estimation. 
#'
#' @export
#' @author Zhian N. Kamvar
#' @references Bruvo, R., Michiels, N. K., D'Sousa, T. G., and Schulenberg, H. (2004) A simple method for calculation of microsatellite genotypes irrespective of ploidy level. \emph{Molecular Ecology} \strong{13}, 2101–2106.
#' @examples
#' # Please note that the data presented is assuming that the nancycat dataset 
#' # contains all dinucleotide repeats, it most likely is not an accurate
#' # representation of the data.
#'
#' # Load the nancycats dataset and construct the repeat vector.
#' data(nancycats)
#' ssr <- rep(2, 9)
#' 
#' # Analyze the 1st population in nancycats
#'
#' bruvo.dist(popsub(nancycats, 1), replen = ssr)
#'
#' # View each population as a heatmap.
#' \dontrun{
#' sapply(nancycats$pop.names, function(x) 
#' heatmap(as.matrix(bruvo.dist(popsub(nancycats, x), replen = ssr)), symm=TRUE))
#' }
#==============================================================================#
#' @useDynLib poppr
bruvo.dist <- function(pop, replen = 1, add = TRUE, loss = TRUE){
  # This attempts to make sure the data is true microsatellite data. It will
  # reject snp and aflp data. 
  if(pop@type != "codom" | all(is.na(unlist(lapply(pop@all.names, as.numeric))))){
    stop("\nThis dataset does not appear to be microsatellite data. Bruvo's Distance can only be applied for true microsatellites.")
  }
  ploid <- ploidy(pop)
  # Bruvo's distance depends on the knowledge of the repeat length. If the user
  # does not provide the repeat length, it can be estimated by the smallest
  # repeat difference greater than 1. This is not a preferred method. 
  if (length(replen) != length(pop@loc.names)){
    replen <- vapply(pop@all.names, function(x) guesslengths(as.numeric(x)), 1)
    #		replen <- rep(replen[1], numLoci)
    warning("\n\nRepeat length vector for loci is not equal to the number of loci represented.\nEstimating repeat lengths from data:\n", immediate.=TRUE)
    cat(replen,"\n\n")
  }
  popcols <- length(pop@loc.names)*ploid
  indnames <- pop@ind.names
  
  # This controlls for the user correcting missing data using "mean". 
  if(any(!pop@tab %in% c(0,((1:ploid)/ploid),1, NA))){
    pop@tab[!pop@tab %in% c(0,((1:ploid)/ploid),1, NA)] <- NA
  }
  if(any(rowSums(pop@tab, na.rm=TRUE) < nLoc(pop))){
    
    pop <- as.matrix(genind2df(pop, sep="/", usepop=FALSE))
    pop[pop %in% c("", NA)] <- paste(rep(0, ploid), collapse="/")
    return(phylo.bruvo.dist(pop, replen=replen, ploid=ploid, add = add, loss = loss))
  }
  else{
    pop <- matrix(as.numeric(as.matrix(genind2df(
  	  pop, oneColPerAll=TRUE, usepop=F))), ncol=popcols)
  }

  # Setting all missing data to 0.
  pop[is.na(pop)] <- 0
  # Dividing the data by the repeat length of each locus.

  pop <- pop / rep(replen, each=ploid*nrow(pop))
  pop <- matrix(round(pop), ncol=popcols)
  # Getting the permutation vector.
  perms <- .Call("permuto", ploid)
  # Calculating bruvo's distance over each locus. 
  distmat <- .Call("bruvo_distance", pop, perms, ploid, add, loss)
  # If there are missing values, the distance returns 100, which means that the
  # comparison is not made. These are changed to NA.
  distmat[distmat == 100] <- NA
  # Obtaining the average distance over all loci.
  avg.dist.vec <- apply(distmat, 1, mean, na.rm=TRUE)
  # presenting the information in a lower triangle distance matrix.
  dist.mat <- matrix(ncol=nrow(pop), nrow=nrow(pop))
  dist.mat[which(lower.tri(dist.mat)==TRUE)] <- avg.dist.vec
  dist.mat <- as.dist(dist.mat)
  attr(dist.mat, "Labels") <- indnames
  attr(dist.mat, "method") <- "Bruvo"
  attr(dist.mat, "call") <- match.call()
  return(dist.mat)
}


#==============================================================================#
#
#' Create a tree using Bruvo's Distance with non-parametric bootstrapping.
#'
#' @param pop a \code{\link{genind}} object
#'
#' @param replen a \code{vector} of \code{integers} indicating the length of the
#' nucleotide repeats for each microsatellite locus.
#' 
#' @param add if \code{TRUE}, genotypes with zero values will be treated under
#' the genome addition model presented in Bruvo et al. 2004.
#' 
#' @param loss if \code{TRUE}, genotypes with zero values will be treated under
#' the genome loss model presented in Bruvo et al. 2004.
#'
#' @param sample an \code{integer} indicated the number of bootstrap replicates
#' desired.
#'
#' @param tree choose between "nj" for neighbor-joining and "upgma" for a upgma
#' tree to be produced.
#'
#' @param showtree \code{logical} if \code{TRUE}, a tree will be plotted with
#' nodelabels.
#' 
#' @param cutoff \code{integer} the cutoff value for bootstrap node label values 
#' (between 0 and 100).
#' 
#' @param quiet \code{logical} defaults to \code{FALSE}. If \code{TRUE}, a
#' progress bar and messages will be supressed.
#' 
#' @param ... any argument to be passed on to \code{\link{boot.phylo}}. eg.
#' \code{quiet = TRUE}.
#'
#' @return a tree with nodelables
#'
#' @seealso \code{\link{nancycats}}, \code{\link{upgma}}, \code{\link{nj}},
#' \code{\link{boot.phylo}}, \code{\link{nodelabels}}, \code{\link{na.replace}},
#' \code{\link{missingno}}.
#'
#' @note This function calculates bruvo's distance for non-special cases (ie.
#' the ploidy and all alleles are known). Currently there is no way to import
#' polyploid partial heterozygote data into adegenet. For Bruvo's Distance 
#' concerning special cases, see the package \code{polysat}. 
#' Missing data is ignored, but be sure that missing data is NOT set to 0 in the
#' genind object. This is not easy to detect and will result in an error. Please
#' use any other method in \code{\link{na.replace}} or \code{\link{missingno}}. 
#'
#' If the user does not provide a vector of appropriate length for \code{replen}
#' , it will be estimated by taking the minimum difference among represented
#' alleles at each locus. It is not recommended to rely on this estimation. 
#'
#' @export
#' @author Javier F. Tabima, Zhian N. Kamvar
#' @references Bruvo, R., Michiels, N. K., D'Sousa, T. G., and Schulenberg, H. (2004) A simple method for calculation of microsatellite genotypes irrespective of ploidy level. \emph{Molecular Ecology} \strong{13}, 2101–2106.
#' @examples
#' # Please note that the data presented is assuming that the nancycat dataset 
#' # contains all dinucleotide repeats, it most likely is not an accurate
#' # representation of the data.
#'
#' # Load the nancycats dataset and construct the repeat vector.
#' data(nancycats)
#' ssr <- rep(1,9)
#' 
#' # Analyze the 1st population in nancycats
#'
#' bruvo.boot(popsub(nancycats, 1), replen=ssr)
#'
#==============================================================================#
#' @importFrom phangorn upgma
#' @importFrom ape nodelabels nj boot.phylo plot.phylo
#   /     \
#   |=(o)=|
#   \     /
bruvo.boot <- function(pop, replen = 1, add = TRUE, loss = TRUE, sample = 100, 
                        tree = "upgma", showtree = TRUE, cutoff = NULL, 
                        quiet = FALSE, ...){
  # This attempts to make sure the data is true microsatellite data. It will
  # reject snp and aflp data. 
  if(pop@type != "codom" | all(is.na(unlist(lapply(pop@all.names, as.numeric))))){
    stop("\nThis dataset does not appear to be microsatellite data. Bruvo's Distance can only be applied for true microsatellites.")
  }
  ploid <- ploidy(pop)
  # Bruvo's distance depends on the knowledge of the repeat length. If the user
  # does not provide the repeat length, it can be estimated by the smallest
  # repeat difference greater than 1. This is not a preferred method. 
  if (length(replen) != length(pop@loc.names)){
    replen <- vapply(pop@all.names, function(x) guesslengths(as.numeric(x)), 1)
    #    replen <- rep(replen[1], numLoci)
    warning("\n\nRepeat length vector for loci is not equal to the number of loci represented.\nEstimating repeat lengths from data:\n", immediate.=TRUE)
    cat(replen,"\n\n")
  }
  # This controlls for the user correcting missing data using "mean". 
  if(any(!pop@tab %in% c(0,((1:ploid)/ploid),1, NA))){
    pop@tab[!pop@tab %in% c(0,((1:ploid)/ploid),1, NA)] <- NA
  }
  # Converting the genind object into a matrix with each allele separated by "/"
  bar <- as.matrix(genind2df(pop, sep="/", usepop=FALSE))
  # The bruvo algorithm will ignore missing data, coded as 0.
  bar[bar %in% c("", NA)] <- paste(rep(0, ploid), collapse="/")
  #stopifnot(require(phangorn))
  # Steps: Create initial tree and then use boot.phylo to perform bootstrap
  # analysis, and then place the support labels on the tree.
  if(tree == "upgma"){
    # require(phangorn)
    newfunk <- match.fun(upgma)
  }
  else if(tree == "nj"){
    newfunk <- match.fun(nj)
  }
  tre <- newfunk(phylo.bruvo.dist(bar, replen=replen, ploid=ploid, add = add, loss = loss))
  if (any (tre$edge.length < 0)){
    warning("The branch lengths of the tree are negative.", immediate.=TRUE)
    us.promp<- as.numeric(readline(prompt="What do you want to do?, You can:\n1. Stop the analysis\n2. Convert negative branch values to Zero (Not biologically relevant)\nEnter your Selection:"))
    if (us.promp == 1){
      cat("Analysis stopped due to negative branches in the resulting tree\n")
      return(tre)
      stop()
    }
    else if (us.promp == 2){
      tre$edge.length[tre$edge.length < 0] <- 0
    }
    else{
      cat("Non-valid option. Returning tree.\n")
      return(tre)
    }
  }
  if(quiet == FALSE){
    cat("\nBootstrapping... (note: calculation of node labels can take a while even after the progress bar is full)\n\n")
  }
  bp <- boot.phylo(tre, bar, FUN = function (x) newfunk(phylo.bruvo.dist(x, replen = replen, ploid = ploid, add = add, loss = loss)), B = sample, quiet=quiet, ...)
  tre$node.labels <- round(((bp / sample)*100))
  if (!is.null(cutoff)){
    if (cutoff < 1 | cutoff > 100){
      cat("Cutoff value must be between 0 and 100.\n")
      cutoff<- as.numeric(readline(prompt = "Choose a new cutoff value between 0 and 100:\n"))
    }
    tre$node.labels[tre$node.labels < cutoff] <- NA
  }
  tre$tip.label <- pop@ind.names
  if(showtree == TRUE){
    plot(tre, show.node.label=TRUE)
  }
#  This works, but unfortunately, it lies ouside of the plot borders.
#  tre$tip.label <- pop@ind.names
#  if(showtree == TRUE){
#    plot(tre, show.tip.label=FALSE)
#    tiplabels(pop@ind.names, adj= c(-0.25, 0.5), frame="n", cex=0.8, font=2)
#    nodelabels(tre$node.label, adj = c(1.3, -0.5), frame="n", cex=0.9, font=3)
#  }

  if(tree=="upgma"){
    axisPhylo(3)
  }
  return(tre)
}
#==============================================================================#
#
#' Create minimum spanning network of selected populations using Brvuo's 
#' distance.
#'
#' @param pop a \code{\link{genind}} object
#'
#' @param replen a \code{vector} of \code{integers} indicating the length of the
#' nucleotide repeats for each microsatellite locus.
#'
#' @param add if \code{TRUE}, genotypes with zero values will be treated under
#' the genome addition model presented in Bruvo et al. 2004.
#' 
#' @param loss if \code{TRUE}, genotypes with zero values will be treated under
#' the genome loss model presented in Bruvo et al. 2004.
#'  
#' @param palette a \code{function} defining the color palette to be used to
#' color the populations on the graph. It defaults to \code{\link{topo.colors}},
#' but you can easily create new schemes by using \code{\link{colorRampPalette}}
#' (see examples for details)
#' 
#' @param sublist a \code{vector} of population names or indexes that the user
#' wishes to keep. Default to "ALL".
#'
#' @param blacklist a \code{vector} of population names or indexes that the user
#' wishes to discard. Default to \code{NULL}
#' 
#' @param vertex.label a \code{vector} of characters to label each vertex. There
#' are two defaults: \code{"MLG"} will label the nodes with the multilocus genotype
#' from the original data set and \code{"inds"} will label the nodes with the 
#' representative individual names.
#' 
#' @param gscale "grey scale". If this is \code{TRUE}, this will scale the color
#' of the edges proportional to Bruvo's distance, with the lines becoming darker
#' for more related nodes. See \code{\link{greycurve}} for details.
#' 
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#' the upper and lower limits for the \code{\link{gray}} function. Default is 0
#' (black) and 0.8 (20\% black). See \code{\link{greycurve}} for details.
#' 
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that
#' will serve as the exponent to the edge weight to scale the grey value to
#' represent that weight. See \code{\link{greycurve}} for details.
#' 
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale
#' will be weighted to emphasize the differences between closely related nodes.
#' If it is 2, the grey scale will be weighted to emphasize the differences
#' between more distantly related nodes. See \code{\link{greycurve}} for details.
#' 
#' @param wscale "width scale". If this is \code{TRUE}, the edge widths will be
#' scaled proportional to the inverse of Bruvo's distance , with the lines 
#' becoming thicker for more related nodes.
#' 
#' @param ... any other arguments that could go into plot.igraph
#'
#' @return 
#' \item{graph}{a minimum spanning network with nodes corresponding to MLGs within
#' the data set. Colors of the nodes represent population membership. Width and
#' color of the edges represent distance.}
#' \item{populations}{a vector of the population names corresponding to the 
#' vertex colors}
#' \item{colors}{a vector of the hexadecimal representations of the colors used
#' in the vertex colors}
#' 
#' @note The edges of these graphs may cross each other if the graph becomes too
#' large. 
#'
#' @seealso \code{\link{nancycats}}, \code{\link{upgma}}, \code{\link{nj}},
#' \code{\link{boot.phylo}}, \code{\link{nodelabels}}, \code{\link{na.replace}},
#' \code{\link{missingno}}, \code{\link{bruvo.boot}}, \code{\link{greycurve}}.
#' 
#' @export
#' @author Javier F. Tabima, Zhian N. Kamvar
#' @references Bruvo, R., Michiels, N. K., D'Sousa, T. G., and Schulenberg, H. (2004) A simple method for calculation of microsatellite genotypes irrespective of ploidy level. \emph{Molecular Ecology} \strong{13}, 2101–2106.
#' @examples
#' 
#' # Load the data set.
#' data(nancycats)
#' 
#' # View populations 8 and 9 with default colors. 
#' bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds", 
#' vertex.label.cex=0.7, vertex.label.dist=0.4)
#' \dontrun{
#' # View heat colors.
#' bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds", 
#' palette=heat.colors, vertex.label.cex=0.7, vertex.label.dist=0.4)
#' 
#' # View custom colors. Here, we use black and orange.
#' bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4)
#' 
#' # View with darker shades of grey (setting the upper limit to 1/2 black 1/2 white).
#' bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4, glim=c(0, 0.5))
#' 
#' # View with no grey scaling.
#' bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4, gscale=FALSE)
#' 
#' # View with no line widths.
#' bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4, wscale=FALSE)
#' 
#' # View with no scaling at all.
#' bruvo.msn(nancycats, replen=rep(1, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black"), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4, vscale=FALSE, gscale=FALSE)
#' 
#' # View the whole population, but without labels.
#' bruvo.msn(nancycats, replen=rep(1, 9), vertex.label=NA)
#' }
#==============================================================================#
#' @importFrom igraph graph.adjacency plot.igraph V E minimum.spanning.tree V<- E<-
bruvo.msn <- function (pop, replen = 1, add = TRUE, loss = TRUE, palette = topo.colors,
                       sublist = "All", blacklist = NULL, vertex.label = "MLG", 
                       gscale = TRUE, glim = c(0,0.8), gadj = 3, gweight = 1, 
                       wscale = TRUE, ...){
  # require(igraph)
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  # Storing the MLG vector into the genind object
  pop$other$mlg.vec <- mlg.vector(pop)
  
  singlepop <- function(pop, vertex.label){
    cpop <- pop[.clonecorrector(pop), ]
    mlg.number <- table(pop$other$mlg.vec)[rank(cpop$other$mlg.vec)]
    bclone <- bruvo.dist(cpop, replen=replen, add = add, loss = loss)
    mclone <- as.dist(bclone)
    #attr(bclone, "Labels") <- paste("MLG.", cpop$other$mlg.vec, sep="")
    g <- graph.adjacency(as.matrix(bclone),weighted=TRUE,mode="undirected")
    mst <- (minimum.spanning.tree(g,algorithm="prim",weights=E(g)$weight))
    if(!is.na(vertex.label[1]) & length(vertex.label) == 1){
      if(toupper(vertex.label) == "MLG"){
        vertex.label <- paste("MLG.", cpop$other$mlg.vec, sep="")
      }
      else if(toupper(vertex.label) == "INDS"){
        vertex.label <- cpop$ind.names
      }
    }
    
    if(gscale == TRUE){
      E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correction=gadj, 
                                       show=FALSE))
    }
    else{
      E(mst)$color <- rep("black", length(E(mst)$weight))
    }
    
    edgewidth <- 2
    if(wscale == TRUE){
      edgewidth <- 1/(E(mst)$weight)
      if(any(E(mst)$weight < 0.08)){
        edgewidth <- 1/(E(mst)$weight + 0.08)
      }
    }
    populations <- ifelse(is.null(pop(pop)), NA, pop$pop.names)
    plot(mst, edge.width = edgewidth, edge.color = E(mst)$color,  
         vertex.label = vertex.label, vertex.size = mlg.number*3, 
         vertex.color = palette(1),  ...)
    legend(-1.55,1,bty = "n", cex = 0.75, 
           legend = populations, title = "Populations", fill = palette(1), 
           border = NULL)
    E(mst)$width <- edgewidth
    V(mst)$size <- mlg.number
    V(mst)$color <- palette(1)
    V(mst)$label <- vertex.label
    return(list(graph = mst, populations = populations, colors = palette(1)))
  }
  if(is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop(pop, vertex.label))
  }
  if(sublist[1] != "ALL" | !is.null(blacklist)){
    pop <- popsub(pop, sublist, blacklist)
  }
  if(is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop(pop, vertex.label))
  }
  # Obtaining population information for all MLGs
  mlg.cp <- mlg.crosspop(pop, mlgsub=1:mlg(pop, quiet=TRUE), quiet=TRUE)
  names(mlg.cp) <- paste("MLG.",sort(unique(pop$other$mlg.vec)),sep="")
  cpop <- pop[.clonecorrector(pop), ]
  # This will determine the size of the nodes based on the number of individuals
  # in the MLG. Subsetting by the MLG vector of the clone corrected set will
  # give us the numbers and the population information in the correct order.
  # Note: rank is used to correctly subset the data
  mlg.number <- table(pop$other$mlg.vec)[rank(cpop$other$mlg.vec)]
  mlg.cp <- mlg.cp[rank(cpop$other$mlg.vec)]
  bclone <- bruvo.dist(cpop, replen=replen, add = add, loss = loss)
  
  ###### Create a graph #######
  g <- graph.adjacency(as.matrix(bclone), weighted = TRUE, mode = "undirected")
  mst <- (minimum.spanning.tree(g, algorithm = "prim", weights = E(g)$weight))
  
  if(!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if(toupper(vertex.label) == "MLG"){
      vertex.label <- paste("MLG.", cpop$other$mlg.vec, sep="")
    }
    else if(toupper(vertex.label) == "INDS"){
      vertex.label <- cpop$ind.names
    }
  }
  ###### Color schemes #######  
  # The pallete is determined by what the user types in the argument. It can be 
  # rainbow, topo.colors, heat.colors ...etc.
  palette <- match.fun(palette)
  color <- palette(length(pop@pop.names))
  if(gscale == TRUE){
    E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correction=gadj, 
                                     show=FALSE))
  }
  else{
    E(mst)$color <- rep("black", length(E(mst)$weight))
  }
  
  edgewidth <- 2
  if(wscale == TRUE){
    edgewidth <- 1/(E(mst)$weight)
    if(any(E(mst)$weight < 0.08)){
      edgewidth <- 1/(E(mst)$weight + 0.08)
    }
  }
  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pop@pop.names %in% names(x)])
  
  plot(mst, edge.width = edgewidth, edge.color = E(mst)$color, 
       vertex.size = mlg.number*3, vertex.shape = "pie", vertex.pie = mlg.cp, 
       vertex.pie.color = mlg.color, vertex.label = vertex.label, ...)
  legend(-1.55 ,1 ,bty = "n", cex = 0.75, legend = pop$pop.names, 
         title = "Populations", fill = color, border = NULL)
  E(mst)$width <- edgewidth
  V(mst)$size <- mlg.number
  V(mst)$shape <- "pie"
  V(mst)$pie <- mlg.cp
  V(mst)$pie.color <- mlg.color
  V(mst)$label <- vertex.label
  return(list(graph = mst, populations = pop$pop.names, colors = color))
}
