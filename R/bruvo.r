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
#' Bruvo's distance for microsatellites
#' 
#' Calculate the average Bruvo's distance over all loci in a population.
#' 
#' @param pop a \code{\link{genind}} object
#'   
#' @param replen a \code{vector} of \code{integers} indicating the length of the
#'   nucleotide repeats for each microsatellite locus.
#'   
#' @param add if \code{TRUE}, genotypes with zero values will be treated under 
#'   the genome addition model presented in Bruvo et al. 2004.
#'   
#' @param loss if \code{TRUE}, genotypes with zero values will be treated under 
#'   the genome loss model presented in Bruvo et al. 2004.
#'   
#' @return a \code{distance matrix}
#'   
#' @note The result of both \code{add = TRUE} and \code{loss = TRUE} is that the
#'   distance is averaged over both values. If both are set to \code{FALSE}, 
#'   then the infinite alleles model is used. For genotypes with all missing 
#'   values, the result will be NA.
#'   
#'   If the user does not provide a vector of appropriate length for 
#'   \code{replen} , it will be estimated by taking the minimum difference among
#'   represented alleles at each locus. IT IS NOT RECOMMENDED TO RELY ON THIS 
#'   ESTIMATION.
#'   
#' @details Ploidy is irrelevant with respect to calculation of Bruvo's 
#'   distance. However, since it makes a comparison between all alleles at a 
#'   locus, it only makes sense that the two loci need to have the same ploidy 
#'   level. Unfortunately for polyploids, it's often difficult to fully separate
#'   distinct alleles at each locus, so you end up with genotypes that appear to
#'   have a lower ploidy level than the organism.
#'   
#'   To help deal with these situations, Bruvo has suggested three methods for 
#'   dealing with these differences in ploidy levels: \itemize{ \item Infinite 
#'   Model - The simplest way to deal with it is to count all missing alleles as
#'   infinitely large so that the distance between it and anything else is 1. 
#'   Aside from this being computationally simple, it will tend to 
#'   \strong{inflate distances between individuals}. \item Genome Addition Model
#'   - If it is suspected that the organism has gone through a recent genome 
#'   expansion, \strong{the missing alleles will be replace with all possible 
#'   combinations of the observed alleles in the shorter genotype}. For example,
#'   if there is a genotype of [69, 70, 0, 0] where 0 is a missing allele, the 
#'   possible combinations are: [69, 70, 69, 69], [69, 70, 69, 70], and [69, 70,
#'   70, 70]. The resulting distances are then averaged over the number of 
#'   comparisons. \item Genome Loss Model - This is similar to the genome 
#'   addition model, except that it assumes that there was a recent genome 
#'   reduction event and uses \strong{the observed values in the full genotype 
#'   to fill the missing values in the short genotype}. As with the Genome 
#'   Addition Model, the resulting distances are averaged over the number of 
#'   comparisons. \item Combination Model - Combine and average the genome 
#'   addition and loss models. } As mentioned above, the infinite model is 
#'   biased, but it is not nearly as computationally intensive as either of the 
#'   other models. The reason for this is that both of the addition and loss 
#'   models requires replacement of alleles and recalculation of Bruvo's 
#'   distance. The number of replacements required is equal to the multiset 
#'   coefficient: \eqn{\left({n \choose k}\right) == {(n+k-1) \choose 
#'   k}}{choose(n+k-1, k)} where \emph{n} is the number of potential 
#'   replacements and \emph{k} is the number of alleles to be replaced. So, for 
#'   the example given above, The genome addition model would require 
#'   \eqn{\left({2 \choose 2}\right) = 3}{choose(2+2-1, 2) == 3} calculations of
#'   Bruvo's distance, whereas the genome loss model would require \eqn{\left({4
#'   \choose 2}\right) = 10}{choose(4+2-1, 2) == 10} calculations.
#'   
#'   To reduce the number of calculations and assumptions otherwise, Bruvo's 
#'   distance will be calculated using the largest observed ploidy in pairwise
#'   comparisons. This means that when comparing [69,70,71,0] and [59,60,0,0],
#'   they will be treated as triploids.
#'   
#' @export
#' @author Zhian N. Kamvar
#'   
#' @references Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza, and 
#'   Hinrich Schulenburg. A simple method for the calculation of microsatellite 
#'   genotype distances irrespective of ploidy level. Molecular Ecology, 
#'   13(7):2101-2106, 2004.
#'   
#' @seealso \code{\link{bruvo.boot}}, \code{\link{bruvo.msn}}
#'   
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
  if (pop@type != "codom" | all(is.na(unlist(lapply(pop@all.names, as.numeric))))){
    stop(non_ssr_data_warning())
  }
  # Bruvo's distance depends on the knowledge of the repeat length. If the user
  # does not provide the repeat length, it can be estimated by the smallest
  # repeat difference greater than 1. This is not a preferred method. 
  if (length(replen) != length(pop@loc.names)){
    replen <- vapply(pop@all.names, function(x) guesslengths(as.numeric(x)), 1)
    warning(repeat_length_warning(replen), immediate. = TRUE)
  }
  bruvomat  <- new('bruvomat', pop, replen)
  funk_call <- match.call()
  dist.mat  <- bruvos_distance(bruvomat, funk_call = funk_call, add, loss)
  return(dist.mat)
}


#==============================================================================#
#
#' Create a tree using Bruvo's Distance with non-parametric bootstrapping.
#' 
#' @param pop a \code{\link{genind}} object
#'   
#' @param replen a \code{vector} of \code{integers} indicating the length of the
#'   nucleotide repeats for each microsatellite locus.
#'   
#' @param add if \code{TRUE}, genotypes with zero values will be treated under 
#'   the genome addition model presented in Bruvo et al. 2004.
#'   
#' @param loss if \code{TRUE}, genotypes with zero values will be treated under 
#'   the genome loss model presented in Bruvo et al. 2004.
#'   
#' @param sample an \code{integer} indicated the number of bootstrap replicates 
#'   desired.
#'   
#' @param tree choose between "nj" for neighbor-joining and "upgma" for a upgma 
#'   tree to be produced.
#'   
#' @param showtree \code{logical} if \code{TRUE}, a tree will be plotted with 
#'   nodelabels.
#'   
#' @param cutoff \code{integer} the cutoff value for bootstrap node label values
#'   (between 0 and 100).
#'   
#' @param quiet \code{logical} defaults to \code{FALSE}. If \code{TRUE}, a 
#'   progress bar and messages will be suppressed.
#'   
#' @param ... any argument to be passed on to \code{\link{boot.phylo}}. eg. 
#'   \code{quiet = TRUE}.
#'   
#' @return a tree of class phylo with nodelables
#'   
#' @seealso \code{\link{bruvo.dist}}, \code{\link{nancycats}}, 
#'   \code{\link{upgma}}, \code{\link{nj}}, \code{\link{boot.phylo}}, 
#'   \code{\link{nodelabels}}, \code{\link{na.replace}}, 
#'   \code{\link{missingno}}.
#'   
#' @note \strong{Please refer to the documentation for bruvo.dist for details on
#'   the algorithm.} If the user does not provide a vector of appropriate length
#'   for \code{replen} , it will be estimated by taking the minimum difference
#'   among represented alleles at each locus. IT IS NOT RECOMMENDED TO RELY ON
#'   THIS ESTIMATION.
#'   
#' @export
#' @author Zhian N. Kamvar, Javier F. Tabima
#'   
#' @references Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza, and
#' Hinrich Schulenburg. A simple method for the calculation of microsatellite
#' genotype distances irrespective of ploidy level. Molecular Ecology,
#' 13(7):2101-2106, 2004.
#' 
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
#' bruvo.boot(popsub(nancycats, 1), replen = ssr)
#' 
#==============================================================================#
#' @importFrom phangorn upgma midpoint
#' @importFrom ape nodelabels nj boot.phylo plot.phylo axisPhylo ladderize 
#' @importFrom ape add.scale.bar nodelabels tiplabels
#   /     \
#   |=(o)=|
#   \     /
bruvo.boot <- function(pop, replen = 1, add = TRUE, loss = TRUE, sample = 100, 
                        tree = "upgma", showtree = TRUE, cutoff = NULL, 
                        quiet = FALSE, ...){
  # This attempts to make sure the data is true microsatellite data. It will
  # reject snp and aflp data. 
  if (pop@type != "codom" | all(is.na(unlist(lapply(pop@all.names, as.numeric))))){
    stop(non_ssr_data_warning())
  }
  # Bruvo's distance depends on the knowledge of the repeat length. If the user
  # does not provide the repeat length, it can be estimated by the smallest
  # repeat difference greater than 1. This is not a preferred method. 
  if (length(replen) != length(pop@loc.names)){
    replen <- vapply(pop@all.names, function(x) guesslengths(as.numeric(x)), 1)
    warning(repeat_length_warning(replen), immediate. = TRUE)
  }
  bootgen <- new('bruvomat', pop, replen)
  # Steps: Create initial tree and then use boot.phylo to perform bootstrap
  # analysis, and then place the support labels on the tree.
  if (tree == "upgma"){
    root <- TRUE
    newfunk <- match.fun(upgma)
  } else if (tree == "nj"){
    root <- FALSE
    newfunk <- match.fun(nj)
  }
  tre <- newfunk(bruvos_distance(bootgen, funk_call = match.call(), add = add, loss = loss))
  if (any (tre$edge.length < 0)){
    warning(negative_branch_warning(), immediate.=TRUE)
	  tre <- fix_negative_branch(tre)
  }
  if (quiet == FALSE){
    cat("\nBootstrapping...\n") 
    cat("(note: calculation of node labels can take a while even after") 
    cat(" the progress bar is full)\n\n")
  }
  bootfun <- function(x){
    return(newfunk(bruvos_distance(x, funk_call = match.call(), add = add, loss = loss)))
  }
  bp <- boot.phylo(tre, bootgen, FUN = bootfun, B = sample, quiet = quiet, 
                   rooted = root, ...)
  tre$node.labels <- round(((bp / sample)*100))
  if (!is.null(cutoff)){
    if (cutoff < 1 | cutoff > 100){
      cat("Cutoff value must be between 0 and 100.\n")
      prompt_msg <- "Choose a new cutoff value between 0 and 100:\n"
      cutoff <- as.numeric(readline(prompt = prompt_msg))
    }
    tre$node.labels[tre$node.labels < cutoff] <- NA
  }
  tre$tip.label <- pop@ind.names
  if (showtree == TRUE){
    poppr.plot.phylo(tre, tree)
  }
  return(tre)
}
#==============================================================================#
#
#' Create minimum spanning network of selected populations using Bruvo's 
#' distance.
#' 
#' @param pop a \code{\link{genind}} object
#'   
#' @param replen a \code{vector} of \code{integers} indicating the length of the
#'   nucleotide repeats for each microsatellite locus.
#'   
#' @param add if \code{TRUE}, genotypes with zero values will be treated under 
#'   the genome addition model presented in Bruvo et al. 2004.
#'   
#' @param loss if \code{TRUE}, genotypes with zero values will be treated under 
#'   the genome loss model presented in Bruvo et al. 2004.
#'   
#' @param palette a \code{function} defining the color palette to be used to 
#'   color the populations on the graph. It defaults to 
#'   \code{\link{topo.colors}}, but you can easily create new schemes by using 
#'   \code{\link{colorRampPalette}} (see examples for details)
#'   
#' @param sublist a \code{vector} of population names or indexes that the user 
#'   wishes to keep. Default to "ALL".
#'   
#' @param blacklist a \code{vector} of population names or indexes that the user
#'   wishes to discard. Default to \code{NULL}
#'   
#' @param vertex.label a \code{vector} of characters to label each vertex. There
#'   are two defaults: \code{"MLG"} will label the nodes with the multilocus 
#'   genotype from the original data set and \code{"inds"} will label the nodes 
#'   with the representative individual names.
#'   
#' @param gscale "grey scale". If this is \code{TRUE}, this will scale the color
#'   of the edges proportional to Bruvo's distance, with the lines becoming 
#'   darker for more related nodes. See \code{\link{greycurve}} for details.
#'   
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#'   the upper and lower limits for the \code{\link{gray}} function. Default is 
#'   0 (black) and 0.8 (20\% black). See \code{\link{greycurve}} for details.
#'   
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that 
#'   will serve as the exponent to the edge weight to scale the grey value to 
#'   represent that weight. See \code{\link{greycurve}} for details.
#'   
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale 
#'   will be weighted to emphasize the differences between closely related 
#'   nodes. If it is 2, the grey scale will be weighted to emphasize the 
#'   differences between more distantly related nodes. See 
#'   \code{\link{greycurve}} for details.
#'   
#' @param wscale "width scale". If this is \code{TRUE}, the edge widths will be 
#'   scaled proportional to Bruvo's distance, with the lines becoming thicker
#'   for more related nodes.
#'   
#' @param showplot logical. If \code{TRUE}, the graph will be plotted. If 
#'   \code{FALSE}, it will simply be returned.
#'
#' @param include.ties logical. If \code{TRUE}, the graph will include all
#'   edges that were arbitrarily passed over in favor of another edge of equal
#'   weight. If \code{FALSE}, which is the default, one edge will be arbitrarily 
#'   selected when two or more edges are tied, resulting in a pure minimum spanning
#'   network. 
#'   
#' @param ... any other arguments that could go into plot.igraph
#'   
#' @return \item{graph}{a minimum spanning network with nodes corresponding to 
#'   MLGs within the data set. Colors of the nodes represent population 
#'   membership. Width and color of the edges represent distance.} 
#'   \item{populations}{a vector of the population names corresponding to the 
#'   vertex colors} \item{colors}{a vector of the hexadecimal representations of
#'   the colors used in the vertex colors}
#'   
#' @note \itemize{ \item \strong{Please see the documentation for
#'   \code{\link{bruvo.dist}} for details on the algorithm}. \item The edges of
#'   these graphs may cross each other if the graph becomes too large. \item The
#'   nodes in the graph represent multilocus genotypes. The colors of the nodes
#'   are representative of population membership. It is not uncommon to see
#'   different populations containing the same multilocus genotype.}
#'   
#' @details The minimum spanning network generated by this function is generated
#'   via igraph's \code{\link[igraph]{minimum.spanning.tree}}. The resultant
#'   graph produced can be plotted using igraph functions, or the entire object
#'   can be plotted using the function \code{\link{plot_poppr_msn}}, which will
#'   give the user a scale bar and the option to layout your data.
#'   
#' @seealso \code{\link{bruvo.dist}}, \code{\link{nancycats}}, 
#'   \code{\link{plot_poppr_msn}}, \code{\link[igraph]{minimum.spanning.tree}} 
#'   \code{\link{bruvo.boot}}, \code{\link{greycurve}}
#'   
#' @export
#' @author Zhian N. Kamvar, Javier F. Tabima
#' @aliases msn.bruvo
#' @references Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza, and 
#'   Hinrich Schulenburg. A simple method for the calculation of microsatellite 
#'   genotype distances irrespective of ploidy level. Molecular Ecology, 
#'   13(7):2101-2106, 2004.
#'   
#' @examples
#' 
#' # Load the data set.
#' data(nancycats)
#' 
#' # View populations 8 and 9 with default colors. 
#' bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
#' vertex.label.cex=0.7, vertex.label.dist=0.4)
#' \dontrun{
#' # View heat colors.
#' bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
#' palette=heat.colors, vertex.label.cex=0.7, vertex.label.dist=0.4)
#' 
#' # View custom colors. Here, we use black and orange.
#' bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4)
#' 
#' # View with darker shades of grey (setting the upper limit to 1/2 black 1/2 white).
#' bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4, glim=c(0, 0.5))
#' 
#' # View with no grey scaling.
#' bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4, gscale=FALSE)
#' 
#' # View with no line widths.
#' bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4, wscale=FALSE)
#' 
#' # View with no scaling at all.
#' bruvo.msn(nancycats, replen=rep(2, 9), sublist=8:9, vertex.label="inds", 
#' palette = colorRampPalette(c("orange", "black")), vertex.label.cex=0.7, 
#' vertex.label.dist=0.4, vscale=FALSE, gscale=FALSE)
#' 
#' # View the whole population, but without labels.
#' bruvo.msn(nancycats, replen=rep(2, 9), vertex.label=NA)
#' }
#==============================================================================#
#' @importFrom igraph graph.adjacency plot.igraph V E minimum.spanning.tree V<- E<- print.igraph add.edges
bruvo.msn <- function (pop, replen = 1, add = TRUE, loss = TRUE, palette = topo.colors,
                       sublist = "All", blacklist = NULL, vertex.label = "MLG", 
                       gscale = TRUE, glim = c(0,0.8), gadj = 3, gweight = 1, 
                       wscale = TRUE, showplot = TRUE,
                       include.ties = FALSE,  ...){

  gadj <- ifelse(gweight == 1, gadj, -gadj)
  if (!is.genclone(pop)){
    # Storing the MLG vector into the genind object
    pop$other$mlg.vec <- mlg.vector(pop)  
  }
  if (is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop_msn(pop, vertex.label, replen = replen, gscale = gscale, 
                         glim = glim, gadj = gadj, wscale = wscale, 
                         palette = palette, include.ties = include.ties))
  }
  if (sublist[1] != "ALL" | !is.null(blacklist)){
    pop <- popsub(pop, sublist, blacklist)
  }
  if (is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop_msn(pop, vertex.label, replen = replen, gscale = gscale, 
                         glim = glim, gadj = gadj, wscale = wscale, 
                         palette = palette, showplot = showplot, include.ties = include.ties, ...))
  }
  # Obtaining population information for all MLGs
  cpop <- pop[.clonecorrector(pop), ]
  if (is.genclone(pop)){
    subs <- sort(unique(pop@mlg))
  } else {
    subs <- 1:mlg(pop, quiet = TRUE)
  }
  mlg.cp <- mlg.crosspop(pop, mlgsub = subs, quiet=TRUE)
  if (is.genclone(pop)){
    mlgs <- pop@mlg
    cmlg <- cpop@mlg
  } else {
    mlgs <- pop$other$mlg.vec
    cmlg <- cpop$other$mlg.vec
  }
  names(mlg.cp) <- paste0("MLG.", sort(unique(mlgs)))
  
  # This will determine the size of the nodes based on the number of individuals
  # in the MLG. Subsetting by the MLG vector of the clone corrected set will
  # give us the numbers and the population information in the correct order.
  # Note: rank is used to correctly subset the data
  mlg.number <- table(mlgs)[rank(cmlg)]
  mlg.cp     <- mlg.cp[rank(cmlg)]
  bclone     <- as.matrix(bruvo.dist(cpop, replen=replen))
  rownames(bclone) <- cpop$pop
  colnames(bclone) <- cpop$pop
  
  ###### Create a graph #######
  g   <- graph.adjacency(as.matrix(bclone), weighted = TRUE, mode = "undirected")
  mst <- minimum.spanning.tree(g, algorithm = "prim", weights = E(g)$weight)

  # Add any relevant edges that were cut from the mst while still being tied for the title of optimal edge
  if(include.ties){
    tied_edges <- .Call("msn_tied_edges",as.matrix(mst[]),as.matrix(bclone),(.Machine$double.eps ^ 0.5))
    if(length(tied_edges) > 0){
      mst <- add.edges(mst, dimnames(mst[])[[1]][tied_edges[c(TRUE,TRUE,FALSE)]], weight=tied_edges[c(FALSE,FALSE,TRUE)])
    }
  }

  if (!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if (toupper(vertex.label) == "MLG"){
      vertex.label <- paste0("MLG.", cmlg)
    } else if (toupper(vertex.label) == "INDS"){
      vertex.label <- cpop$ind.names
    }
  }
  ###### Color schemes #######  
  # The palette is determined by what the user types in the argument. It can be 
  # rainbow, topo.colors, heat.colors ...etc.
  palette <- match.fun(palette)
  color   <- palette(length(pop@pop.names))
  mst     <- update_edge_scales(mst, wscale, gscale, glim, gadj)

  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pop@pop.names %in% names(x)])
  if (showplot){
    plot.igraph(mst, edge.width = E(mst)$width, edge.color = E(mst)$color, 
         vertex.size = mlg.number*3, vertex.shape = "pie", vertex.pie = mlg.cp, 
         vertex.pie.color = mlg.color, vertex.label = vertex.label, ...)
    legend(-1.55, 1, bty = "n", cex = 0.75, legend = pop$pop.names, 
           title = "Populations", fill = color, border = NULL)
  }
  V(mst)$size      <- mlg.number
  V(mst)$shape     <- "pie"
  V(mst)$pie       <- mlg.cp
  V(mst)$pie.color <- mlg.color
  V(mst)$label     <- vertex.label
  return(list(graph = mst, populations = pop$pop.names, colors = color))
}
