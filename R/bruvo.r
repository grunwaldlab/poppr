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
#'   nucleotide repeats for each microsatellite locus. E.g. a locus with a (CAT) 
#'   repeat would have a replen value of 3. (Also see \code{\link{fix_replen}})
#'   
#' @param add if \code{TRUE}, genotypes with zero values will be treated under 
#'   the genome addition model presented in Bruvo et al. 2004. See the
#'   \strong{Note} section for options.
#'   
#' @param loss if \code{TRUE}, genotypes with zero values will be treated under 
#'   the genome loss model presented in Bruvo et al. 2004. See the
#'   \strong{Note} section for options.
#'   
#' @return an object of class \code{\link{dist}}
#'   
#' @details 
#'   Bruvo's distance between two alleles is calculated as 
#'   \deqn{d = 1 - 2^{-\mid x \mid}}{d = 1 - (2^(-abs(x)))}, where \strong{x}
#'   is the number of repeat units between the two alleles (see the Algorithms 
#'   and Equations vignette for more details). These distances are calculated 
#'   over all combinations of alleles at a locus and then the minimum average
#'   distance between allele combinations is taken as the distance for that 
#'   locus. All loci are then averaged over to obtain the distance between two
#'   samples. Missing data is ignored (in the same fashion as 
#'   \code{mean(c(1:9, NA), na.rm = TRUE)}) if all alleles are missing. See the
#'   next section for other cases.
#'   
#'   \subsection{Polyploids}{
#'   Ploidy is irrelevant with respect to calculation of Bruvo's 
#'   distance. However, since it makes a comparison between all alleles at a 
#'   locus, it only makes sense that the two loci need to have the same ploidy 
#'   level. Unfortunately for polyploids, it's often difficult to fully separate
#'   distinct alleles at each locus, so you end up with genotypes that appear to
#'   have a lower ploidy level than the organism.
#'   
#'   To help deal with these situations, Bruvo has suggested three methods for 
#'   dealing with these differences in ploidy levels: \itemize{ \item
#'   \strong{Infinite Model} - The simplest way to deal with it is to count all
#'   missing alleles as infinitely large so that the distance between it and
#'   anything else is 1. Aside from this being computationally simple, it will
#'   tend to \strong{inflate distances between individuals}. \item
#'   \strong{Genome Addition Model} - If it is suspected that the organism has
#'   gone through a recent genome expansion, \strong{the missing alleles will be
#'   replace with all possible combinations of the observed alleles in the
#'   shorter genotype}. For example, if there is a genotype of [69, 70, 0, 0]
#'   where 0 is a missing allele, the possible combinations are: [69, 70, 69,
#'   69], [69, 70, 69, 70], and [69, 70, 70, 70]. The resulting distances are
#'   then averaged over the number of comparisons. \item \strong{Genome Loss
#'   Model} - This is similar to the genome addition model, except that it
#'   assumes that there was a recent genome reduction event and uses \strong{the
#'   observed values in the full genotype to fill the missing values in the
#'   short genotype}. As with the Genome Addition Model, the resulting distances
#'   are averaged over the number of comparisons. \item \strong{Combination
#'   Model} - Combine and average the genome addition and loss models. } As
#'   mentioned above, the infinite model is biased, but it is not nearly as
#'   computationally intensive as either of the other models. The reason for
#'   this is that both of the addition and loss models requires replacement of
#'   alleles and recalculation of Bruvo's distance. The number of replacements
#'   required is equal to the multiset coefficient: \eqn{\left({n \choose
#'   k}\right) == {(n+k-1) \choose k}}{choose(n+k-1, k)} where \emph{n} is the
#'   number of potential replacements and \emph{k} is the number of alleles to
#'   be replaced. So, for the example given above, The genome addition model
#'   would require \eqn{\left({2 \choose 2}\right) = 3}{choose(2+2-1, 2) == 3}
#'   calculations of Bruvo's distance, whereas the genome loss model would
#'   require \eqn{\left({4 \choose 2}\right) = 10}{choose(4+2-1, 2) == 10}
#'   calculations.
#'   
#'   To reduce the number of calculations and assumptions otherwise, Bruvo's 
#'   distance will be calculated using the largest observed ploidy in pairwise 
#'   comparisons. This means that when comparing [69,70,71,0] and [59,60,0,0], 
#'   they will be treated as triploids.
#'   }
#'   
#' @note Do not use missingno with this function. 
#'   \subsection{Repeat Lengths (replen)}{
#'   The \code{replen} argument is crucial for proper analysis of Bruvo's
#'   distance since the calculation relies on the knowledge of the number of
#'   steps between alleles. To calculate Bruvo's distance, your raw allele calls
#'   are first divided by the repeat lengths and then rounded. This can create a
#'   problem with repeat lengths of even size due to the IEC 60559 standard that
#'   says rounding at 0.5 is to the nearest even number, meaning that it is
#'   possible for two alleles that are one step apart may appear to be exactly
#'   the same. This can be fixed by subtracting a tiny number from the repeat
#'   length with the function \code{\link{fix_replen}}. Please consider using
#'   this before running Bruvo's distance.
#'   }
#'   \subsection{Model Choice}{ The \code{add} and \code{loss} arguments 
#'   modify the model choice accordingly: \itemize{ \item \strong{Infinite 
#'   Model:}  \code{add = FALSE, loss = FALSE} \item \strong{Genome Addition 
#'   Model:}  \code{add = TRUE, loss = FALSE} \item \strong{Genome Loss Model:} 
#'   \code{add = FALSE, loss = TRUE} \item \strong{Combination Model}
#'   \emph{(DEFAULT):}  \code{add = TRUE, loss = TRUE} } Details of each model
#'   choice are described in the \strong{Details} section, above. Additionally,
#'   genotypes containing all missing values at a locus will return a value of
#'   \code{NA} and not contribute to the average across loci. }
#'   \subsection{Repeat Lengths}{ If the user does not provide a vector of 
#'   appropriate length for \code{replen} , it will be estimated by taking the 
#'   minimum difference among represented alleles at each locus. IT IS NOT 
#'   RECOMMENDED TO RELY ON THIS ESTIMATION. }
#'   
#' @export
#' @author Zhian N. Kamvar
#'   
#' @references Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza, and 
#'   Hinrich Schulenburg. A simple method for the calculation of microsatellite 
#'   genotype distances irrespective of ploidy level. Molecular Ecology, 
#'   13(7):2101-2106, 2004.
#'   
#' @seealso \code{\link{fix_replen}}, \code{\link{test_replen}},
#'   \code{\link{bruvo.boot}}, \code{\link{bruvo.msn}}
#'   
#' @examples
#' # Please note that the data presented is assuming that the nancycat dataset 
#' # contains all dinucleotide repeats, it most likely is not an accurate
#' # representation of the data.
#' 
#' # Load the nancycats dataset and construct the repeat vector.
#' data(nancycats)
#' 
#' # Assume the alleles are all dinucleotide repeats.
#' ssr <- rep(2, nLoc(nancycats))
#' test_replen(nancycats, ssr)         # Are the repeat lengths consistent?
#' (ssr <- fix_replen(nancycats, ssr)) # Nope. We need to fix them.
#' 
#' # Analyze the first population in nancycats
#' bruvo.dist(popsub(nancycats, 1), replen = ssr)
#' 
#' # View each population as a heatmap.
#' \dontrun{
#' sapply(popNames(nancycats), function(x) 
#' heatmap(as.matrix(bruvo.dist(popsub(nancycats, x), replen = ssr)), symm=TRUE))
#' }
#==============================================================================#
#' @useDynLib poppr
bruvo.dist <- function(pop, replen = 1, add = TRUE, loss = TRUE){
  # This attempts to make sure the data is true microsatellite data. It will
  # reject snp and aflp data. 
  if (pop@type != "codom" || all(is.na(unlist(lapply(alleles(pop), as.numeric))))){
    stop(non_ssr_data_warning())
  }
  # Bruvo's distance depends on the knowledge of the repeat length. If the user
  # does not provide the repeat length, it can be estimated by the smallest
  # repeat difference greater than 1. This is not a preferred method. 
  if (length(replen) != length(locNames(pop))){
    replen <- vapply(alleles(pop), function(x) guesslengths(as.numeric(x)), 1)
    warning(repeat_length_warning(replen), immediate. = TRUE)
  }
  bruvomat  <- new('bruvomat', pop, replen)
  funk_call <- match.call()
  if (length(add) != 1 || !is.logical(add) || length(loss) != 1 || !is.logical(loss)){
    stop("add and loss flags must be either TRUE or FALSE. Please check your input.")
  }
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
#' @param tree any function that can generate a tree from a distance matrix.
#'   Default is \code{\link[phangorn]{upgma}}.
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
#' @param root \code{logical} This is a parameter passed on to
#'   \code{\link[ape]{boot.phylo}}. If the \code{tree} argument produces a
#'   rooted tree (e.g. "upgma"), then this value should be \code{TRUE}. If it
#'   produces an unrooted tree (e.g. "nj"), then the value should be
#'   \code{FALSE}. By default, it is set to \code{NULL}, which will assume an
#'   unrooted phylogeny unless the function name contains "upgma".
#' 
#' @param ... any argument to be passed on to \code{\link{boot.phylo}}. eg. 
#'   \code{quiet = TRUE}.
#'   
#'   
#' @return a tree of class phylo with nodelables
#'   
#' @seealso \code{\link{bruvo.dist}}, \code{\link{nancycats}}, 
#'   \code{\link{upgma}}, \code{\link{nj}}, \code{\link{boot.phylo}}, 
#'   \code{\link{nodelabels}}, \code{\link{tab}}, 
#'   \code{\link{missingno}}.
#'   
#' @details This function will calculate a tree based off of Bruvo's distance
#'   and then utilize \code{\link[ape]{boot.phylo}} to randomly sample loci with
#'   replacement, recalculate the tree, and tally up the bootstrap support
#'   (measured in percent success). While this function can take any tree
#'   function, it has native support for two algorithms: \code{\link[ape]{nj}}
#'   and \code{\link[phangorn]{upgma}}. If you want to use any other functions,
#'   you must load the package before you use them (see examples).
#' 
#' @note \strong{Please refer to the documentation for bruvo.dist for details on
#'   the algorithm.} If the user does not provide a vector of appropriate length
#'   for \code{replen} , it will be estimated by taking the minimum difference
#'   among represented alleles at each locus. IT IS NOT RECOMMENDED TO RELY ON
#'   THIS ESTIMATION.
#'   
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
#' \dontrun{
#' 
#' # Always load the library before you specify the function.
#' library("ape")
#' 
#' # Estimate the tree based off of the BIONJ algorithm.
#' 
#' bruvo.boot(popsub(nancycats, 9), replen = ssr, tree = bionj)
#' 
#' # Utilizing  balanced FastME
#' bruvo.boot(popsub(nancycats, 9), replen = ssr, tree = fastme.bal)
#' 
#' # To change parameters for the tree, wrap it in a function.
#' # For example, let's build the tree without utilizing subtree-prune-regraft
#' 
#' myFastME <- function(x) fastme.bal(x, nni = TRUE, spr = FALSE, tbr = TRUE)
#' bruvo.boot(popsub(nancycats, 9), replen = ssr, tree = myFastME)
#' 
#' }
#' 
#==============================================================================#
#' @importFrom phangorn upgma midpoint
#' @importFrom ape nodelabels nj boot.phylo plot.phylo axisPhylo ladderize 
#' @importFrom ape add.scale.bar nodelabels tiplabels is.ultrametric
#   /     \
#   |=(o)=|
#   \     /
bruvo.boot <- function(pop, replen = 1, add = TRUE, loss = TRUE, sample = 100, 
                        tree = "upgma", showtree = TRUE, cutoff = NULL, 
                        quiet = FALSE, root = NULL, ...){
  # This attempts to make sure the data is true microsatellite data. It will
  # reject snp and aflp data. 
  if (pop@type != "codom" || all(is.na(unlist(lapply(alleles(pop), as.numeric))))){
    stop(non_ssr_data_warning())
  }
  # Bruvo's distance depends on the knowledge of the repeat length. If the user
  # does not provide the repeat length, it can be estimated by the smallest
  # repeat difference greater than 1. This is not a preferred method. 
  if (length(replen) != length(locNames(pop))){
    replen <- vapply(alleles(pop), function(x) guesslengths(as.numeric(x)), 1)
    warning(repeat_length_warning(replen), immediate. = TRUE)
  }
  bootgen <- new('bruvomat', pop, replen)
  # Steps: Create initial tree and then use boot.phylo to perform bootstrap
  # analysis, and then place the support labels on the tree.
  treechar <- paste(as.character(substitute(tree)), collapse = "")
  if ("upgma" %in% treechar){
    treefun <- upgma
  } else if ("nj" %in% treechar){
    treefun <- nj
  } else {
    treefun <- match.fun(tree)    
  }
  bootfun <- function(x){
    treefun(bruvos_distance(x, funk_call = match.call(), add = add, 
                            loss = loss))
  }

  tre <- bootfun(bootgen)
  if (is.null(root)){
    root <- ape::is.ultrametric(tre)
  }
  if (any (tre$edge.length < 0)){
    warning(negative_branch_warning(), immediate.=TRUE)
	  tre <- fix_negative_branch(tre)
  }
  if (quiet == FALSE){
    cat("\nBootstrapping...\n") 
    cat("(note: calculation of node labels can take a while even after") 
    cat(" the progress bar is full)\n\n")
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
  tre$tip.label <- indNames(pop)
  if (showtree){
    poppr.plot.phylo(tre, treechar, root)
  }
  return(tre)
}
#==============================================================================#
#
#' Create minimum spanning network of selected populations using Bruvo's 
#' distance.
#' 
#' @param gid a \code{\link{genind}} object
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
#' @param mlg.compute if the multilocus genotypes are set to "custom" (see
#'   \code{\link{mll.custom}} for details) in your genclone object, this will
#'   specify which mlg level to calculate the nodes from. See details.
#'   
#' @param palette a \code{vector} or \code{function} defining the color palette 
#'   to be used to color the populations on the graph. It defaults to 
#'   \code{\link{topo.colors}}. See examples for details.
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
#' @param threshold numeric. If greater than the default value of 0.0, this will
#'   be passed to \code{\link{mlg.filter}} prior to creating the msn.
#'
#' @param clustering.algorithm string. If \code{threshold} is greater than 0, this
#'   this will also be passed to \code{\link{mlg.filter}} prior to creating the msn.
#'   For both of these arguments, see \code{\link{mlg.filter}} for more details.
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
#'   \subsection{mlg.compute}{
#'   Each node on the graph represents a different multilocus genotype. 
#'   The edges on the graph represent genetic distances that connect the
#'   multilocus genotypes. In genclone objects, it is possible to set the
#'   multilocus genotypes to a custom definition. This creates a problem for
#'   clone correction, however, as it is very possible to define custom lineages
#'   that are not monophyletic. When clone correction is performed on these
#'   definitions, information is lost from the graph. To circumvent this, The
#'   clone correction will be done via the computed multilocus genotypes, either
#'   "original" or "contracted". This is specified in the \code{mlg.compute}
#'   argument, above.}
#'   
#' @seealso \code{\link{bruvo.dist}}, \code{\link{nancycats}}, 
#'   \code{\link{plot_poppr_msn}}, \code{\link[igraph]{minimum.spanning.tree}} 
#'   \code{\link{bruvo.boot}}, \code{\link{greycurve}} \code{\link{poppr.msn}}
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
#' vertex.label.dist=0.4, gscale=FALSE)
#' 
#' # View the whole population, but without labels.
#' bruvo.msn(nancycats, replen=rep(2, 9), vertex.label=NA)
#' }
#==============================================================================#
#' @importFrom igraph graph.adjacency plot.igraph V E minimum.spanning.tree V<- E<- print.igraph add.edges
bruvo.msn <- function (gid, replen = 1, add = TRUE, loss = TRUE, 
                       mlg.compute = "original", 
                       palette = topo.colors,
                       sublist = "All", blacklist = NULL, vertex.label = "MLG", 
                       gscale = TRUE, glim = c(0,0.8), gadj = 3, gweight = 1, 
                       wscale = TRUE, showplot = TRUE, 
                       include.ties = FALSE, threshold = 0.0, 
                       clustering.algorithm="farthest_neighbor", ...){
  if (!is.genclone(gid)){
    # convert to genclone object
    gid <- as.genclone(gid)  
  }
  # if (is.null(pop(gid)) | nPop(gid) == 1){
  #   return(singlepop_msn(gid, vertex.label, add = add, loss = loss, 
  #                        replen = replen, gscale = gscale, 
  #                        glim = glim, gadj = gadj, wscale = wscale, 
  #                        palette = palette))
  gadj <- ifelse(gweight == 1, gadj, -gadj)

  if (toupper(sublist[1]) != "ALL" | !is.null(blacklist)){
    gid <- popsub(gid, sublist, blacklist)
  }
  if (is.null(pop(gid)) | nPop(gid) == 1){
    return(singlepop_msn(gid, vertex.label, add = add, loss = loss, 
                         replen = replen, gscale = gscale, mlg.compute = mlg.compute,
                         glim = glim, gadj = gadj, wscale = wscale, 
                         palette = palette, include.ties = include.ties,
                         threshold=threshold, clustering.algorithm=clustering.algorithm, ...))
  }
  if (class(gid$mlg) != "MLG"){
    # Froce gid$mlg to be class MLG
    gid$mlg <- new("MLG", gid$mlg)
  }
  # Obtaining population information for all MLGs
  classstat <- (is.genclone(gid) | is(gid, "snpclone")) && is(gid@mlg, "MLG")
  if (classstat){
    visible <- visible(gid@mlg)
    mll(gid)  <- mlg.compute
  }
  # Updating the MLG with filtered data
  if(threshold > 0){
    filter.stats <- mlg.filter(gid,threshold,distance=bruvo.dist,algorithm=clustering.algorithm,replen=replen,stats="ALL", add = add, loss = loss)
    # TODO: The following two lines should be a product of mlg.filter
    visible(gid$mlg) <- "contracted"
    gid$mlg[] <- filter.stats[[1]]  
    # Obtaining population information for all MLGs
    cgid <- gid[if(length(-which(duplicated(gid$mlg[]))==0)) which(!duplicated(gid$mlg[])) else -which(duplicated(gid$mlg[])) ,]
    bclone <- filter.stats[[3]]
    if (!is.matrix(bclone)) bclone <- as.matrix(bclone)
  } else {
    cgid <- gid[.clonecorrector(gid), ]
    bclone <- as.matrix(bruvo.dist(cgid, replen=replen, add = add, loss = loss))
  }
  if (is.genclone(gid)){
    mlgs <- mll(gid)
    cmlg <- mll(cgid)
  } else {
    mlgs <- gid$other$mlg.vec
    cmlg <- cgid$other$mlg.vec
  }
  subs <- sort(unique(mlgs))
  mlg.cp <- mlg.crosspop(gid, mlgsub = subs, quiet=TRUE)
  if (is.numeric(mlgs)){
    names(mlg.cp) <- paste0("MLG.", sort(unique(mlgs)))    
  }

  
  # This will determine the size of the nodes based on the number of individuals
  # in the MLG. Subsetting by the MLG vector of the clone corrected set will
  # give us the numbers and the population information in the correct order.
  # Note: rank is used to correctly subset the data
  mlg.number <- table(mlgs)[rank(cmlg)]
  mlg.cp     <- mlg.cp[rank(cmlg)]
  
  ###### Create a graph #######
  g   <- graph.adjacency(as.matrix(bclone), weighted = TRUE, mode = "undirected")
  if(length(cgid@mlg[]) > 1){ 
    mst <- minimum.spanning.tree(g, algorithm = "prim", weights = E(g)$weight)
    # Add any relevant edges that were cut from the mst while still being tied for the title of optimal edge
    if(include.ties){
      tied_edges <- .Call("msn_tied_edges",as.matrix(mst[]),as.matrix(bclone),(.Machine$double.eps ^ 0.5))
      if(length(tied_edges) > 0){
        mst <- add.edges(mst, dimnames(mst[])[[1]][tied_edges[c(TRUE,TRUE,FALSE)]], weight=tied_edges[c(FALSE,FALSE,TRUE)])
      }
    }
  } else {
    mst <- minimum.spanning.tree(g)
  }

  if (!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if (toupper(vertex.label) == "MLG"){
      if (is.numeric(cmlg) && !classstat){
        vertex.label <- paste0("MLG.", cmlg)
      } else if (visible == "custom"){
        mll(gid) <- visible
        vertex.label <- correlate_custom_mlgs(gid, mlg.compute)
      } else {
        vertex.label <- paste0("MLG.", cmlg)
      }
    } else if (toupper(vertex.label) == "INDS"){
      vertex.label <- indNames(cgid)
    }
  }
  ###### Color schemes #######  
  # The pallete is determined by what the user types in the argument. It can be 
  # rainbow, topo.colors, heat.colors ...etc.
  npop   <- nPop(gid)
  pnames <- popNames(gid)
  color  <- palette_parser(palette, npop, pnames)
  
  if (length(mll(cgid)) > 1){ 
    mst <- update_edge_scales(mst, wscale, gscale, glim, gadj)
  }

  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pnames %in% names(x)])
  if (showplot){
    plot.igraph(mst, edge.width = E(mst)$width, edge.color = E(mst)$color, 
         vertex.size = mlg.number*3, vertex.shape = "pie", vertex.pie = mlg.cp, 
         vertex.pie.color = mlg.color, vertex.label = vertex.label, ...)
    graphics::legend(-1.55, 1, bty = "n", cex = 0.75, legend = pnames, 
           title = "Populations", fill = color, border = NULL)
  }
  V(mst)$size      <- mlg.number
  V(mst)$shape     <- "pie"
  V(mst)$pie       <- mlg.cp
  V(mst)$pie.color <- mlg.color
  V(mst)$label     <- vertex.label
  return(list(graph = mst, populations = pnames, colors = color))
}
#' Test repeat length consistency.
#' 
#' This function will test for consistency in the sense that all alleles are 
#' able to be represented as discrete units after division and rounding.
#' @param gid a genind object
#' @param replen a numeric vector of repeat motif lengths.
#' @return a logical vector indicating whether or not the repeat motif length is
#'   consistent.
#'   
#' @details This function is modified from the version used in 
#'   \url{http://dx.doi.org/10.5281/zenodo.13007}.
#'   
#' @references Zhian N. Kamvar, Meg M. Larsen, Alan M. Kanaskie, Everett M. 
#'   Hansen, & Niklaus J. Grünwald. Sudden_Oak_Death_in_Oregon_Forests: Spatial
#'   and temporal population dynamics of the sudden oak death epidemic in Oregon
#'   Forests. ZENODO, http://doi.org/10.5281/zenodo.13007, 2014.
#'   
#'   Kamvar, Z. N., Larsen, M. M., Kanaskie, A. M., Hansen, E. M., & Grünwald,
#'   N. J. (2015). Spatial and temporal analysis of populations of the sudden
#'   oak death pathogen in Oregon forests. Phytopathology 105:982-989.
#'   doi: \href{http://dx.doi.org/10.1094/PHYTO-12-14-0350-FI}{10.1094/PHYTO-12-14-0350-FI}
#'   
#'   Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza, and Hinrich 
#'   Schulenburg. A simple method for the calculation of microsatellite genotype
#'   distances irrespective of ploidy level. Molecular Ecology, 13(7):2101-2106,
#'   2004.
#'   
#' @export
#' @seealso \code{\link{fix_replen}} \code{\link{bruvo.dist}}
#'   \code{\link{bruvo.msn}} \code{\link{bruvo.boot}}
#' @author Zhian N. Kamvar
#' @examples
#' data(nancycats)
#' test_replen(nancycats, rep(2, 9))
test_replen <- function(gid, replen){
  alleles <- lapply(gid@all.names, as.numeric)
  are_consistent <- vapply(1:nLoc(gid), consistent_replen, logical(1), 
                           alleles, replen)
  names(are_consistent) <- locNames(gid)
  return(are_consistent)
}

consistent_replen <- function(index, alleles, replen){
  !any(duplicated(round(alleles[[index]]/replen[index])))
}

#' Find and fix inconsistent repeat lengths
#' 
#' Attempts to fix inconsistent repeat lengths found by \code{test_replen}
#' 
#' @param gid a genind object
#' @param replen a numeric vector of repeat motif lengths.
#' @param e a number to be subtracted or added to inconsistent repeat lengths to
#'   allow for proper rounding.
#' @param fix_some if \code{TRUE} (default), when there are inconsistent repeat 
#'   lengths that cannot be fixed by subtracting or adding e, those than can be 
#'   fixed will. If \code{FALSE}, the original repeat lengths will not be fixed.
#' @return a numeric vector of corrected repeat motif lengths.
#' @details This function is modified from the version used in 
#'   \url{http://dx.doi.org/10.5281/zenodo.13007}.\cr Before being fed into the
#'   algorithm to calculate Bruvo's distance, the amplicon length is divided by
#'   the repeat unit length. Because of the amplified primer sequence attached
#'   to sequence repeat, this division does not always result in an integer and
#'   so the resulting numbers are rounded. The rounding also protects against
#'   slight mis-calls of alleles. Because we know that \deqn{\frac{(A - e) - (B
#'   - e)}{r}}{((A - e) - (B - e))/r} is equivalent to \deqn{\frac{A - B}{r}}{(A
#'   - B)/r}, we know that the primer sequence will not alter the relationships
#'   between the alleles. Unfortunately for nucleotide repeats that have powers
#'   of 2, rounding in R is based off of the IEC 60559 standard (see
#'   \code{\link{round}}), that means that any number ending in 5 is rounded to
#'   the nearest \emph{even} digit. This function will attempt to alleviate this
#'   problem by adding a very small amount to the repeat length so that division
#'   will not result in a 0.5. If this fails, the same amount will be
#'   subtracted. If neither of these work, a warning will be issued and it is up
#'   to the user to determine if the fault is in the allele calls or the repeat
#'   lengths.
#' @export
#' 
#' @references Zhian N. Kamvar, Meg M. Larsen, Alan M. Kanaskie, Everett M. 
#'   Hansen, & Niklaus J. Grünwald. Sudden_Oak_Death_in_Oregon_Forests: Spatial
#'   and temporal population dynamics of the sudden oak death epidemic in Oregon
#'   Forests. ZENODO, http://doi.org/10.5281/zenodo.13007, 2014.
#'   
#'   Kamvar, Z. N., Larsen, M. M., Kanaskie, A. M., Hansen, E. M., & Grünwald,
#'   N. J. (2015). Spatial and temporal analysis of populations of the sudden
#'   oak death pathogen in Oregon forests. Phytopathology 105:982-989.
#'   doi: \href{http://dx.doi.org/10.1094/PHYTO-12-14-0350-FI}{10.1094/PHYTO-12-14-0350-FI}
#'   
#'   Ruzica Bruvo, Nicolaas K. Michiels, Thomas G. D'Souza, and Hinrich 
#'   Schulenburg. A simple method for the calculation of microsatellite genotype
#'   distances irrespective of ploidy level. Molecular Ecology, 13(7):2101-2106,
#'   2004.
#'   
#' @author Zhian N. Kamvar
#' @seealso \code{\link{test_replen}} \code{\link{bruvo.dist}}
#'   \code{\link{bruvo.msn}} \code{\link{bruvo.boot}}
#' @examples
#' 
#' data(nancycats)
#' fix_replen(nancycats, rep(2, 9))
#' # Let's start with an example of a tetranucleotide repeat motif and imagine
#' # that there are twenty alleles all 1 step apart:
#' (x <- 1:20L * 4L)
#' # These are the true lengths of the different alleles. Now, let's add the
#' # primer sequence to them. 
#' (PxP <- x + 21 + 21)
#' # Now we make sure that x / 4 is equal to 1:20, which we know each have
#' # 1 difference.
#' x/4
#' # Now, we divide the sequence with the primers by 4 and see what happens.
#' (PxPc <- PxP/4)
#' (PxPcr <- round(PxPc))
#' diff(PxPcr) # we expect all 1s
#' 
#' # Let's try that again by subtracting a tiny amount from 4
#' (PxPc <- PxP/(4 - 1e-5))
#' (PxPcr <- round(PxPc))
#' diff(PxPcr)
fix_replen <- function(gid, replen, e = 1e-5, fix_some = TRUE){
  if (length(replen) != nLoc(gid)) {
    stop(paste0("length of repeats (", length(replen), ") does not equal",
                " the number of loci (", nLoc(gid), ")."))
  }
  consistent_reps <- test_replen(gid, replen)
  names(replen)   <- locNames(gid)
  ADD <- FALSE
  SUB <- FALSE
  newReps <- replen
  while (any(!consistent_reps)){
    if (!SUB){
      newReps[!consistent_reps] <- newReps[!consistent_reps] - e
      SUB <- TRUE
    } else {
      newReps[!consistent_reps] <- newReps[!consistent_reps] + (2*e)
      ADD <- TRUE
    }
    consistent_reps <- test_replen(gid, newReps)
    if (any(!consistent_reps) & ADD & SUB){
      inconsistent <- paste(names(replen[!consistent_reps]), collapse = ", ")
      msg <- paste("The repeat lengths for", inconsistent, 
                   "are not consistent.\n\n",
                   "This might be due to inconsistent allele calls or repeat",
                   "lengths that are too large.\n",
                   "Check the alleles to make sure there are no duplicated",
                   "or similar alleles that might end up being the same after",
                   "division.\n")
      if (fix_some){
        original <- test_replen(gid, replen)
        fixed    <- paste(names(replen[!original & consistent_reps]), collapse = ", ")
        msg <- paste(msg, "\nRepeat lengths with some modification are",
                     "being returned:", fixed)
        newReps[!consistent_reps] <- replen[!consistent_reps]
      } else {
        msg <- paste(msg, "\nOriginal repeat lengths are being returned.")
        newReps <- replen
      }
      warning(msg, immediate. = TRUE)
      consistent_reps <- TRUE
    }
  }
  return(newReps)
}
