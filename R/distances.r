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
#' Calculate a distance matrix based on relative dissimilarity
#' 
#' diss.dist uses the same discreet dissimilarity matrix utilized by the index
#' of association (see \code{\link{ia}} for details). It returns a distance
#' reflecting a ratio of the number of observed differences by the number of
#' possible differences. Eg. two individuals who share half of the same alleles
#' will have a distance of 0.5. This function can analyze distances for any
#' marker system.
#' 
#' @param x a \code{\link{genind}} object. 
#' 
#' @param diff \code{logical}. If \code{TRUE}, this will count the number of
#' differences between individuals. \code{FALSE} will count the number of 
#' similarities. Default set to \code{TRUE}
#' 
#' @param frac \code{logical}. Should the distance be represented from 0 to 1?
#' Default set to \code{TRUE}. \code{FALSE} will return the distance represented
#' as integers from 1 to n where n is the number of loci*ploidy.
#' 
#' @param mat \code{logical}. Return a matrix object. Default set to 
#' \code{FALSE}, returning a dist object. \code{TRUE} returns a matrix object.
#' 
#' @return Pairwise distances between individuals present in the genind object.
#' @author Zhian N. Kamvar
#' 
#' @note This is exactly the same as \code{\link{provesti.dist}}, except that
#' it performs better for large numbers of individuals (n > 125). 
#'
#' @examples
#' 
#' # A simple example. Let's analyze the mean distance among populations of A.
#' # euteiches.
#' 
#' data(Aeut)
#' mean(diss.dist(popsub(Aeut, 1)))
#' \dontrun{
#' mean(diss.dist(popsub(Aeut, 2)))
#' mean(diss.dist(Aeut))
#' }
#' @export
#==============================================================================#

diss.dist <- function(x, diff=TRUE, frac=TRUE, mat=FALSE){
  stopifnot(is.genind(x))
  ploid     <- ploidy(x)
  ind.names <- x@ind.names
  inds      <- nrow(x@tab)
  np        <- choose(inds, 2)
  dist.mat  <- matrix(data = 0, nrow = inds, ncol = inds)
  numLoci   <- length(x@loc.names)
  type      <- x@type
  if(type == "PA"){
    dist_by_locus <- matrix(.Call("pairdiffs", x@tab))
    ploid <- 1
  } else {
    x <- seploc(x)
    dist_by_locus <- vapply(x, function(x) .Call("pairdiffs", x@tab)*(ploid/2),
                            numeric(np))
  }
#   if (!alleles & type == "codom"){
#     dist_by_locus[dist_by_locus > 0] <- 1
#     ploid <- 1
#   }
  if (!diff){
    max_dist      <- numLoci*ploid
    dist_by_locus <- max_dist - dist_by_locus
  }
  dist.mat[lower.tri(dist.mat)] <- rowSums(dist_by_locus)
  colnames(dist.mat)            <- ind.names
  rownames(dist.mat)            <- ind.names
  if (frac){
    dist.mat <- dist.mat/(numLoci*ploid)
  }
  dist.mat <- as.dist(dist.mat)
  if (mat == TRUE){
    dist.mat <- as.matrix(dist.mat)
  }
  return(dist.mat)
}


#==============================================================================#
# Calculating Nei's distance for a genind object. Lifted from dist.genpop with
# modifications.
#' Calculate Genetic Distance for a genind or genclone object.
#' 
#' These functions are modified from the function \link[adegenet]{dist.genpop} to
#' be applicable for distances between individuals.
#' 
#' @param x a \linkS4class{genind}, \linkS4class{genclone}, or matrix object.
#'   
#' @param warning If \code{TRUE}, a warning will be printed if any infinite 
#'   values are detected and replaced. If \code{FALSE}, these values will be 
#'   replaced without warning. See Details below.
#'   
#' @return an object of class dist with the same number of observations as the 
#'   number of individuals in your data.
#'   
#' @details It is important to be careful with the interpretation of these 
#'   distances as they were originally intended for calculation of 
#'   between-population distance. As Nei's distance is the negative log of 0:1, 
#'   this means that it is very possible to obtain distances of infinity. When 
#'   this happens, infinite values are corrected to be 10 * max(D) where D is
#'   the distance matrix without infinite values.
#'   
#' @note Provesti's distance is identical to \code{\link{diss.dist}}, except
#'   that \code{\link{diss.dist}} is optimized for a larger number of
#'   individuals (n > 125) at the cost of required memory.
#' 
#' @rdname genetic_distance
#' @export
#' @examples
#' 
#' data(nancycats)
#' nan9 <- popsub(nancycats, 9)
#' neinan <- nei.dist(nan9)
#' ednan <- edwards.dist(nan9)
#' rodnan <- rogers.dist(nan9)
#' reynan <- reynolds.dist(nan9)
#' pronan <- provesti.dist(nan9)
#' 
#==============================================================================#
nei.dist <- function(x, warning = TRUE){
  if (is(x, "gen"))
    MAT    <- get_gen_mat(x)
  else if (length(dim(x)) == 2)
    MAT <- x
  else
    stop("Object must be a matrix or genind object")
  IDMAT <- MAT %*% t(MAT)
  vec   <- sqrt(diag(IDMAT))
  IDMAT <- IDMAT/vec[col(IDMAT)]
  IDMAT <- IDMAT/vec[row(IDMAT)]
  D     <- -log(IDMAT)
  if (any(D == Inf)){
    D <- infinite_vals_replacement(D, warning)
  }
  D     <- as.dist(D)
  labs  <- get_gen_dist_labs(x)
  D     <- make_attributes(D, length(labs), labs, "Nei", match.call())
  return(D)
}

# modified from adegenet dist.genpop

#' @rdname genetic_distance
#' @export
edwards.dist <- function(x){
  if (is(x, "gen")){ 
    MAT  <- get_gen_mat(x)
    nloc <- length(x@loc.names)
  } else if (length(dim(x)) == 2){
    MAT  <- x
    nloc <- ncol(x)
  } else{
    stop("Object must be a matrix or genind object")
  }
  MAT     <- sqrt(MAT)
  D       <- MAT%*%t(MAT)
  D       <- 1 - D/nloc # Negative number are generated here.
  diag(D) <- 0
  D       <- sqrt(D)
  D       <- as.dist(D)
  labs    <- get_gen_dist_labs(x)
  D       <- make_attributes(D, length(labs), labs, "Edwards", match.call())
  return(D)
}


#' @rdname genetic_distance
#' @export
rogers.dist <- function(x){
  if (is(x, "gen")){ 
    if (is.genind(x) && x@type == "PA"){
      MAT     <- x@tab
      nloc    <- length(x@loc.names)
      loc.fac <- factor(x@loc.names, levels = x@loc.names)
      nlig    <- nrow(x@tab)
    } else {
      MAT     <- get_gen_mat(x)
      nloc    <- length(x@loc.names)
      loc.fac <- x@loc.fac
      nlig    <- nrow(x@tab)      
    }
  }
  else if (length(dim(x)) == 2){
    MAT     <- x
    nloc    <- ncol(x)
    loc.fac <- factor(colnames(x), levels = colnames(x))
    nlig    <- nrow(x)
  }
  else{
    stop("Object must be a matrix or genind object")
  }
  # kX is a list of K=nloc matrices
  kX <- lapply(split(MAT, loc.fac[col(MAT)]), matrix, nrow = nlig)
  D  <- matrix(0, nlig, nlig)
  for(i in 1:length(kX)){
    D <- D + dcano(kX[[i]])
  }
  D    <- D/length(kX)
  D    <- as.dist(D)
  labs <- get_gen_dist_labs(x)
  D    <- make_attributes(D, nlig, labs, "Rogers", match.call())
  return(D)
}

#' @rdname genetic_distance
#' @export
reynolds.dist <- function(x){
  if (is(x, "gen")){ 
    MAT    <- get_gen_mat(x)
    nloc   <- length(x@loc.names)
  }
  else if (length(dim(x)) == 2){
    MAT  <- x
    nloc <- ncol(x)
  }
  else{
    stop("Object must be a matrix or genind object")
  }
  denomi       <- MAT %*% t(MAT)
  vec          <- apply(MAT, 1, function(x) sum(x*x))
  D            <- -2*denomi + vec[col(denomi)] + vec[row(denomi)]
  diag(D)      <- 0
  denomi       <- 2*nloc - 2*denomi
  diag(denomi) <- 1
  D            <- D/denomi
  D            <- sqrt(D)
  D            <- as.dist(D)
  labs         <- get_gen_dist_labs(x)
  D            <- make_attributes(D, length(labs), labs, "Reynolds", 
                                  match.call())
  return(D)
}

#' @rdname genetic_distance
#' @export
provesti.dist <- function(x){
  if (is(x, "gen")){
    MAT   <- get_gen_mat(x)
    nlig  <- nrow(x@tab)
    nloc  <- length(x@loc.names)
    ploid <- x@ploidy
  } else if (length(dim(x)) == 2){
    MAT  <- x
    nlig <- nrow(x)
    nloc <- ncol(x)
  } else {
    stop("Object must be a matrix or genind object")
  }
  w0   <- 1:(nlig-1)
  loca <- function(k, nlig, MAT, nLoc){
    w1     <- (k+1):nlig
    resloc <- vapply(w1, function(y) sum(abs(MAT[k, ] - MAT[y, ]), na.rm = TRUE), numeric(1))
    if (x@type == "codom"){
      # This only applies to codominant data because dominant data can only take
      # on a single state. Dividing by two indicates that the observations can
      # occupy co-occurring states.
      resloc <- resloc/2
    }
    return(resloc/nloc)
  }

  d    <- unlist(lapply(w0, loca, nlig, MAT, nLoc))

  labs <- get_gen_dist_labs(x)
  d    <- make_attributes(d, nlig, labs, "Provesti", match.call())

  # resmat <- matrix(numeric(0), nlig, nlig)
  # resmat[lower.tri(resmat)] <- d
  # d <- as.dist(resmat)
  # attr(d, "Labels") <- labs
  return(d)
}


#==============================================================================#
#' Calculate a dendrogram with bootstrap support using any distance applicable
#' to genind or genclone objects.
#' 
#' @param x a \linkS4class{genind}, \linkS4class{genclone}, or matrix object.
#'   
#' @param tree one of "upgma" (Default) or "nj" defining the type of dendrogram 
#'   to be produced, UPGMA or Neighbor-Joining.
#'   
#' @param distance a character or function defining the distance to be applied 
#'   to x. Defaults to \code{\link{nei.dist}}.
#'   
#' @param sample An integer representing the number of bootstrap replicates
#'   Default is 100.
#'   
#' @param cutoff An integer from 0 to 100 setting the cutoff value to return the
#'   bootstrap values on the nodes. Default is 0.
#'   
#' @param showtree If \code{TRUE} (Default), a dendrogram will be plotted. If
#'   \code{FALSE}, nothing will be plotted.
#'   
#' @param missing any method to be used by \code{\link{missingno}}: "mean"
#'   (default), "zero", "loci", "genotype", or "ignore".
#'
#' @param quiet if \code{FALSE} (Default), a progress bar will be printed to 
#' screen. 
#'   
#' @param ... any parameters to be passed off to the distance method.
#'   
#' @return an object of class \code{\link[ape]{phylo}}.
#'   
#' @details This function utilizes an internal class called
#' \code{\linkS4class{bootgen}} that allows bootstrapping of objects that
#' inherit the genind class. This is necessary due to the fact that columns in
#' the genind matrix are defined as alleles and are thus interrelated. This
#' function will specifically bootstrap loci so that results are biologically
#' relevant. With this function, the user can also define a custom distance to
#' be performed on the genind or genclone object.
#' 
#' @note \code{\link{provesti.dist}} and \code{\link{diss.dist}} are exactly the
#'   same, but \code{\link{diss.dist}} scales better for large numbers of
#'   individuals (n > 125) at the cost of required memory.
#' 
#' @seealso \code{\link{nei.dist}} \code{\link{edwards.dist}}
#'   \code{\link{rogers.dist}} \code{\link{reynolds.dist}}
#'   \code{\link{provesti.dist}} \code{\link{diss.dist}}
#'   \code{\link{bruvo.boot}} \code{\link[ape]{boot.phylo}}
#'   \code{\link[adegenet]{dist.genpop}} \code{\link{dist}}
#'   
#' @export
#' @examples
#' 
#' data(nancycats)
#' nan9 <- popsub(nancycats, 9)
#' 
#' set.seed(9999)
#' # Generate a tree using nei's distance
#' neinan <- aboot(nan9, dist = nei.dist)
#' 
#' set.seed(9999)
#' # Generate a tree using custom distance
#' bindist <- function(x) dist(x$tab, method = "binary")
#' binnan <- aboot(nan9, dist = bindist)
#' 
#' \dontrun{
#' # AFLP data
#' data(Aeut)
#' 
#' # Nei's distance
#' anei <- aboot(Aeut, dist = nei.dist, sample = 1000, cutoff = 50)
#' 
#' # Rogers' distance
#' arog <- aboot(Aeut, dist = roger.dist, sample = 1000, cutoff = 50)
#' 
#' }
#==============================================================================#
aboot <- function(x, tree = "upgma", distance = "nei.dist", sample = 100,
                     cutoff = 0, showtree = TRUE, missing = "mean", quiet = FALSE,
                     ...){
  if (is.genind(x)){
    x <- missingno(x, missing)
  }
  if (x@type == "PA"){
    xboot           <- x@tab
    colnames(xboot) <- locNames(x)
    if (is.genpop(x)){
      rownames(xboot) <- x@pop.names
    } else {
      rownames(xboot) <- indNames(x)
    }
  } else {
    xboot <- new("bootgen", x)
  }
  ARGS     <- c("nj", "upgma")
  treearg  <- match.arg(tree, ARGS)
  treefunk <- tree_generator(treearg, distance, ...)
  xtree    <- treefunk(xboot)
  if (any(xtree$edge.len < 0)){
    xtree <- fix_negative_branch(xtree)
  }
  root     <- ifelse(treearg == "nj", FALSE, TRUE)
  nodelabs <- boot.phylo(xtree, xboot, treefunk, B = sample, rooted = root, quiet = quiet)
  nodelabs <- (nodelabs/sample)*100
  nodelabs <- ifelse(nodelabs >= cutoff, nodelabs, NA)
  if (is.genind(x)){
    xtree$tip.label <- indNames(x)
  } else {
    xtree$tip.label <- x@pop.names
  }
  xtree$node.label <- nodelabs
  if (showtree){
    poppr.plot.phylo(xtree, tree)
  }
  return(xtree)
}