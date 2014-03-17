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
  inds      <- nInd(x)
  np        <- choose(inds, 2)
  dist.mat  <- matrix(data = 0, nrow = inds, ncol = inds)
  numLoci   <- nLoc(x)
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


