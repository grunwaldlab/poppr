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
#' Display a greyscale gradient adjusted to specific parameters
#' 
#' This function has one purpose. It is for deciding the appropriate scaling for
#' a grey palette to be used for edge weights of a minimum spanning network.
#' 
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#' the upper and lower limits for the \code{\link{gray}} function. Default is 0
#' (black) and 0.8 (20\% black). 
#' 
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that
#' will serve as the exponent to the edge weight to scale the grey value to
#' represent that weight.
#' 
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale
#' will be weighted to emphasize the differences between closely related nodes.
#' If it is 2, the grey scale will be weighted to emphasize the differences
#' between more distantly related nodes. 
#' 
#' @return A plot displaying a grey gradient from 0.001 to 1 with minimum and 
#' maximum values displayed as yellow lines, and an equation for the correction 
#' displayed in red. 
#'
#' @examples
#' # Normal grey curve with an adjustment of 3, an upper limit of 0.8, and
#' # weighted towards smaller values.
#' greycurve()
#' 
#' # 1:1 relationship grey curve.
#' greycurve(gadj=1, glim=1:0)
#' 
#' # Grey curve weighted towards larger values.
#' greycurve(gweight=2)
#' 
#' # Same as the first, but the limit is 1.
#' greycurve(glim=1:0)
#' 
#' # Setting the lower limit to 0.1 and weighting towards larger values.
#' greycurve(glim=c(0.1,0.8), gweight=2)
#'
#' @export
#==============================================================================#
greycurve <- function(glim = c(0,0.8), gadj = 3, gweight = 1){
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  adjustcurve(seq(0.001, 1, 0.001), glim, correction=gadj, show=TRUE)
}

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
#' @param pop a \code{\link{genind}} object. 
#' 
#' @return A distance object.
#'
#' @examples
#' 
#' # A simple example. Let's analyze the mean distance among populations of A.
#' # euteiches.
#' 
#' data(Aeut)
#' mean(diss.dist(popsub(Aeut, 1)))
#' mean(diss.dist(popsub(Aeut, 2)))
#' mean(diss.dist(Aeut))
#'
#' @export
#==============================================================================#

diss.dist <- function(pop){
  ploid <- ploidy(pop)
  ind.names <- pop@ind.names
  inds <- nInd(pop)
  np <- choose(inds, 2)
  dist.vec <- matrix(data = 0, nrow=inds, ncol=inds)
  if(pop@type == "PA"){
    dist.vec[lower.tri(dist.vec)] <- .Call("pairdiffs",pop@tab*ploid)/ploid
  }
  else{
    pop <- seploc(pop)
    numLoci <- length(pop)
    temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
    temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs",x@tab*ploid)/ploid, 
                            temp.d.vector[, 1])
    dist.vec[lower.tri(dist.vec)] <- rowSums(temp.d.vector)
  }
  colnames(dist.vec) <- ind.names
  rownames(dist.vec) <- ind.names
  loci <- ifelse(is.list(pop), length(pop), nLoc(pop))
  return(as.dist(dist.vec/(loci*ploid)))
}

