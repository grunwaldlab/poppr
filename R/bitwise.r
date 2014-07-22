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
#' Calculate a distance matrix comparing samples based on the number of alleles
#' that differ in zygosity.
#' 
#' This function does pairwise comparisons between diploid samples in a genlight
#' object. The number representing the distance between two samples is equal to
#' the number of alleles in the samples that do not have the same zygosity.
#'
#' @name bitwise
#'
#' @param x a genlight object. 
#'
#' @param percent \code{logical}. Should the distance be represented from 0 to 1? 
#' Default set to \code{TRUE}. \code{FALSE} will return the distance represented 
#' as integers from 1 to n where n is the number of loci. 
#'  
#' @param mat \code{logical}. Return a matrix object. Default set to  
#' \code{FALSE}, returning a dist object. \code{TRUE} returns a matrix object. 
#'
#' @param missing_match \code{logical}. Determines whether two samples differing
#' by missing data in a location should be counted as matching at that location.
#' Default set to \code{TRUE}, which forces missing data to match with anything. 
#' \code{FALSE} forces missing data to not match with any other information. 
#'
#' @param threads The maximum number of parallel threads to be used within this
#'   function. A value of 0 (default) will attempt to use as many threads as there
#'   are available cores/CPUs. In most cases this is ideal. A value of 1 will force
#'   the function to run serially, which may increase stability on some systems.
#'   Other values may be specified, but should be used with caution.
#'
#' 
#' @return Pairwise distances between individuals present in the genlight object.
#' @author Zhian N. Kamvar, Jonah Brooks
#' 
#' @export
#==============================================================================#
bitwise.dist <- function(x, percent=TRUE, mat=FALSE, missing_match=TRUE, threads=0){
  stopifnot(class(x)[1] == "genlight")
  # Stop if the ploidy of the genlight object is not consistent
  stopifnot(min(ploidy(x)) == max(ploidy(x))) 
  # Stop if the ploidy of the genlight object is not haploid or diploid
  stopifnot(min(ploidy(x)) == 2 || min(ploidy(x)) == 1)
  ploid     <- min(ploidy(x))
  ind.names <- indNames(x)
  inds      <- nInd(x)
  numPairs   <- nLoc(x)

  # Threads must be something that can cast to integer
  if(!is.numeric(threads) && !is.integer(threads) && threads >= 0)
  {
    stop("Threads must be a non-negative numeric or integer value")
  }

  # Cast parameters to proper types before passing them to C
  threads <- as.integer(threads)

  if (ploid == 1)
  {
    pairwise_dist <- .Call("bitwise_distance_haploid", x, missing_match, threads)
  }
  else
  {
    pairwise_dist <- .Call("bitwise_distance_diploid", x, missing_match, threads)
  }
  dist.mat <- pairwise_dist
  dim(dist.mat) <- c(inds,inds)
  colnames(dist.mat)            <- ind.names
  rownames(dist.mat)            <- ind.names
  if (percent){
    dist.mat <- dist.mat/(numPairs)
  }
  if (mat == FALSE){
    dist.mat <- as.dist(dist.mat)
  }
  return(dist.mat)
}



#==============================================================================#
#' Calculates and returns a vector of Pgen values for the given genlight object.
#' Each element represents the probability that the individual at that element
#' would have been produced via random mating using estimates derived from the
#' genotypes present in the genlight object.
#'
#' @rdname bitwise
#'
#' @param x a genlight object. 
#' 
#' @param percent a \code{logical} to determine whether the values should be returned
#'  as percentages or logarithms of percentages. \code{FALSE} is the default, and 
#'  returns the logarithmic values rather than the percentage values. This option has
#'  a much larger range and is highly recommended. \code{TRUE} returns the percentage
#'  chance for each genotype to be produced via random mating, rather than the log
#'  equivalent.
#' 
#' @return A vector containing one Pgen value for each genotype in the genlight object.
#' @author Zhian N. Kamvar, Jonah Brooks
#' 
#' @export
#==============================================================================#
bitwise.pgen <- function(x, percent) {
  stopifnot(class(x)[1] == "genlight")

  pgen_vector <- .Call("get_pgen_vector",x)

  if(percent)
  {
    pgen_vector <- exp(pgen_vector)
  }

  return(pgen_vector);

}
