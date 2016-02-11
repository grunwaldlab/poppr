#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#
# This software was authored by Zhian N. Kamvar and Javier F. Tabima, graduate 
# students at Oregon State University; Jonah C. Brooks, undergraduate student at
# Oregon State University; and Dr. Nik Grünwald, an employee of USDA-ARS.
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
#' Calculate a dissimilarity distance matrix for SNP data.
#' 
#' This function performs the same task as \code{\link{diss.dist}}, calculating
#' the number of allelic differences between two samples.
#' 
#' @param x a genlight, genind, genclone, or snpclone object.
#'   
#' @param percent \code{logical}. Should the distance be represented from 0 to
#'   1? Default set to \code{TRUE}. \code{FALSE} will return the distance
#'   represented as integers from 1 to n where n is the number of loci.
#'   
#' @param mat \code{logical}. Return a matrix object. Default set to 
#'   \code{FALSE}, returning a dist object. \code{TRUE} returns a matrix object.
#'   
#' @param missing_match \code{logical}. Determines whether two samples differing
#'   by missing data in a location should be counted as matching at that
#'   location. Default set to \code{TRUE}, which forces missing data to match
#'   with anything. \code{FALSE} forces missing data to not match with any other
#'   information.
#'   
#' @param differences_only \code{logical}. Determines whether the matrix should 
#'   count differences or distances. For instance, 0 to 2 would be a distance of
#'   2 but a difference of 1.
#'   
#' @param threads The maximum number of parallel threads to be used within this 
#'   function. A value of 0 (default) will attempt to use as many threads as
#'   there are available cores/CPUs. In most cases this is ideal. A value of 1
#'   will force the function to run serially, which may increase stability on
#'   some systems. Other values may be specified, but should be used with
#'   caution.
#'
#' 
#' @details The distance calculated here is quite simple and goes by many names,
#'   depending on its application. The most familiar name might be the Hamming
#'   distance, or the number of differences between two strings.
#'   
#' @return A dist object containing pairwise distances between samples.
#'   
#' @author Zhian N. Kamvar, Jonah C. Brooks
#' 
#' @export
#' @seealso \code{\link{diss.dist}},
#'    \code{\link{snpclone}},
#'    \code{\link[adegenet]{genlight}},
#'    \code{\link{win.ia}}, 
#'    \code{\link{samp.ia}}
#' @examples
#' set.seed(999)
#' x <- glSim(n.ind = 10, n.snp.nonstruc = 5e2, n.snp.struc = 5e2, ploidy = 2)
#' x
#' system.time(xd <- bitwise.dist(x))
#' xd
#==============================================================================#
bitwise.dist <- function(x, percent = TRUE, mat = FALSE, missing_match = TRUE, 
                         differences_only = FALSE, threads = 0){
  stopifnot(inherits(x, c("genlight", "genclone", "genind", "snpclone")))
  # Stop if the ploidy of the genlight object is not consistent
  stopifnot(min(ploidy(x)) == max(ploidy(x))) 
  # Stop if the ploidy of the genlight object is not haploid or diploid
  stopifnot(min(ploidy(x)) == 2 || min(ploidy(x)) == 1)

  ploid     <- min(ploidy(x))
  ind.names <- indNames(x)
  inds      <- nInd(x)
  numPairs  <- nLoc(x)

  # Use Prevosti if this is a genclone or genind object
  if(!is(x, "genlight")){
    dist.mat <- prevosti.dist(x)
    if (percent == FALSE){
      dist.mat <- dist.mat*ploid*numPairs
    }
    if (mat == TRUE){
      dist.mat <- as.matrix(dist.mat)
    }
    # Return this matrix and exit function
    return(dist.mat)
  }

  # Continue function for genlight objects

  # Ensure that every SNPbin object has data for all chromosomes
  if(ploid == 2){
    for(i in 1:length(x$gen)){
      if(length(x$gen[[i]]$snp) == 1){
        x$gen[[i]]$snp <- append(x$gen[[i]]$snp, list(as.raw(rep(0,length(x$gen[[i]]$snp[[1]])))))
      }
    }
  }
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
    pairwise_dist <- .Call("bitwise_distance_diploid", x, missing_match, differences_only, threads)
  }
  dist.mat <- pairwise_dist
  dim(dist.mat) <- c(inds,inds)
  colnames(dist.mat) <- ind.names
  rownames(dist.mat) <- ind.names
  if (percent){
    if(differences_only)
    {
      dist.mat <- dist.mat/(numPairs)
    }
    else
    {
      dist.mat <- dist.mat/(numPairs*ploid)
    }
  }
  if (mat == FALSE){
    dist.mat <- as.dist(dist.mat)
  }
  return(dist.mat)
}

#==============================================================================#
#' Determines whether openMP is support on this system.
#'
#' @return FALSE if openMP is not supported, TRUE if it is
#' @author Zhian N. Kamvar, Jonah C. Brooks
#' 
#' @export
#' @examples
#' poppr_has_parallel()
#==============================================================================#
poppr_has_parallel <- function(){

  supported <- .Call("omp_test", PACKAGE = "poppr")

  if(supported == 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }

}



#==============================================================================#
#' Calculate the index of association between samples in a genlight object.
#' 
#' This function parses over a genlight object to calculate and return the
#' index of association for those samples.
#'
#' @param x a genlight object. 
#'
#' @param missing_match a boolean determining whether missing data should be
#'   considered a match. If TRUE (default) missing data at a locus will match
#'   with any data at that locus in each comparison. If FALSE, missing data at
#'   a locus will cause all comparisons to return the maximum possible distance
#'   at that locus (ie, if sample 1 has missing data at locus 1, and sample 2 
#'   is heterozygous at locus 1, the distance at that locus will be 1. If sample
#'   2 was heterozygous or missing at locus 1, the distance would be 2.
#'
#' @param differences_only a boolean determining how distance should be counted
#'   for diploids. Whether TRUE or FALSE the distance between a heterozygous locus
#'   and a homozygous locus is 1. If FALSE (default) the distance between opposite
#'   homozygous loci is 2. If TRUE that distance counts as 1, indicating only that
#'   the two samples differ at that locus.
#'
#' @param threads The maximum number of parallel threads to be used within this
#'   function. A value of 0 (default) will attempt to use as many threads as there
#'   are available cores/CPUs. In most cases this is ideal. A value of 1 will force
#'   the function to run serially, which may increase stability on some systems.
#'   Other values may be specified, but should be used with caution.
#'
#' @return Index of association representing the samples in this genlight object.
#' @author Zhian N. Kamvar, Jonah C. Brooks
#' 
#' @export
#' @seealso \code{\link{win.ia}}, \code{\link{samp.ia}}
#' @keywords internal
#==============================================================================#
bitwise.ia <- function(x, missing_match=TRUE, differences_only=FALSE, threads=0){
  stopifnot(class(x)[1] %in% c("genlight", "snpclone"))
  # Stop if the ploidy of the genlight object is not consistent
  stopifnot(min(ploidy(x)) == max(ploidy(x))) 
  # Stop if the ploidy of the genlight object is not haploid or diploid
  stopifnot(min(ploidy(x)) == 2 || min(ploidy(x)) == 1)

  ploid     <- min(ploidy(x))

  # Threads must be something that can cast to integer
  if(!is.numeric(threads) && !is.integer(threads) && threads >= 0)
  {
    stop("Threads must be a non-negative numeric or integer value")
  }
  # Cast parameters to proper types before passing them to C
  threads <- as.integer(threads)
  # Ensure that every SNPbin object has data for all chromosomes
  if(ploid == 2){
    for(i in 1:length(x$gen)){
      if(length(x$gen[[i]]$snp) == 1){
        x$gen[[i]]$snp <- append(x$gen[[i]]$snp, list(as.raw(rep(0,length(x$gen[[i]]$snp[[1]])))))
      }
    }
    IA <- .Call("association_index_diploid", x, missing_match, differences_only, 
                threads, PACKAGE = "poppr")
  }
  else if(ploid == 1)
  {
    IA <- .Call("association_index_haploid", x, missing_match, threads, 
                PACKAGE = "poppr")
  }
  else
  {
    stop("bitwise.ia only supports haploids and diploids")
  }
  #TODO: Allow for automated index generation, such as random or window based

  #TODO: Call C function and return
  

  return(IA)

}

#==============================================================================#
#' Calculate windows of the index of association for genlight objects.
#' 
#' Genlight objects can contain millions of loci. Since it does not make much 
#' sense to calculate the index of association over that many loci, this
#' function will scan windows across the loci positions and calculate the index
#' of association.
#' 
#' @param x a genlight object.
#'   
#' @param window an integer specifying the size of the window.
#'   
#' @param min.snps an integer specifying the minimum number of snps allowed per 
#'   window. If a window does not meet this criteria, the value will return as
#'   NA.
#'   
#' @param threads The maximum number of parallel threads to be used within this 
#'   function. A value of 0 (default) will attempt to use as many threads as
#'   there are available cores/CPUs. In most cases this is ideal. A value of 1
#'   will force the function to run serially, which may increase stability on
#'   some systems. Other values may be specified, but should be used with
#'   caution.
#'   
#' @param quiet if \code{FALSE}, a progress bar will be printed to the screen.
#'   
#' @return Index of association representing the samples in this genlight
#'   object.
#'   
#' @note this will calculate the standardized index of association from Agapow
#' 2001. See \code{\link{ia}} for details.
#' 
#' @author Zhian N. Kamvar, Jonah C. Brooks
#'   
#' @export
#' @seealso \code{\link[adegenet]{genlight}},
#'    \code{\link{snpclone}},
#'    \code{\link{samp.ia}},
#'    \code{\link{ia}},
#'    \code{\link{bitwise.dist}}
#' @examples
#' 
#' # with structured snps assuming 1e4 positions
#' set.seed(999)
#' x <- glSim(n.ind = 10, n.snp.nonstruc = 5e2, n.snp.struc = 5e2, ploidy = 2)
#' position(x) <- sort(sample(1e4, 1e3))
#' res <- win.ia(x, window = 300L) # Calculate for windows of size 300
#' plot(res, type = "l")
#' 
#' \dontrun{
#' # unstructured snps
#' set.seed(999)
#' x <- glSim(n.ind = 10, n.snp.nonstruc = 1e3, ploidy = 2)
#' position(x) <- sort(sample(1e4, 1e3))
#' res <- win.ia(x, window = 300L) # Calculate for windows of size 300
#' plot(res, type = "l")
#' }
#' 
#==============================================================================#
win.ia <- function(x, window = 100L, min.snps = 3L, threads = 1L, quiet = FALSE){
  stopifnot(is(x, "genlight"))
  if (!is.null(position(x))){
    xpos <- position(x)
  } else {
    xpos <- seq(nLoc(x))
  }
  winmat  <- make_windows(maxp = max(xpos), minp = min(xpos), window = window)
  nwin    <- nrow(winmat)
  res_mat <- vector(mode = "numeric", length = nwin)
  if (!quiet) progbar <- txtProgressBar(style = 3)
  for (i in seq(nwin)){
    posns <- which(xpos %in% winmat[i, 1]:winmat[i, 2])
    if (length(posns) < min.snps){
      res_mat[i] <- NA
    } else {
      res_mat[i] <- bitwise.ia(x[, posns], threads = threads)
    }
    if (!quiet){
      setTxtProgressBar(progbar, i/nwin)
    }
  }
  if (!quiet) cat("\n")
  return(res_mat)
}

#==============================================================================#
#' Calculate random samples of the index of association for genlight objects.
#' 
#' Genlight objects can contain millions of loci. Since it does not make much 
#' sense to calculate the index of association over that many loci, this
#' function will randomly sample sites to calculate the index of association.
#' 
#' @param x a genlight object.
#'   
#' @param n.snp the number of snps to be used to calculate standardized index
#' of association.
#' 
#' @param reps the number of times to perform the calculation.
#'   
#' @param threads The maximum number of parallel threads to be used within this 
#'   function. A value of 0 (default) will attempt to use as many threads as
#'   there are available cores/CPUs. In most cases this is ideal. A value of 1
#'   will force the function to run serially, which may increase stability on
#'   some systems. Other values may be specified, but should be used with
#'   caution.
#' 
#' @param quiet if \code{FALSE}, a progress bar will be printed to the screen.
#'   
#'   
#' @note this will calculate the standardized index of association from Agapow
#' 2001. See \code{\link{ia}} for details.
#' 
#' @return Index of association representing the samples in this genlight
#'   object.
#' @author Zhian N. Kamvar, Jonah C. Brooks
#'   
#' @export
#' @seealso \code{\link[adegenet]{genlight}},
#'    \code{\link{snpclone}},
#'    \code{\link{win.ia}},
#'    \code{\link{ia}},
#'    \code{\link{bitwise.dist}}
#' @examples
#' # with structured snps assuming 1e4 positions
#' set.seed(999)
#' x <- glSim(n.ind = 10, n.snp.nonstruc = 5e2, 
#'            n.snp.struc = 5e2, ploidy = 2,
#'            parallel = FALSE)
#' position(x) <- sort(sample(1e4, 1e3))
#' res <- samp.ia(x)
#' hist(res, breaks = "fd")
#' 
#' # with unstructured snps assuming 1e4 positions
#' set.seed(999)
#' x <- glSim(n.ind = 10, n.snp.nonstruc = 1e3, ploidy = 2)
#' position(x) <- sort(sample(1e4, 1e3))
#' res <- samp.ia(x)
#' hist(res, breaks = "fd")
#==============================================================================#
samp.ia <- function(x, n.snp = 100L, reps = 100L, threads = 1L, quiet = FALSE){
  stopifnot(is(x, "genlight"))
  nloc <- nLoc(x)
  res_mat <- vector(mode = "numeric", length = reps)
  if (!quiet) progbar <- txtProgressBar(style = 3)
  for (i in seq(reps)){
    posns <- sample(nloc, n.snp)
    res_mat[i] <- bitwise.ia(x[, posns], threads = threads)
    if (!quiet){
      setTxtProgressBar(progbar, i/reps)
    }
  }
  if (!quiet) cat("\n")
  return(res_mat)
}
# Sat Aug 15 20:02:40 2015 ------------------------------
# 
# This function was used in place of bitwise.ia before it
# was fixed. Since it has no purpose now, it is being 
# commented out, but kept here for reference.
# 
# snpia <- function(x, threads = 1L){
#   nloc <- nLoc(x)
#   nind <- nInd(x)
#   np <- choose(nind, 2)
#   d_mat <- vapply(seq(nloc), function(i) as.vector(bitwise.dist(x[, i], percent = FALSE, threads = threads)), integer(np))
#   D <- rowSums(d_mat) 
#   SD <- sum(D)        
#   Sd <- colSums(d_mat)
#   Sd2 <- colSums(d_mat*d_mat)
#   Vo <- (sum(D*D) - (SD*SD)/np)/np
#   varj <- (Sd2 - (Sd*Sd)/np)/np
#   Ve <- sum(varj)
#   Svarij <- .Call("pairwise_covar", varj, PACKAGE = "poppr")
#   return((Vo - Ve)/(2 * sum(Svarij)))
# }
