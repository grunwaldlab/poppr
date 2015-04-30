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
#' @param x a genlight, genind, or genclone object. 
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
#' @param differences_only \code{logical}. Determines whether the matrix should
#' count differences or distances. For instance, 0 to 2 would be a distance of 2
#' but a difference of 1.
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
bitwise.dist <- function(x, percent=TRUE, mat=FALSE, missing_match=TRUE, differences_only=FALSE, threads=0){
  stopifnot(class(x)[1] %in% c("genlight", "genclone", "genind", "snpclone"))
  # Stop if the ploidy of the genlight object is not consistent
  stopifnot(min(ploidy(x)) == max(ploidy(x))) 
  # Stop if the ploidy of the genlight object is not haploid or diploid
  stopifnot(min(ploidy(x)) == 2 || min(ploidy(x)) == 1)

  ploid     <- min(ploidy(x))
  ind.names <- indNames(x)
  inds      <- nInd(x)
  numPairs   <- nLoc(x)

  # Use Provesti if this is a genclone or genind object
  if(!is(x, "genlight")){
    dist.mat <- provesti.dist(x)
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
  colnames(dist.mat)            <- ind.names
  rownames(dist.mat)            <- ind.names
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
#' Calculates and returns a matrix of Pgen values for the given genlight object.
#' Each column represents a genotype in the genlight objects, and each row
#' represents a specific, sequential set of 8 base pairs. The values in each cell
#' represent the Pgen value of the corresponding set of base pairs. These values
#' indicate the probability of observing these alleles in a randomly mating 
#' population using estimates derived from the genotypes present in the genlight 
#' object.
#'
#' @param x a genlight object. 
#'
#' @param log a \code{logical} to determine whether the values should be returned 
#'  as percentages or logarithms of percentages. \code{TRUE} is the default, and 
#'  returns the logarithmic values rather than the percentage values. This option has 
#'  a much larger range and is highly recommended. \code{FALSE} returns the percentage 
#'  chance for each genotype to be produced via random mating, rather than the log 
#'  equivalent.
#'
#' @param by.pop a \code{logical} to determine whether allelic frequencies should
#'  be calculated per population (\code{TRUE}, default) or across all populations
#'  in the data (\code{FALSE}).
#'
#' @param window.size an \code{integer} to determine how many SNPs should be
#'  included in each pgen calculation. The default is 1, causing every SNP to
#'  have its own pgen value in the result matrix. Higher values can be used
#'  to reduce matrix size, but may result in precision errors if pgen values are
#'  too small. This argument only affects processing of genlight objects.
#'
#' @return A vector containing one Pgen value for each genotype in the genlight object.
#' @author Zhian N. Kamvar, Jonah Brooks
#' 
#' @export
#==============================================================================#
pgen <- function(x, log=TRUE, by.pop=TRUE, window.size=1) {
  stopifnot(class(x)[1] == "genlight" || is.genind(x))
  # Stop if the ploidy of the genlight object is not consistent
  stopifnot(min(ploidy(x)) == max(ploidy(x))) 
  # Stop if the ploidy of the genlight object is not diploid
  stopifnot(min(ploidy(x)) == 2)
  # Warn if window size is unusable
  if(!is.integer(window.size) && !is.numeric(window.size) && window.size <= 0){
    warning("window.size must be a positive integer or numeric. Using default value of 1.")
    window.size = as.integer(1)
  }
  else {
    window.size = as.integer(window.size)
  }

  if(is(x, "genlight")){
    # Ensure that every SNPbin object has data for both chromosomes
    for(i in 1:length(x$gen)){
      if(length(x$gen[[i]]$snp) == 1){
        x$gen[[i]]$snp <- append(x$gen[[i]]$snp, list(as.raw(rep(0,length(x$gen[[i]]$snp[[1]])))))
      }
    }

    # Ensure there is a population assignment for every genotype
    if(is.null(x$pop) || by.pop == FALSE){
      x$pop <- as.factor(rep(1,length(x$gen)))
    }   
 
    # TODO: Support haploids

    pgen_matrix <- .Call("get_pgen_matrix_genlight",x,window.size)

    if(!log)
    {
      pgen_matrix <- exp(pgen_matrix)
    }

    dim(pgen_matrix) <- c(length(x$gen), ceiling(x$n.loc/window.size))
    rownames(pgen_matrix)            <- x$ind.names
    if(window.size == 1) {
      colnames(pgen_matrix)          <- x$loc.names
    }
  }
  else if(is.genind(x)){
    # Ensure there is a population assignment for every genotype
    if(is.null(x$pop) || by.pop == FALSE){
      x$pop <- as.factor(rep(1,dim(x$tab)[1]))
      x$pop.names <- as.character(1)
    }   
    pops <- pop(x)
    freqs <- makefreq(genind2genpop(x,quiet=TRUE),quiet=TRUE)$tab
    pgen_matrix <- .Call("get_pgen_matrix_genind",x,freqs,pops)
    # TODO: calculate in log form in C, then convert back here.
    if(!log)
    {
      pgen_matrix <- exp(pgen_matrix)
    }
    dim(pgen_matrix) <- c(dim(x$tab)[1],length(x$loc.names))
    colnames(pgen_matrix)            <- x$loc.names
    rownames(pgen_matrix)            <- x$ind.names
  }
  else{
    stop("Object provided must be a genlight, genind, or genclone object.")
  }
  return(pgen_matrix);
}



#==============================================================================#
#' Determines whether openMP is support on this system.
#'
#' @return FALSE if openMP is not supported, TRUE if it is
#' @author Zhian N. Kamvar, Jonah Brooks
#' 
#' @export
#==============================================================================#
poppr_has_parallel <- function(){

  supported <- .Call("omp_test")

  if(supported == 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }

}



#==============================================================================#
#' Calculate the index of association between samples in a genlight object.
#' 
#' TODO: Add description of method
#'
#' @param x a genlight object. 
#'
#' @param indices A vector of integers indicating which loci should be
#'   sampled while calculating the index of association. If NULL (default),
#'   all loci will be used.
#'
#' @param threads The maximum number of parallel threads to be used within this
#'   function. A value of 0 (default) will attempt to use as many threads as there
#'   are available cores/CPUs. In most cases this is ideal. A value of 1 will force
#'   the function to run serially, which may increase stability on some systems.
#'   Other values may be specified, but should be used with caution.
#'
#' @return Index of association representing the samples in this genlight object.
#' @author Zhian N. Kamvar, Jonah Brooks
#' 
#' @export
#' @keywords internal
#==============================================================================#
bitwise.IA <- function(x, indices=NULL, threads=0){
  stopifnot(is(x, "genlight"))
  #TODO: Validate the arguments

  #TODO: Allow for automated index generation, such as random or window based

  #TODO: Call C function and return

  return(0)

}

#==============================================================================#
#' Calculate windows of the index of association for genlight objects.
#' 
#' Genlight objects can contain millions of loci. Since it does not make much 
#' sense to caluclate the index of association over that many loci, this
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
#' @note this will calculate the standardized index of assocation from Agapow
#' 2001. See \code{\link{ia}} for details.
#' 
#' @author Zhian N. Kamvar, Jonah Brooks
#'   
#' @export
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
  nwin <- ceiling(max(xpos)/window)
  winmat <- matrix(window*1:nwin, nrow = nwin, ncol = 2)
  winmat[, 1] <- winmat[, 1] - window + 1
  res_mat <- vector(mode = "numeric", length = nwin)
  if (!quiet) progbar <- txtProgressBar(style = 3)
  for (i in seq(nwin)){
    posns <- which(xpos %in% winmat[i, 1]:winmat[i, 2])
    if (length(posns) < min.snps){
      res_mat[i] <- NA
    } else {
      res_mat[i] <- snpia(x[, posns], threads = threads)
    }
    if (!quiet){
      setTxtProgressBar(progbar, i/nwin)
    }
  }
  return(res_mat)
}

#==============================================================================#
#' Calculate random samples of the index of association for genlight objects.
#' 
#' Genlight objects can contain millions of loci. Since it does not make much 
#' sense to caluclate the index of association over that many loci, this
#' function will randomly sample sites to calculate the index of association.
#' 
#' @param x a genlight object.
#'   
#' @param n.snp the number of snps to be used to calcuate standardized index
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
#' @note this will calculate the standardized index of assocation from Agapow
#' 2001. See \code{\link{ia}} for details.
#' 
#' @return Index of association representing the samples in this genlight
#'   object.
#' @author Zhian N. Kamvar, Jonah Brooks
#'   
#' @export
#' 
#' @examples
#' \dontrun{
#' # with structured snps assuming 1e4 positions
#' set.seed(999)
#' x <- glSim(n.ind = 10, n.snp.nonstruc = 5e2, n.snp.struc = 5e2, ploidy = 2)
#' position(x) <- sort(sample(1e4, 1e3))
#' res <- samp.ia(x)
#' hist(res)
#' 
#' # with unstructured snps assuming 1e4 positions
#' set.seed(999)
#' x <- glSim(n.ind = 10, n.snp.nonstruc = 1e3, ploidy = 2)
#' position(x) <- sort(sample(1e4, 1e3))
#' res <- samp.ia(x)
#' hist(res)
#' }
#==============================================================================#
samp.ia <- function(x, n.snp = 100L, reps = 100L, threads = 1L, quiet = FALSE){
  stopifnot(is(x, "genlight"))
  nloc <- nLoc(x)
  res_mat <- vector(mode = "numeric", length = reps)
  if (!quiet) progbar <- txtProgressBar(style = 3)
  for (i in seq(reps)){
    posns <- sample(nloc, n.snp)
    res_mat[i] <- snpia(x[, posns], threads = threads)
    if (!quiet){
      setTxtProgressBar(progbar, i/reps)
    }
  }
  return(res_mat)
}

snpia <- function(x, threads = 1L){
  nloc <- nLoc(x)
  nind <- nInd(x)
  np <- choose(nind, 2)
  d_mat <- vapply(seq(nloc), function(i) as.vector(bitwise.dist(x[, i], percent = FALSE, threads = threads)), integer(np))
  D <- rowSums(d_mat)
  SD <- sum(D)
  Sd <- colSums(d_mat)
  Sd2 <- colSums(d_mat*d_mat)
  Vo <- (sum(D*D) - (SD*SD)/np)/np
  varj <- (Sd2 - (Sd*Sd)/np)/np
  Ve <- sum(varj)
  Svarij <- .Call("pairwise_covar", varj, PACKAGE = "poppr")
  return((Vo - Ve)/(2 * sum(Svarij)))
}
