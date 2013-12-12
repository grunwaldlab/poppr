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
#' Shuffle individuals in a \code{\link{genind}} object independently over each
#' locus.
#'
#' @param pop a \code{\link{genind}} object
#'
#' @param method an integer between 1 and 4. See details below.
#'
#' @return a \code{\link{genind}} object shuffled by a specified method
#'
#' @section Details: This function will shuffle each locus in the data set independently 
#' of one another, rendering them essentially unlinked. The following methods
#' are available to shuffle your data:
#' \enumerate{ 
#'  \item \strong{Permute Alleles} This will redistribute all alleles in the
#' sample throughout the locus. Missing data is fixed in place. This maintains
#' allelic structure, but heterozygosity is variable.
#'  \item \strong{Parametric Bootstrap} This will redistribute available alleles
#' within the locus based on their allelic frequencies. This means that both the
#' allelic state and heterozygosity will vary. The resulting data set will not
#' have missing data.
#'  \item \strong{Non-Parametric Bootstrap} This will shuffle the
#' allelic state for each individual. Missing data is fixed in place.
#'  \item \strong{Multilocus Style Permutation} This will shuffle the genotypes 
#' at each locus, maintaining the heterozygosity and allelic structure.
#' } 
#'
#'
#' @export
#' @author Zhian N. Kamvar
#'
#' @references  Paul-Michael Agapow and Austin Burt. 2001. Indices of multilocus 
#' linkage disequilibrium. \emph{Molecular Ecology Notes}, 1(1-2):101-102 
#'
#' @examples
#' # load the microbov dataset
#' data(microbov)
#' # Let's look at a single population for now. Howsabout Zebu
#' Zebu <- popsub(microbov, "Zebu")
#' summary(Zebu)
#'
#' # Take note of the Number of alleles per population and the Observed
#' # heterozygosity as we go through each method.
#' 
#' # Multilocus Style: maintain allelic state and heterozygosity.
#' summary(shufflepop(Zebu, method=1))
#' \dontrun{
#' # Permute Alleles: maintain allelic state; heterozygosity varies.
#' summary(shufflepop(Zebu, method=2))
#'
#' # Parametric Bootstrap: do not maintain allelic state or heterozygosity
#' summary(shufflepop(Zebu, method=3))
#' 
#' # Non-Parametric Bootstrap: do not maintain allelic state or heterozygosity.
#' summary(shufflepop(Zebu, method=4))
#' }
#==============================================================================#
shufflepop <- function(pop, method=1){
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  if (all((1:4)!=method)) {
    cat("1 = Permute Alleles (maintain allelic structure)\n")
    cat("2 = Parametric Bootstrap (simulate new population based on allelic frequency)\n")
    cat("3 = Non-Parametric Bootstrap (simulate new population)\n")
    cat("4 = Multilocus style (maintain heterozygosity and allelic structure)\n")
    cat("Select an integer (1, 2, 3, or 4): ")
    method <- as.integer(readLines(n = 1))
  }
  if (all((1:4)!=method)){
    stop ("Non convenient method number")
  }
  if(pop@type == "PA"){
    if(method == 1 | method == 4){
      pop@tab <- vapply(1:ncol(pop@tab),
                        function(x) sample(pop@tab[, x]), pop@tab[, 1])
    }
    else if(method == 2){
      paramboot <- function(x){
        one <- mean(pop@tab[, x], na.rm=TRUE)
        zero <- 1-one
        return(sample(c(1,0), length(pop@tab[, x]), prob=c(one, zero), replace=TRUE))
      }
      pop@tab <- vapply(1:ncol(pop@tab), paramboot, pop@tab[, 1])
    }
    else if(method == 3){
      pop@tab <- vapply(1:ncol(pop@tab),
                        function(x) sample(pop@tab[, x], replace=TRUE), pop@tab[, 1])
    }
  } else {
    addpop <- function(locus="L1", pop, method=method){
      pop@tab[, pop@loc.fac %in% locus] <<- .locus.shuffler(pop[, loc=locus], method=method)@tab
    }
    invisible(lapply(names(pop@loc.names), addpop, pop, method))
  }
  return(pop)
}
shufflefunk <- function(pop, FUN, sample=1, method=1, ...){
  FUN <- match.fun(FUN)
  lapply(1:sample, function(x) FUN(shufflepop(pop, method=method), ...))
}
#==============================================================================#
# .sampling will reshuffle the alleles per individual, per locus via the 
# .single.sampler function, which is described below. It will then calculate the
# Index of Association for the resampled population for the number of times
# indicated in "iterations".
#==============================================================================#
.sampling <- function(pop, iterations, quiet=FALSE, missing="ignore", type=type, 
                      method=1){ 
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  if(!is.list(pop)){
    if(type=="PA"){
      .Ia.Rd <- .PA.Ia.Rd
    }
  }
	sample.data <- data.frame(list(Ia = vector(mode = "numeric", 
                                             length = iterations),
                                 rbarD = vector(mode = "numeric", 
                                                length = iterations)
                                 )
                            )
  if(!quiet) progbar <- txtProgressBar(style = 3)
	for (c in 1:iterations){
    IarD <- .Ia.Rd(.all.shuffler(pop, type, method=method), missing=missing)
    sample.data$Ia[c]    <- IarD[1]
    sample.data$rbarD[c] <- IarD[2]
    if (!quiet){
      setTxtProgressBar(progbar, c/iterations)
    }
  }
  if(!quiet) close(progbar)
	return(sample.data)
}

#==============================================================================#
# pop = a list of genind objects with one locus each.
# 
# This function shuffles alleles at each locus over all loci. 
#==============================================================================#

.all.shuffler <- function(pop, type=type, method=1){
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  if(type=="PA"){
    if(method == 1 | method == 4){
      pop@tab <- vapply(1:ncol(pop@tab),
                        function(x) sample(pop@tab[, x]), pop@tab[, 1])
    }
    else if(method == 2){
      paramboot <- function(x){
        one <- mean(pop@tab[, x], na.rm=TRUE)
        zero <- 1-one
        return(sample(c(1,0), length(pop@tab[, x]), prob=c(one, zero), replace=TRUE))
      }
      pop@tab <- vapply(1:ncol(pop@tab), paramboot, pop@tab[, 1])
    }
    else if(method == 3){
      pop@tab <- vapply(1:ncol(pop@tab),
                        function(x) sample(pop@tab[, x], replace=TRUE), pop@tab[, 1])
    }
  } 
  else {
    pop <- lapply(pop, .locus.shuffler, method=method)
  }
  return(pop)
}

#==============================================================================#
# This function shuffles the alleles at each locus per individual.
#==============================================================================# 

.locus.shuffler <- function(pop, method=1){
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  # if the locus is homozygous for all individuals, shuffling will not occur.
  if( ncol(pop@tab) > 1 ){

# Maintenence of the heterozygotic and allelic structure.
    if(method == 4)
      pop@tab <- .both.shuff(pop@tab)

# Maintenece of only the allelic structure. Heterozygosity can fluctuate.
    else if(method == 1)
      pop@tab <- .permut.shuff(pop@tab)

# Parametric Bootstraping where both heterozygosity and allelic structure can
# change based on allelic frequency. 
    else if(method == 2){
      weights <- colMeans(pop@tab, na.rm = TRUE)
      #pop@tab  <- t(apply(pop@tab, 1, .diploid.shuff, weights))
      pop@tab <- t(rmultinom(nrow(pop@tab), size = ploidy(pop), prob = weights))/ ploidy(pop)
    }
# Non-Parametric Bootstrap.
    else if(method == 3){
      weights <- rep(1, ncol(pop@tab))
#       pop@tab  <- t(apply(pop@tab, 1, .diploid.shuff, weights))
      pop@tab <- t(rmultinom(nrow(pop@tab), size = ploidy(pop), prob = weights))/ ploidy(pop)
    }
# Maintaining heterozygosity.    
#    if(method == 5){
#      temp <- vapply(1:nrow(pop@tab), function(x) sample(pop@tab[x,]), pop@tab[1,])
#      pop@tab <- t(temp)
#    }
  }
  return(pop)
}



#==============================================================================# 
# This is a parametric boostrap resampling method, It works on one single
# genotype at each individual. In order to make this work you need to take your
# separate loci, then apply this function over all individuals at that locus. 
# this way, it is independent for each individual and each locus. 
#
# vec is a vector of possible alleles for the locus.
#
# weights is a corresponding vector giving the allelic frequency for each allele
# at that locus in that population. The sum of the frequencies should be 1. 
#==============================================================================# 
.diploid.shuff <- function(vec, weights){
  # Dealing with missing values is probably not necessary for a parametric
  # bootstrap resampling. 
  #
  #
  # Missing values need to be handled by preserving the structure.   
  #if(any(is.na(vec))){
  #  if(any(!is.na(vec))){
  #    vec <- rep(NA, length(vec))
  #    vec[sample(length(vec), 1, prob = weights)] <- 0.5
  #  }
  #  return(vec)
  #}   
  # samples potential alleles from the locus
  temp <- sample(length(vec), 2, replace=TRUE, prob = weights)
  #cat("Weights: ",weights, "\n")
  #cat("Sampled: ",temp, "\n")     
  vec <- vector(mode="numeric", length=length(vec))
  # the first allele assigned as heterozygote
  vec[temp[1]] <- 0.5
  # if the second allele is the same as the first, it will become homozygous
  vec[temp[2]] <- vec[temp[2]] + 0.5
  return(vec)
}

#==============================================================================#
#
# A function for true permutation shuffling. This will retain the same number of
# observed allelic states as the original sample, but differ in the number of
# heterozygous genotypes. 
#
#==============================================================================#
  
.permut.shuff <- function(mat){
  # gathering the counts of allelic states (represented as column numbers)
  bucket     <- as.vector(colSums(mat, na.rm=T) * 2)
  bucketlist <- rep(1:length(bucket), bucket)
  # reshuffling the allelic states
  samplist <- bucketlist[sample(length(bucketlist))]
  newmat   <- mat  
  newmat[which(!is.na(mat))] <- 0
  
  # Treating Missing Data
  #
  # go on to the next locus if all are NA
  if(all(is.na(mat))){
    next
  }
  # placing the individuals with missing data in the sample list so that it is
  # the same length as the number of rows * 2 to avoid errors.
  if(any(is.na(mat))){
    temp <- vector(mode="numeric", length=nrow(mat)*2)
    miss <- which(is.na(mat[, 1]))
    temp[c(-miss*2, -((miss*2)-1))] <- samplist
    samplist <- temp
  }

  # A function that will repopulate the new matrix with the shuffled alleles
  # x is a list of even numbers, newmat is an empty matrix, and samplist is 
  # a list of shuffled alleleic states.
  populate.mat <- function(x, newmat, samplist){
    if (samplist[x] == samplist[x-1]){
      newmat[x/2, samplist[x-1]] <<- 1
    }
    else {
      newmat[x/2, samplist[x]]   <<- 0.5
      newmat[x/2, samplist[x-1]] <<- 0.5
    }
  }
  x <- which(1:(nrow(mat)*2)%%2==0)
  lapply(x, populate.mat, newmat, samplist)
  return(newmat)
}

#==============================================================================#
# shuffling both alleles at a locus, thus preserving the number of allelic
# states as well as the heterozygosity. 
#==============================================================================#
.both.shuff <- function(mat){
  mat <- mat[sample(nrow(mat)), ]
  return(mat)
}
