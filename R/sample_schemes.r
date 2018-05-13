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
#' Shuffle individuals in a \code{\linkS4class{genclone}} or
#' \code{\linkS4class{genind}} object independently over each locus.
#' 
#' @param pop a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}}
#'   object
#'   
#' @param method an integer between 1 and 4. See details below.
#'   
#' @return a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} object
#'   shuffled by a specified method
#'   
#' @details This function will shuffle each locus in the data set independently 
#'   of one another, rendering them essentially unlinked. The following methods 
#'   are available to shuffle your data: \enumerate{ \item \strong{Permute 
#'   Alleles} This will redistribute all alleles in the sample throughout the 
#'   locus. Missing data is fixed in place. This maintains allelic structure, 
#'   but heterozygosity is variable. \item \strong{Parametric Bootstrap} This 
#'   will redistribute available alleles within the locus based on their allelic
#'   frequencies. This means that both the allelic state and heterozygosity will
#'   vary. The resulting data set will not have missing data. \item 
#'   \strong{Non-Parametric Bootstrap} This will shuffle the allelic state for 
#'   each individual. Missing data is fixed in place. \item \strong{Multilocus 
#'   Style Permutation} This will shuffle the genotypes at each locus, 
#'   maintaining the heterozygosity and allelic structure. }
#'   
#'   
#' @export
#' @author Zhian N. Kamvar
#'   
#' @references  Paul-Michael Agapow and Austin Burt. 2001. Indices of multilocus
#'   linkage disequilibrium. \emph{Molecular Ecology Notes}, 1(1-2):101-102
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
#' # Permute Alleles: maintain allelic state; heterozygosity varies.
#' summary(shufflepop(Zebu, method=1))
#' \dontrun{
#' # Parametric Bootstrap: do not maintain allelic state or heterozygosity
#' summary(shufflepop(Zebu, method=2))
#' 
#' # Non-Parametric Bootstrap: do not maintain allelic state or heterozygosity.
#' summary(shufflepop(Zebu, method=3))
#' 
#' # Multilocus Style: maintain allelic state and heterozygosity.
#' summary(shufflepop(Zebu, method=4))
#' }
#==============================================================================#
shufflepop <- function(pop, method=1){
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  if (all((1:4)!=method)) {
    msg <- paste("Method", method, "is not defined. Please choose a defined method:\n\n",
                 "1 = Permute Alleles (maintain allelic structure)\n",
                 "2 = Parametric Bootstrap (simulate new population based on allelic frequency)\n",
                 "3 = Non-Parametric Bootstrap (simulate new population)\n",
                 "4 = Multilocus style (maintain heterozygosity and allelic structure)\n")
    stop(msg)
  }
  if(pop@type == "PA"){
    if(method == 1 | method == 4){
      pop@tab <- vapply(1:ncol(tab(pop)),
                        function(x) sample(tab(pop)[, x]), pop@tab[, 1])
    }
    else if(method == 2){
      paramboot <- function(x){
        one <- mean(tab(pop)[, x], na.rm=TRUE)
        zero <- 1 - one
        return(sample(c(1L ,0L), length(tab(pop)[, x]), prob=c(one, zero), replace=TRUE))
      }
      pop@tab <- vapply(1:ncol(tab(pop)), paramboot, tab(pop)[, 1])
    }
    else if(method == 3){
      pop@tab <- vapply(1:ncol(tab(pop)),
                        function(x) sample(tab(pop)[, x], replace=TRUE), pop@tab[, 1])
    }
  } else {
    addpop <- function(locus = "L1", pop, method=method){
      pop@tab[, locFac(pop) %in% locus] <<- .locus.shuffler(pop[, loc = locus], method=method)@tab
    }
    invisible(lapply(locNames(pop), addpop, pop, method))
  }
  return(pop)
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
  if(!quiet) progbar <- dplyr::progress_estimated(iterations)
  for (c in seq(iterations)){
    IarD <- .Ia.Rd(.all.shuffler(pop, type, method=method), missing=missing)   
    sample.data$Ia[c]    <- IarD[1]
    sample.data$rbarD[c] <- IarD[2]
    if (!quiet) print(progbar$tick())
  }

  if(!quiet){
    print(progbar$stop())
    cat("\n")
  }
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
      pop@tab <- vapply(1:ncol(tab(pop)),
                        function(x) sample(tab(pop)[, x]), pop@tab[, 1])
    } else if(method == 2) {
      paramboot <- function(x){
        one <- mean(tab(pop)[, x], na.rm=TRUE)
        zero <- 1-one
        return(sample(c(1L ,0L), length(tab(pop)[, x]), prob=c(one, zero), replace=TRUE))
      }
      pop@tab <- vapply(1:ncol(tab(pop)), paramboot, pop@tab[, 1])
    } else if(method == 3) {
      pop@tab <- vapply(1:ncol(tab(pop)),
                        function(x) sample(tab(pop)[, x], replace=TRUE), pop@tab[, 1])
    }
  } else {
    pop <- lapply(pop, .locus.shuffler, method = method)
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
  if( ncol(tab(pop)) > 1 ){

    if (method == 4){

      # Maintenence of the heterozygotic and allelic structure.
      pop@tab <- .both.shuff(tab(pop))

    } else if(method == 1){

      # Maintenece of only the allelic structure. Heterozygosity can fluctuate.
      pop@tab <- .permut.shuff(tab(pop))

    } else if(method == 2){

      # Parametric Bootstraping where both heterozygosity and allelic structure
      # can change based on allelic frequency.
      weights <- colMeans(tab(pop), na.rm = TRUE)
      theSize <- sample(ploidy(pop), 1)
      pop@tab <- t(stats::rmultinom(nrow(tab(pop)), size = theSize, prob = weights))

    } else if(method == 3){

      # Non-Parametric Bootstrap.
      weights <- rep(1, ncol(tab(pop)))
      theSize <- sample(ploidy(pop), 1)
      pop@tab <- t(stats::rmultinom(nrow(tab(pop)), size = theSize, prob = weights))

    }
  }
  return(pop)
}

#==============================================================================#
#
# A function for true permutation shuffling. This will retain the same number of
# observed allelic states as the original sample, but differ in the number of
# heterozygous genotypes. 
#
#==============================================================================#
  
.permut.shuff <- function(mat){
  bucket     <- colSums(mat, na.rm = TRUE)
  ploidy     <- rowSums(mat, na.rm = TRUE)
  bucketlist <- as.integer(sample(rep(1:length(bucket), bucket))) - 1L
  mat        <- .Call("permute_shuff", mat, bucketlist, ploidy, PACKAGE = "poppr")
  return(mat)
}

#==============================================================================#
# shuffling both alleles at a locus, thus preserving the number of allelic
# states as well as the heterozygosity. 
#==============================================================================#
.both.shuff <- function(mat){
  typed <- which(!is.na(rowSums(mat)))
  mat[typed, ] <- mat[sample(typed), ]
  return(mat)
}

