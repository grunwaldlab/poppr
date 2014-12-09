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
# Implementation of ade4's AMOVA function. Note that this cannot be used at the
# moment to calculate within individual variances. It will either compute a
# distance matrix or take in a distance matrix. Missing data must be treated
# here and is currently treated as extra alleles, but this can be modified by
# the user. Since ade4 needs a euclidean matrix, this function will, by default,
# correct the distance via the cailliez correction, which will add a number to 
# all the distances to satisfy euclidean nature. 
# Clone correction at the lowest level of the hierarchy is possible. 
#
# Note: This takes a nested formula argument. If you want to analyze the
# hierarchy of Year to Population to Subpopulation, you should make sure you
# have the right data frame in your "other" slot and then write the formula
# thusly: ~ Year/Population/Subpopulation
#
# arguments:
#
# x            = a genind object
# hier         = a formula such as ~Pop/Subpop
# clonecorrect = This refers to clone correction of the final output relative to
#                The lowest hierarchical level. FALSE
# within       = should within individual variation be calculated? TRUE
# dist         = A user provided distance matrix NULL
# squared      = Is the distance matrix squared? TRUE
# correction   = A correction for non-euclidean distances provided by ade4.
#                The default, "quasieuclid", seems to give the best results.
# dfname       = the data frame containing the population hierarchy.
#                "population_hierarchy"
# sep          = the separator for the population hierarchy levels. "_"
# missing      = how to deal with missing data. Default is "loci".
# cutoff       = a cutoff for percent missing data to tolerate. 0.05
# quiet        = Should messages be printed? TRUE
#==============================================================================#
#' Perform Analysis of Molecular Variance (AMOVA) on genind or genclone objects.
#' 
#' This function utilizes the ade4 implementation of AMOVA. See 
#' \code{\link[ade4]{amova}} for details on the specific implementation.
#' 
#' @param x a \code{\linkS4class{genind}} or \code{\linkS4class{genclone}}
#'   object
#'   
#' @param hier a hierarchical \code{\link{formula}} that defines your population
#'   hierarchy. (e.g.: ~Population/Subpopulation). \strong{See Details below.}
#'   
#' @param clonecorrect \code{logical} if \code{TRUE}, the data set will be clone
#'   corrected with respect to the lowest level of the hierarchy. The default is
#'   set to \code{FALSE}. See \code{\link{clonecorrect}} for details.
#'   
#' @param within \code{logical}. When this is set to \code{TRUE} (Default), 
#'   variance within individuals are calculated as well. If this is set to
#'   \code{FALSE}, The lowest level of the hierarchy will be the sample level.
#'   See Details below.
#'   
#' @param dist an optional distance matrix calculated on your data.
#'   
#' @param squared if a distance matrix is supplied, this indicates whether or
#'   not it represents squared distances.
#'   
#' @param correction a \code{character} defining the correction method for 
#'   non-euclidean distances. Options are \code{\link[ade4]{quasieuclid}}
#'   (Default), \code{\link[ade4]{lingoes}}, and \code{\link[ade4]{cailliez}}.
#'   See Details below.
#'   
#' @param dfname if the input data set is a \code{\linkS4class{genind}} object, 
#'   specify the name of the data frame in the \code{\link[adegenet]{other}}
#'   slot defining the population hierarchy. Defaults to
#'   \code{"population_hierarchy"}
#'   
#' @param sep A single character used to separate the hierarchical levels. This
#' defaults to "_".
#'   
#' @param missing specify method of correcting for missing data utilizing
#'   options given in the function \code{\link{missingno}}. Default is
#'   \code{"loci"}.
#'   
#' @param cutoff specify the level at which missing data should be
#'   removed/modified. See \code{\link{missingno}} for details.
#'   
#' @param quiet \code{logical} If \code{FALSE} (Default), messages regarding any
#'   corrections will be printed to the screen. If \code{TRUE}, no messages will
#'   be printed.
#'   
#' @return a list of class \code{amova} from the ade4 package. See 
#'   \code{\link[ade4]{amova}} for details.
#' 
#' @details The poppr implementation of AMOVA is a very detailed wrapper for the
#'   ade4 implementation. The output is an \code{\link[ade4]{amova}} class list 
#'   that contains the results in the first four elements. The inputs are contained in the 
#'   last three elements. The inputs required for the ade4 implementation are: 
#'   \enumerate{
#'   \item a distance matrix on all unique genotypes (haplotypes)
#'   \item a data frame defining the hierarchy of the distance matrix 
#'   \item  a genotype (haplotype) frequency table.} 
#'   All of this data can be constructed from a 
#'   \code{\linkS4class{genind}} object, but can be daunting for a novice R 
#'   user. \emph{This function automates the entire process}. Since there are many 
#'   variables regarding genetic data, some points need to be highlighted: 
#'
#'   \subsection{On Hierarchies:}{The hierarchy is defined by different hierarchical
#'   levels that separate your data. In a \code{\linkS4class{genclone}} object,
#'   these levels are inherently defined in the \code{hierarchy} slot. For
#'   \code{\linkS4class{genind}} objects, these levels must be defined in a data
#'   frame located within the \code{\link[adegenet]{other}} slot. It is best
#'   practice to name this data frame \code{"population_hierarchy"}.}
#'
#'   \subsection{On Within Individual Variance:}{ Heterozygosities within diploid
#'   genotypes are sources of variation from within individuals and can be
#'   quantified in AMOVA. When \code{within = TRUE}, poppr will split diploid
#'   genotypes into haplotypes and use those to calculate within-individual
#'   variance. No estimation of phase is made. This acts much like the default
#'   settings for AMOVA in the Arlequin software package. Within individual
#'   variance will not be calculated for haploid individuals or dominant
#'   markers.} 
#'
#'   \subsection{On Euclidean Distances:}{ AMOVA, as defined by
#'   Excoffier et al., utilizes an absolute genetic distance measured in the
#'   number of differences between two samples across all loci. With the ade4
#'   implementation of AMOVA (utilized by poppr), distances must be Euclidean
#'   (due to the nature of the calculations). Unfortunately, many genetic
#'   distance measures are not always euclidean and must be corrected for before
#'   being analyzed. Poppr automates this with three methods implemented in 
#'   ade4, \code{\link{quasieuclid}}, \code{\link{lingoes}}, and 
#'   \code{\link{cailliez}}. The correction of these distances should not 
#'   adversely affect the outcome of the analysis.}
#'   
#' @keywords amova
#' @aliases amova
#' 
#' @references Excoffier, L., Smouse, P.E. and Quattro, J.M. (1992) Analysis of
#' molecular variance inferred from metric distances among DNA haplotypes:
#' application to human mitochondrial DNA restriction data. \emph{Genetics},
#' \strong{131}, 479-491.
#' 
#' @seealso \code{\link[ade4]{amova}} \code{\link{clonecorrect}}
#'   \code{\link{diss.dist}} \code{\link{missingno}}
#'   \code{\link[ade4]{is.euclid}} \code{\link{sethierarchy}}
#' @export
#' @examples
#' data(Aeut)
#' agc <- as.genclone(Aeut)
#' agc
#' amova.result <- poppr.amova(agc, ~Pop/Subpop)
#' amova.result
#' amova.test <- randtest(amova.result) # Test for significance
#' plot(amova.test)
#' amova.test
#' \dontrun{
#' amova.cc.result <- poppr.amova(agc, ~Pop/Subpop, clonecorrect = TRUE)
#' amova.cc.result
#' amova.cc.test <- randtest(amova.cc.result)
#' plot(amova.cc.test)
#' amova.cc.test
#' }
#==============================================================================#
#' @importFrom ade4 amova is.euclid cailliez quasieuclid lingoes
poppr.amova <- function(x, hier = NULL, clonecorrect = FALSE, within = TRUE, 
                        dist = NULL, squared = TRUE, correction = "quasieuclid", 
                        dfname = "population_hierarchy", sep = "_", 
                        missing = "loci", cutoff = 0.05, quiet = FALSE){
  if (!is.genind(x)) stop(paste(substitute(x), "must be a genind object."))
  if (is.null(hier)) stop("A population hierarchy must be specified")
  parsed_hier <- gsub(":", sep, attr(terms(hier), "term.labels"))
  full_hier <- parsed_hier[length(parsed_hier)]
  
  if (is.genclone(x)){
    setpop(x) <- hier
    other(x)[[dfname]] <- gethierarchy(x, hier, combine = FALSE)
  } else {
    if (!dfname %in% names(other(x))){
      stop(paste(dfname, "is not present in the 'other' slot"))
    }
    if (!full_hier %in% names(other(x)[[dfname]])){
      hiers <- all.vars(hier)
      if (!all(hiers %in% names(other(x)[[dfname]]))){
        hier_incompatible_warning(hiers, df)
      }
      x <- splitcombine(x, hier = hiers, dfname = dfname, method = 2)
    } else {
      pop(x) <- other(x)[[dfname]][[full_hier]]
    }
  }
  # Treat missing data. This is a part I do not particularly like. The distance
  # matrix must be euclidean, but the dissimilarity distance will not allow
  # missing data to contribute to the distance. 
  #
  # This is corrected by setting missing to zero: more diversity
  # missing to mean of the columns: indiscreet distances.
  # remove loci at cutoff
  # remove individuals at cutoff
  if (clonecorrect){
    x <- clonecorrect(x, hier = hier, keep = 1:length(all.vars(hier)))
  }
  if (within & ploidy(x) == 2 & check_Hs(x)){
    hier <- update(hier, ~./Individual)
    x    <- pool_haplotypes(x, dfname = dfname)
  }
  x       <- missingno(x, type = missing, cutoff = cutoff, quiet = quiet)
  hierdf  <- make_hierarchy(hier, other(x)[[dfname]])
  xstruct <- make_ade_df(hier, hierdf)
  if (is.null(dist)){
    xdist <- sqrt(diss.dist(clonecorrect(x, hier = NA), percent = FALSE))
  } else {
    datalength <- choose(nInd(x), 2)
    mlgs       <- mlg(x, quiet = TRUE)
    mlglength  <- choose(mlgs, 2)
    if (length(dist) > mlglength & length(dist) == datalength){
      corrected <- .clonecorrector(x)
      xdist     <- as.dist(as.matrix(dist)[corrected, corrected])
    } else if(length(dist) == mlglength){
      xdist <- dist
    } else {
      msg <- paste("\nDistance matrix does not match the data.",
      "\nUncorrected observations expected..........", nInd(x),
      "\nClone corrected observations expected......", mlgs,
      "\nObservations in provided distance matrix...", ceiling(sqrt(length(dist)*2)),
      ifelse(within == TRUE, "\nTry setting within = FALSE.", "\n"))
      stop(msg)
    }
    if (squared){
      xdist <- sqrt(xdist)
    }
  }
  if (!is.euclid(xdist)){
    CORRECTIONS <- c("cailliez", "quasieuclid", "lingoes")
    try(correct <- match.arg(correction, CORRECTIONS))
    if (!exists("correct")){
      stop(not_euclid_msg(correction))
    } else {
      correct_fun <- match.fun(correct)
      if (correct == CORRECTIONS[2]){
        cat("Distance matrix is non-euclidean.\n")
        cat("Utilizing quasieuclid correction method. See ?quasieuclid for details.\n")
        xdist <- correct_fun(xdist)        
      } else {
        xdist <- correct_fun(xdist, print = TRUE, cor.zero = FALSE)        
      }
    }
  }
  xtab    <- t(mlg.table(x, bar = FALSE, quiet = TRUE, mlgsub = unique(mlg.vector(x))))
  xtab    <- as.data.frame(xtab)
  return(ade4::amova(samples = xtab, distances = xdist, structures = xstruct))
}
