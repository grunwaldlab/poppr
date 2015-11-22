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
#' Perform Analysis of Molecular Variance (AMOVA) on genind or genclone objects.
#' 
#' This function simplifies the process necessary for performing AMOVA in R. It
#' gives user the choice of utilizing either the \pkg{ade4} or the \pkg{pegas}
#' implementation of AMOVA. See \code{\link[ade4]{amova}} (ade4) and
#' \code{\link[pegas]{amova}} (pegas) for details on the specific
#' implementation.
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
#' @param dist an optional distance matrix calculated on your data. If this is
#'   set to \code{NULL} (default), the raw pairwise distances will be calculated
#'   via \code{\link{diss.dist}}.
#'   
#' @param squared if a distance matrix is supplied, this indicates whether or
#'   not it represents squared distances.
#'   
#' @param correction a \code{character} defining the correction method for 
#'   non-euclidean distances. Options are \code{\link[ade4]{quasieuclid}}
#'   (Default), \code{\link[ade4]{lingoes}}, and \code{\link[ade4]{cailliez}}.
#'   See Details below.
#'   
#' @param sep Deprecated. As of poppr version 2, this argument serves no purpose.
#' 
#' @param filter \code{logical} When set to \code{TRUE}, mlg.filter will be run 
#'   to determine genotypes from the distance matrix. It defaults to 
#'   \code{FALSE}. You can set the parameters with \code{algorithm} and
#'   \code{threshold} arguments. Note that this will not be performed when 
#'   \code{within = TRUE}. Note that the threshold should be the number of
#'   allowable substitutions if you don't supply a distance matrix.
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
#' @param method Which method for calculating AMOVA should be used? Choices 
#'   refer to package implementations: "ade4" (default) or "pegas". See details
#'   for differences.
#'   
#' @param nperm the number of permutations passed to the pegas implementation of
#'   amova.
#'   
#' @inheritParams mlg.filter
#'   
#' @return a list of class \code{amova} from the ade4 package. See 
#'   \code{\link[ade4]{amova}} for details.
#' 
#' @details The poppr implementation of AMOVA is a very detailed wrapper for the
#'   ade4 implementation. The output is an \code{\link[ade4]{amova}} class list 
#'   that contains the results in the first four elements. The inputs are
#'   contained in the last three elements. The inputs required for the ade4
#'   implementation are:
#'   \enumerate{
#'   \item a distance matrix on all unique genotypes (haplotypes)
#'   \item a data frame defining the hierarchy of the distance matrix 
#'   \item  a genotype (haplotype) frequency table.} 
#'   All of this data can be constructed from a \code{\linkS4class{genind}}
#'   object, but can be daunting for a novice R user. \emph{This function
#'   automates the entire process}. Since there are many variables regarding
#'   genetic data, some points need to be highlighted:
#'   
#'   \subsection{On Hierarchies:}{The hierarchy is defined by different 
#'   population strata that separate your data hierarchically. These strata are
#'   defined in the \strong{strata} slot of \code{\linkS4class{genind}} and
#'   \code{\linkS4class{genclone}}} objects. They are useful for defining the
#'   population factor for your data. See the function \code{\link{strata}} for
#'   details on how to properly define these strata.
#'
#'   \subsection{On Within Individual Variance:}{ Heterozygosities within
#'   diploid genotypes are sources of variation from within individuals and can
#'   be quantified in AMOVA. When \code{within = TRUE}, poppr will split diploid
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
#'   \subsection{On Filtering:}{ Filtering multilocus genotypes is performed by
#'   \code{\link{mlg.filter}}. This can necessarily only be done AMOVA tests 
#'   that do not account for within-individual variance. The distance matrix used
#'   to calculate the amova is derived from using \code{\link{mlg.filter}} with
#'   the option \code{stats = "distance"}, which reports the distance between 
#'   multilocus genotype clusters. One useful way to utilize this feature is to
#'   correct for genotypes that have equivalent distance due to missing data. 
#'   (See example below.)}
#'   
#'   \subsection{On Methods:}{ Both \pkg{ade4} and \pkg{pegas} have 
#'   implementations of AMOVA, both of which are appropriately called "amova". 
#'   The ade4 version is faster, but there have been questions raised as to the 
#'   validity of the code utilized. The pegas version is slower, but careful 
#'   measures have been implemented as to the accuracy of the method. It must be
#'   noted that there appears to be a bug regarding permuting analyses where 
#'   within individual variance is accounted for (\code{within = TRUE}) in the 
#'   pegas implementation. If you want to perform permutation analyses on the 
#'   pegas implementation, you must set \code{within = FALSE}. In addition,
#'   while clone correction is implemented for both methods, filtering is only
#'   implemented for the ade4 version.}
#'   
#' @keywords amova
#' @aliases amova
#' 
#' @references Excoffier, L., Smouse, P.E. and Quattro, J.M. (1992) Analysis of
#' molecular variance inferred from metric distances among DNA haplotypes:
#' application to human mitochondrial DNA restriction data. \emph{Genetics},
#' \strong{131}, 479-491.
#' 
#' @seealso \code{\link[ade4]{amova}} (ade4) \code{\link[pegas]{amova}} (pegas) 
#'   \code{\link{clonecorrect}} \code{\link{diss.dist}} \code{\link{missingno}} 
#'   \code{\link[ade4]{is.euclid}} \code{\link{strata}}
#' @export
#' @examples
#' data(Aeut)
#' strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
#' agc <- as.genclone(Aeut)
#' agc
#' amova.result <- poppr.amova(agc, ~Pop/Subpop)
#' amova.result
#' amova.test <- randtest(amova.result) # Test for significance
#' plot(amova.test)
#' amova.test
#' 
#' \dontrun{
#' 
#' # You can get the same results with the pegas implementation
#' amova.pegas <- poppr.amova(agc, ~Pop/Subpop, method = "pegas")
#' amova.pegas
#' amova.pegas$varcomp/sum(amova.pegas$varcomp)
#' 
#' # Clone correction is possible
#' amova.cc.result <- poppr.amova(agc, ~Pop/Subpop, clonecorrect = TRUE)
#' amova.cc.result
#' amova.cc.test <- randtest(amova.cc.result)
#' plot(amova.cc.test)
#' amova.cc.test
#' 
#' 
#' # Example with filtering
#' data(monpop)
#' splitStrata(monpop) <- ~Tree/Year/Symptom
#' poppr.amova(monpop, ~Symptom/Year) # gets a warning of zero distances
#' poppr.amova(monpop, ~Symptom/Year, filter = TRUE, threshold = 0.1) # no warning
#' }
#==============================================================================#
#' @importFrom ade4 amova is.euclid cailliez quasieuclid lingoes
poppr.amova <- function(x, hier = NULL, clonecorrect = FALSE, within = TRUE, 
                        dist = NULL, squared = TRUE, correction = "quasieuclid", 
                        sep = "_", filter = FALSE, threshold = 0, 
                        algorithm = "farthest_neighbor", missing = "loci", 
                        cutoff = 0.05, quiet = FALSE, 
                        method = c("ade4", "pegas"), nperm = 0){
  if (!is.genind(x)) stop(paste(substitute(x), "must be a genind object."))
  if (is.null(hier)) stop("A population hierarchy must be specified")
  methods <- c("ade4", "pegas")
  method <- match.arg(method, methods)
  setPop(x) <- hier
  if (filter && (!within | all(ploidy(x) == 1) | !check_Hs(x) | x@type == "PA")){
    
    if (!is.genclone(x)){
      x <- as.genclone(x)
    }
    if (!is(x@mlg, "MLG")){
      x@mlg <- new("MLG", x@mlg)
    }
    if (!quiet){
      message("Filtering ...")
      message("Original multilocus genotypes ... ", nmll(x, "original"))
    }
    if (is.null(dist)){
      nulldist <- TRUE
      dist <- diss.dist(x, percent = FALSE)
    } else {
      nulldist <- FALSE
    }
    filt_stats <- mlg.filter(x, threshold = threshold, algorithm = algorithm,
                             distance = dist, stats = "ALL")
    if (nulldist){
      filt_stats$DISTANCE <- sqrt(filt_stats$DISTANCE)
      squared <- FALSE
    }
    dist <- as.dist(filt_stats$DISTANCE)
    # Forcing this. Probably should make an explicit method for this.
    x@mlg@mlg["contracted"]     <- filt_stats$MLGS
    distname(x@mlg)             <- "diss.dist"
    distalgo(x@mlg)             <- algorithm
    cutoff(x@mlg)["contracted"] <- threshold
    mll(x) <- "contracted"
    if (!quiet) message("Contracted multilocus genotypes ... ", nmll(x))
  }
  if (clonecorrect){
    x <- clonecorrect(x, strata = hier, keep = 1:length(all.vars(hier)))
  }
  if (within & all(ploidy(x) == 2) & check_Hs(x) & x@type != "PA"){
    hier <- update(hier, ~./Individual)
    x    <- pool_haplotypes(x)
  }
  # Treat missing data. This is a part I do not particularly like. The distance
  # matrix must be euclidean, but the dissimilarity distance will not allow
  # missing data to contribute to the distance. 
  #
  # This is corrected by setting missing to zero: more diversity
  # missing to mean of the columns: indiscreet distances.
  # remove loci at cutoff
  # remove individuals at cutoff
  x       <- missingno(x, type = missing, cutoff = cutoff, quiet = quiet)
  hierdf  <- strata(x, formula = hier)
  if (method == "ade4") xstruct <- make_ade_df(hier, hierdf)
  if (is.null(dist)){
    if (method == "ade4"){
      xdist <- sqrt(diss.dist(clonecorrect(x, strata = NA), percent = FALSE))
    } else {
      xdist <- sqrt(diss.dist(x, percent = FALSE))
    }
  } else {
    datalength <- choose(nInd(x), 2)
    mlgs       <- mlg(x, quiet = TRUE)
    mlglength  <- choose(mlgs, 2)
    if (length(dist) > mlglength & length(dist) == datalength){
      corrected <- TRUE
      if (method == "ade4"){
        corrected <- .clonecorrector(x)
      }
      xdist     <- as.dist(as.matrix(dist)[corrected, corrected])
    } else if(length(dist) == mlglength){
      xdist <- dist
    } else {
      distobs <- ceiling(sqrt(length(dist)*2))
      msg <- paste("\nDistance matrix does not match the data.\n",
      "\n\tUncorrected observations expected..........", nInd(x),
      "\n\tClone corrected observations expected......", mlgs,
      "\n\tObservations in provided distance matrix...", distobs,
      ifelse(within == TRUE, "\n\n\tTry setting within = FALSE.", "\n"))
      stop(msg)
    }
    if (squared){
      xdist <- sqrt(xdist)
    }
  }
  if (!is.euclid(xdist)){
    CORRECTIONS <- c("cailliez", "quasieuclid", "lingoes")
    try(correct <- match.arg(correction, CORRECTIONS), silent = TRUE)
    if (!exists("correct")){
      stop(not_euclid_msg(correction))
    } else {
      correct_fun <- match.fun(correct)
      if (correct == CORRECTIONS[2]){
        message("Distance matrix is non-euclidean.")
        message(c("Utilizing quasieuclid correction method.",
                  " See ?quasieuclid for details."))
        xdist <- correct_fun(xdist)        
      } else {
        xdist <- correct_fun(xdist, print = TRUE, cor.zero = FALSE)        
      }
    }
  }
  if (method == "ade4"){
    allmlgs <- unique(mlg.vector(x))
    xtab    <- t(mlg.table(x, plot = FALSE, quiet = TRUE, mlgsub = allmlgs))
    xtab    <- as.data.frame(xtab)
    return(ade4::amova(samples = xtab, distances = xdist, structures = xstruct))
  } else {
    form <- paste(all.vars(hier), collapse = "/")
    hier <- as.formula(paste("xdist ~", form))
    return(pegas::amova(hier, data = hierdf, nperm = nperm, is.squared = FALSE))
  }
}
