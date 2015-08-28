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
#' Round Robin Multilocus Genotypes
#' 
#' This function will mask each locus one by one and then calculate multilocus
#' genotypes from the remaining loci in a round-robin fashion. This is used for
#' calculating the round robin allele frequencies for pgen and psex.
#' 
#' @param gid a genind, genclone, or loci object.
#' 
#' @author Zhian N. Kamvar, Jonah C. Brooks, Stacey Hatfield
#' 
#' @return a matrix of multilocus genotype assignments by masked locus. There 
#'   will be n rows and m columns where n = number of samples and m = number of
#'   loci.
#'
#' @export
#' @examples
#' 
#' # Find out the round-robin multilocus genotype assignments for P. ramorum
#' data(Pram)
#' pmlg_rr <- rrmlg(Pram)
#' head(pmlg_rr)
#' \dontrun{
#' # You can find out how many unique genotypes are found without each locus:
#' 
#' colSums(!apply(pmlg_rr, 2, duplicated))
#' }
#==============================================================================#
rrmlg <- function(gid){
  if (is.genind(gid)){
    gid <- pegas::as.loci(gid)
  }
  the_loci <- attr(gid, "locicol")
  res      <- integer(nrow(gid))
  suppressWarnings(gid <- vapply(gid[the_loci], as.integer, res))
  out <- .Call("mlg_round_robin", gid, PACKAGE = "poppr")
  dimnames(out) <- dimnames(gid)
  return(out)
}

#==============================================================================#
#' Round Robin Allele Frequencies
#' 
#' This function utilizes \code{\link{rrmlg}} to calculate multilocus genotypes 
#' and then subsets each locus by the resulting MLGs to calculate the
#' round-robin allele frequencies used for pgen and psex.
#' 
#' @param gid a genind or genclone object
#' @param res Either "list" (default) or "vector".
#' 
#' @return a vector or list of allele frequencies
#' 
#' @author Zhian N. Kamvar, Jonah C. Brooks, Stacey Hatfield
#' @export
#==============================================================================#
rraf <- function(gid, res = "list"){
  RES     <- c("list", "vector")
  res     <- match.arg(res, RES)
  loclist <- seploc(gid, res = "matrix")
  mlgs    <- rrmlg(gid)
  out     <- lapply(locNames(gid), rrcc, loclist, mlgs)
  names(out) <- locNames(gid)
  if (res == "vector"){
    out <- unlist(out, use.names = FALSE)
    names(out) <- colnames(tab(gid))
  }
  return(out)
}


rrcc <- function(i, loclist, mlgs){
  cc  <- !duplicated(mlgs[, i])
  res <- colMeans(loclist[[i]][cc, , drop = FALSE], na.rm = TRUE)
  res[res < .Machine$double.eps^0.5] <- 1/length(cc)
  return(res)
}