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
# #==============================================================================#
# #~ Bootgen object
# #~ 
# #~ An internal object used for bootstrapping. Not intended for user interaction.
# #~ 
# #~ @section Extends: 
# #~ Class \code{"\linkS4class{genind}"}, directly.
# #~ 
# #~ @name bootgen-class
# #~ @rdname bootgen-class
# #~ @export
# #~ @slot replen repeat length of microsatellite loci
# #~ @author Zhian Kamvar
# #~ @import methods
# #==============================================================================#
# setClass("bootgen", 
#          contains = c("genind"),
#          representation = representation(replen = "numeric"),
# )
# 
# #==============================================================================#
# #~ Methods used for the bootgen object. 
# #~ 
# #~ This is not designed for user interaction.
# #~ 
# #~ @rdname bootgen-methods
# #~ @param x a \code{"\linkS4class{bootgen}"} object
# #~ @param i vector of numerics indicating number of individuals desired
# #~ @param j a vector of numerics corresponding to the loci desired.
# #~ @param ... unused.
# #~ @param drop set to \code{FALSE}
# #==============================================================================#
# setMethod(
#   f = "[",
#   signature(x = "bootgen"),
#   definition = function(x, i, j, ..., drop = FALSE){
#     if (missing(i)) i <- TRUE
#     if (missing(j)) j <- TRUE
#     
#     ## Information Gathering
#     j_length <- 1:length(j)
#     if (length(j) > nLoc(x) | any(j > nLoc(x))){
#       stop('subscript out of bounds')
#     }
#     loci <- levels(x@loc.fac)
#     # numbers of alleles per locus
#     locnall <- x@loc.nall[j]
#     names(locnall) <- names(x@loc.nall[j_length])
#     # names of all the alleles. Length of each list element will be equal to locnall
#     allnames <- x@all.names[j]
#     names(allnames) <- names(x@all.names)[j_length]
#     
#     ## Subsetting table columns
#     indices <- unlist(lapply(loci[j], function(locus){res <- which(x@loc.fac %in% locus)}))
#     locnames <- rep(names(allnames), locnall)
#     tabnames <- paste(locnames, unlist(allnames), sep = ".")
#     if (!is.null(pop(x))){
#       res <- truenames(x)$tab[i, indices, drop = drop]
#     } else {
#       res <- truenames(x)[i, indices, drop = drop]
#     }
#     colnames(res) <- tabnames
#     
#     ## Resetting all factors that need to be set. 
#     x@tab <- res
#     x@loc.fac <- factor(rep(names(locnall), locnall))
#     x@loc.names <- names(locnall)
#     x@loc.nall <- locnall
#     x@all.names <- allnames
#     x@replen <- x@replen[j]
#     return(x)
#   }
# )
# 
# #==============================================================================#
# #~ @rdname bootgen-methods
# #==============================================================================#
# setMethod(
#   f = "dim",
#   signature(x = "bootgen"),
#   definition = function(x){
#     return(c(nInd(x), nLoc(x)))
#   }
# )
# 
# 
# #==============================================================================#
# #~ @rdname bootgen-methods
# #~ @param .Object a character, "bootgen"
# #~ @param gen \code{"\linkS4class{genind}"} object
# #~ @param replen a vector of numbers indicating the repeat length for each 
# #~ microsatellite locus.
# #==============================================================================#
# setMethod(
#   f = "initialize",
#   signature = "bootgen",
#   definition = function(.Object, gen, replen){
#     if (missing(gen)) gen <- new("genind")
#     if (missing(replen)){
#       replen <- vapply(gen@all.names, function(y) guesslengths(as.numeric(y)), 1)
#     }
#     lapply(names(gen), function(y) slot(.Object, y) <<- slot(gen, y))
#     slot(.Object, "replen") <- replen
#     return(.Object)
#   })
# 
# 

#==============================================================================#
#' bruvomat object
#' 
#' An internal object used for bruvo's distance. 
#' Not intended for user interaction.
#' 
#' 
#' @name bruvomat-class
#' @rdname bruvomat-class
#' @export
#' @slot mat a matrix of genotypes with one allele per locus. Number of rows will
#' be equal to (ploidy)*(number of loci)
#' @slot replen repeat length of microsatellite loci
#' @slot ploidy the ploidy of the data set
#' @slot ind.names names of individuals in matrix rows.
#' @author Zhian Kamvar
#' @import methods
#==============================================================================#
setClass(
  Class = "bruvomat", 
  representation = representation(
    mat = "matrix", 
    replen = "numeric",
    ploidy = "numeric",
    ind.names = "character"
  ),
   prototype = prototype(
    mat = matrix(ncol = 0, nrow = 0),
    replen = 0,
    ploidy = 2,
    ind.names = "none"
  )
)

#==============================================================================#
#' @rdname bruvomat-methods
#' @param .Object a character, "bruvomat"
#' @param gen \code{"\linkS4class{genind}"} object
#' @param replen a vector of numbers indicating the repeat length for each 
#' microsatellite locus.
#==============================================================================#
setMethod(
  f = "initialize",
  signature = "bruvomat",
  definition = function(.Object, gen, replen){
    if (missing(gen)) gen <- new("genind")
    if (missing(replen)){
      replen <- vapply(gen@all.names, function(y) guesslengths(as.numeric(y)), 1)
    }
    ploid <- ploidy(gen)
    # This controlls for the user correcting missing data using "mean". 
    if (any(!gen@tab %in% c((0:ploid)/ploid, NA))){
      gen@tab[!gen@tab %in% c((0:ploid)/ploid, NA)] <- NA
    }
    # This will check for data that has missing scored as "zero".
    # Data such as this cannot easily be treated by this function and will be
    # sent to the internal function phylo.bruvo.boot
    popcols <- ploid*nLoc(gen)
    if(!any(is.na(gen@tab)) & any(rowSums(gen@tab, na.rm=TRUE) < nLoc(gen))){
      mat1 <- as.matrix.data.frame(genind2df(gen, sep="/", usepop=FALSE))
      mat1[mat1 %in% c("", NA)] <- paste(rep(0, ploid), collapse="/")
      mat2 <- apply(mat1, 1, strsplit, "/")
      mat3 <- apply(as.matrix(t(sapply(mat2, unlist))), 2, as.numeric)
      vec1 <- suppressWarnings(as.numeric(unlist(mat3)))
      pop  <- matrix(vec1, nrow=nInd(gen), ncol=popcols)
    } else {
      popdf <- genind2df(gen, oneColPerAll=TRUE, usepop=FALSE)
      mat1  <- as.matrix.data.frame(popdf)
      pop   <- suppressWarnings(matrix(as.numeric(mat1), ncol=popcols))
    }
    slot(.Object, "mat")       <- pop
    slot(.Object, "replen")    <- replen
    slot(.Object, "ploidy")    <- ploid
    slot(.Object, "ind.names") <- indNames(gen)
    return(.Object)
  }
)

#==============================================================================#
#' @rdname bruvomat-methods
#==============================================================================#
setMethod(
  f = "dim",
  signature(x = "bruvomat"),
  definition = function(x){
    return(c(nrow(x@mat), ncol(x@mat)/x@ploidy))
  }
)

#==============================================================================#
#' Methods used for the bruvomat object. 
#' 
#' This is not designed for user interaction.
#' 
#' @rdname bruvomat-methods
#' @param x a \code{"\linkS4class{bruvomat}"} object
#' @param i vector of numerics indicating number of individuals desired
#' @param j a vector of numerics corresponding to the loci desired.
#' @param ... unused.
#' @param drop set to \code{FALSE}
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "bruvomat"),
  definition = function(x, i, j, ..., drop = FALSE){
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    x@replen <- x@replen[j]
    x@ind.names <- x@ind.names[i]
    cols <- rep(1:ncol(x), each = x@ploidy)
    replacement <- vapply(j, function(ind) which(cols == ind), 1:x@ploidy)
    x@mat <- x@mat[i, as.vector(replacement), drop = FALSE]
    return(x)
  }
)