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
#==============================================================================#
#' genclone object
#' 
#' An object based on the genind object that retains multilocus genotype and 
#' population hierarchical information.
#' 
#' @section Extends: 
#' Class \code{"\linkS4class{genind}"}, directly.
#' 
#' @name genclone-class
#' @rdname genclone-class
#' @export
#' @slot mlg a vector representing multilocus genotypes for the data set.
#' @slot hierarchy a data frame containing hierarchical levels.
#' @author Zhian N. Kamvar
#' @import methods
#==============================================================================#
setClass("genclone", 
         contains = c("genind"),
         representation = representation(mlg = "numeric", hierarchy = "data.frame"),
)

#==============================================================================#
#' Methods used for the genclone object
#' 
#' words.
#' 
#' @rdname genclone-methods
#' @param x a genclone object
#' @param i vector of numerics indicating number of individuals desired
#' @param j a vector of numerics corresponding to the loci desired.
#' @param ... passed on to the \code{\linkS4class{genind}} object.
#' @param drop set to \code{FALSE}
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "genclone"),
  definition = function(x, i, j, ..., drop = FALSE){
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    slot(x, "mlg")       <- slot(x, "mlg")[i]
    slot(x, "hierarchy") <- slot(x, "hierarchy")[i, , drop = FALSE]
    callNextMethod()
    return(x)
  }
)

#==============================================================================#
#' @rdname genclone-methods
#' @param .Object a character, "genclone"
#' @param gen \code{"\linkS4class{genind}"} object
#' @param hierarchy a data frame where each row i represents the different
#' population assignments of individual i in the data set. If this is empty, the
#' hierarchy will be created from the population factor.
#==============================================================================#
setMethod(      
  f = "initialize",
  signature = "genclone",
  definition = function(.Object, gen, hierarchy){
    if (missing(gen)) gen <- new("genind")
    if (missing(hier)){
      if (is.null(pop(gen))){
        hier <- data.frame()
      } else {
        hier <- data.frame(Pop = pop(gen))
      }
    }
    lapply(names(gen), function(y) slot(.Object, y) <<- slot(gen, y))
    slot(.Object, "mlg")       <- mlg.vector(gen)
    slot(.Object, "hierarchy") <- hier
    return(.Object)
  }
)