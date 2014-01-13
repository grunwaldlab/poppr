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
#' Genclone class
#' 
#' Genclone is an S4 class that extends the \code{\linkS4class{genind}}
#' from the \emph{\link{adegenet}} package. It will have all of the same
#' attributes as the \code{\linkS4class{genind}}, but it will contain two
#' extra slots that will help retain information about population hierarchies
#' and multilocus genotypes.
#' 
#' @section Extends: 
#' Class \code{"\linkS4class{genind}"}, directly.
#' 
#' @details The genclone class will allow for more optimized methods of clone
#' correcting and analyzing data over multiple levels of population hierarchy.
#' 
#' Previously, for hierarchical analysis to work in a \code{\link{genind}}
#' object, the user had to place a data frame in the \code{\link{other}} slot of
#' the object. The suggested name of the data frame was 
#' \code{population_hierarchy}, and this was used to be able to store the
#' hierarchical information inside the object so that the user did not have to
#' keep track of that information. This method worked, but it became apparent
#' that this method of doing things became a bit confusing to the user as the
#' method for changing the population of an object became:
#' 
#' \code{pop(object) <- other(object)$population_hierarchy$population_name}
#' 
#' That is a lot to keep track of. The new hierarchy slot will allow the user
#' to be able to insert a population hierarchy into the slot and then change it
#' with one function and a formula:
#' 
#' \code{sethierarchy(object, ~Population/Subpopulation)}
#' 
#' making this become slightly more intuitive and tractable.
#' 
#' Calculations of multilocus genotypes is rapid, but unfortunately, the
#' assignments can change if the number of individuals in the population is
#' different. This means that if you analyze the multilocus genotypes for two
#' subpopulations, they will have different multilocus genotype assignments even
#' though they may share multilocus genotypes. The new slot allows us to be able
#' to assign the multilocus genotypes and retain that information no matter how 
#' we subset the data set.
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
#' @rdname genclone-method
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
#' @rdname genclone-method
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
    if (missing(hierarchy)){
      if (is.null(pop(gen))){
        hierarchy <- data.frame()
      } else {
        hierarchy <- data.frame(Pop = pop(gen))
      }
    }
    lapply(names(gen), function(y) slot(.Object, y) <<- slot(gen, y))
    slot(.Object, "mlg")       <- mlg.vector(gen)
    slot(.Object, "hierarchy") <- hierarchy
    return(.Object)
  }
)

#==============================================================================#
#' @rdname genclone-method
#==============================================================================#
setMethod(
  f = "show",
  signature = "genclone",
  definition = function(object){
    callNextMethod(object)
    cat("\rPoppr-specfic elements")
    cat("\n----------------------\n")
    cat("@mlg: a numeric vector designating each of the", length(unique(object@mlg)),
        "multilocus genotypes.\n")
    if (length(object@hierarchy) > 0){
      cat("@hierarchy: a data frame identifying", length(object@hierarchy),
          "hierarchical levels:", names(object@hierarchy), "\n", fill = 80)
    } else {
      cat("@hierarchy: empty")
    }
    
  })

#==============================================================================#
#' Access and manipulate the population hierarchy for genclone objects.
#' 
#' The following methods allow the user to quickly change the hierarchy or
#' population of a genclone object. 
#' 
#' @export 
#' @rdname hierarchy-methods
#' @aliases gethierarchy,genclone-method
#' @param x a genclone object
#' @param formula a nested formula indicating the order of the population
#' hierarchy.
#' @param combine if \code{TRUE}, the levels will be combined according to the
#' formula argument. If it is \code{FALSE}, the levels will not be combined.
#' @docType methods
#==============================================================================#
gethierarchy <- function(x, formula = NULL, combine = TRUE) standardGeneric("gethierarchy")

#' @export
setGeneric("gethierarchy")


setMethod(
  f = "gethierarchy",
  signature = c(x = "genclone"),
  definition = function(x, formula = NULL, combine = TRUE){
    if (is.null(formula)) return(x@hierarchy)
    vars <- all.vars(formula)
    if (any(!vars %in% names(x@hierarchy))){
      stop(hier_incompatible_warning(vars, x@hierarchy))
    }
    if (combine){
      hier <- make_hierarchy(formula, x@hierarchy)
    } else {
      hier <- x@hierarchy[all.vars(formula)]
    }
    return(hier)
  })

#==============================================================================#
#' @export
#' @rdname hierarchy-methods
#' @aliases sethierarchy<-,genclone-method
#' @param value a data frame giving the population hierarchy of each individual.
#' @docType methods
#==============================================================================#
"sethierarchy<-" <- function(x, value) standardGeneric("sethierarchy<-")

#' @export
setGeneric("sethierarchy<-")


setMethod(
  f = "sethierarchy<-",
  signature = c(x = "genclone"),
  definition = function(x, value){
    if (!inherits(value, "data.frame")){
      stop(paste(substitute(value), "is not a data frame"))
    }
    if (nrow(value) != nInd(x)){
      stop("Number of rows in data frame not equal to number of individuals in object.")
    }
    x@hierarchy <- value
    return(x)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases sethierarchy,genclone-method
#' @docType methods
#==============================================================================#
"sethierarchy" <- function(x, value) standardGeneric("sethierarchy<-")

#' @export
setGeneric("sethierarchy")

setMethod(
  f = "sethierarchy",
  signature = c(x = "genclone"),
  definition = function(x, value){
    if (!inherits(value, "data.frame")){
      stop(paste(substitute(value), "is not a data frame"))
    }
    if (nrow(value) != nInd(x)){
      stop("Number of rows in data frame not equal to number of individuals in object.")
    }
    x@hierarchy <- value
    return(x)
  })


#==============================================================================#
#' Manipulate the population factor of genclone objects.
#' 
#' The following methods allow the user to quickly change the population of a 
#' genclone object. 
#' 
#' @export 
#' @rdname population-methods
#' @param x a genclone object
#' @param formula a nested formula indicating the order of the population
#' hierarchy.
#' @param value same as formula
#' @aliases setpop,genclone-method
#' @docType methods
#==============================================================================#
"setpop" <- function(x, formula = NULL) standardGeneric("setpop")

#' @export
setGeneric("setpop")


setMethod(
  f = "setpop",
  signature = c(x = "genclone"),
  definition = function(x, formula = NULL){
    if (is.null(formula) | !is.language(formula)){
      stop(paste(substitute(formula), "must be a valid formula object."))
    }
    vars <- all.vars(formula)
    if (!all(vars %in% names(x@hierarchy))){
      stop(hier_incompatible_warning(vars, x@hierarchy))
    }
    pop(x) <- make_hierarchy(formula, x@hierarchy)[[length(vars)]]
    return(x)
  })

#==============================================================================#
#' @export
#' @rdname population-methods
#' @aliases setpop<-,genclone-method
#' @docType methods
#==============================================================================#
"setpop<-" <- function(x, value) standardGeneric("setpop<-")

#' @export
setGeneric("setpop<-")


setMethod(
  f = "setpop<-",
  signature = c(x = "genclone"),
  definition = function(x, value){
    if (is.null(value) | !is.language(value)){
      stop(paste(substitute(value), "must be a valid value object."))
    }
    vars <- all.vars(value)
    if (!all(vars %in% names(x@hierarchy))){
      stop(hier_incompatible_warning(vars, x@hierarchy))
    }
    pop(x) <- make_hierarchy(value, x@hierarchy)[[length(vars)]]
    return(x)
  })
# #==============================================================================#
# #' Multilocus genotype functions
# #' 
# #' words
# #' 
# #' @export
# #' @rdname mlg-methods
# #' @param pop a genclone object
# #' @param ... other things quiet perhaps
# #' @docType methods
# #==============================================================================#
# setGeneric("mlg", function(pop, ...) standardGeneric("mlg"))
# 
# #' @rdname mlg-methods
# setMethod(
#   f = "mlg",
#   signature = "genclone",
#   definition = function(pop, quiet = FALSE){
#     if (!quiet){
#       cat("#############################\n")
#       cat("# Number of Individuals: ", nInd(pop), "\n")
#       cat("# Number of MLG: ", unique(pop@mlg), "\n")
#       cat("#############################\n")
#     }
#     return(unique(pop@mlg))
#   })


# #==============================================================================#
# #' @export
# #' @rdname mlg-methods
# #' @docType methods
# #==============================================================================#
# setGeneric("mlg.vector", function(pop) standardGeneric("mlg.vector"))
# 
# #' @rdname mlg-methods
# setMethod(
#   f = "mlg.vector",
#   signature = "genclone",
#   definition = function(pop){
#     return(pop@mlg)
#   })