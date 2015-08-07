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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Below are methods for the MLG class. I will break up the sections into generic
# methods that are needed for them to make sense in the R world and custom
# generic methods that aren't generic in R, but need to be to make these work.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#==============================================================================#



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GENERIC METHODS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#==============================================================================#
#' Methods used for MLG objects
#' 
#' Default methods for accessing and subsetting MLG objects.
#' 
#' @rdname MLG-method
#' @param .Object a character, "MLG"
#' @param mlg a vector where each element assigns the multilocus genotype of
#' that individual in the data set. 
#' @keywords internal
#' @docType methods
#' @author Zhian N. Kamvar
#==============================================================================#
setMethod(      
  f = "initialize",
  signature("MLG"),
  definition = function(.Object, mlg){
    if (missing(mlg)){
      return(.Object)
    }
    if (is.vector(mlg)){
      mlg <- data.frame(list(expanded = mlg, original = mlg, 
                  contracted = mlg, custom = factor(mlg)))
    } else if (class(mlg) %in% "MLG"){
      return(mlg)  
    }
    slot(.Object, "cutoff")   <- c(expanded = 0.0, contracted = 0.0) 
    slot(.Object, "mlg")      <- mlg
    slot(.Object, "visible")  <- "original"
    slot(.Object, "distname") <- "nei.dist"
    slot(.Object, "distargs") <- list()
    slot(.Object, "distalgo") <- "farthest_neighbor"
    return(.Object)
  }
)

#==============================================================================#
#' @rdname MLG-method
#' @param x an MLG object
#' @param i a vector of integers or logical values to index the MLG vector.
#' @param all a logical value indicating whether or not to return the subset of
#' all MLG values or only the numeric.
#==============================================================================#
setMethod(
  f = "[",
  signature = "MLG",
  definition = function(x, i, j, ..., all = FALSE, drop = TRUE){
    if (missing(i)) i <- TRUE

    if (missing(j) & all){ # Retain the state of the MLG object, 
      x@mlg <- x@mlg[i, ]  # but return the subset rows.
      return(x)
    } else if (missing(j)){
      j <- x@visible
    }
    
    return(x@mlg[i, j, ..., drop = drop])
  }
)

#==============================================================================#
#' @rdname MLG-method
#' @param j One of 
#' \itemize{
#'  \item "original"
#'  \item "expanded"
#'  \item "contracted"
#'  \item "custom"
#' }
#' @param value the vector of MLGs to replace. \strong{For internal use only}.
#==============================================================================#
setReplaceMethod(
  f = "[",
  signature = "MLG",
  definition = function(x, i, j, value){
    if (missing(j)) j <- x@visible
    if (missing(i)) i <- TRUE
    x@mlg[i, j] <- value
    return(x)
  }
)

#==============================================================================#
#' @rdname MLG-method
#' @param object an MLG object.
#==============================================================================#
setMethod(
  f = "show",
  signature("MLG"),
  definition = function(object){
    if (nrow(object@mlg) == 0){
      print("an empty MLG object")
      return()
    }
    cutoff <- ifelse(!object@visible %in% c("original", "custom"), 
                      paste(" mlgs with a cutoff of", object@cutoff[object@visible]),
                      " mlgs")
    dist <- ifelse(object@visible != "contracted", 
                   ".", 
                   paste(" based on the function", object@distname))
    cat(length(object), " ", object@visible, cutoff, dist, "\n", sep = "")
    print(object@mlg[[object@visible]])
  }
)
#==============================================================================#
#' @rdname MLG-method
#==============================================================================#
setMethod(
  f = "length",
  signature("MLG"),
  definition = function(x){
    callGeneric(x@mlg[[x@visible]])
  })




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# MATH METHODS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#==============================================================================#
#' @rdname MLG-method
#' @param e1 an MLG object
#' @param e2 a number
#==============================================================================#
setMethod(
  f = "Ops", 
  signature = "MLG",
  definition = function(e1, e2){
    callGeneric(e1@mlg[[e1@visible]], e2)
  }
)
#==============================================================================#
#' @rdname MLG-method
#==============================================================================#
setMethod(
  f = "Math", 
  signature = "MLG",
  definition = function(x){
    callGeneric(x@mlg[[x@visible]])
  }
)
#==============================================================================#
#' @rdname MLG-method
#' @param digits the number of digits to retain
#==============================================================================#
setMethod(
  f = "Math2", 
  signature = "MLG",
  definition = function(x, digits){
    callGeneric(x@mlg[[x@visible]], digits)
  }
)
#==============================================================================#
#' @rdname MLG-method
#' @param ... passed on to summary methods
#' @param na.rm passed on to summary methods
#==============================================================================#
setMethod(
  f = "Summary", 
  signature = "MLG",
  definition = function(x, ..., na.rm = FALSE){
    callGeneric(x@mlg[[x@visible]], ..., na.rm)
  }
)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SHOULD BE BUT AREN'T GENERIC METHODS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setGeneric("unique")
#==============================================================================#
#' Unique and Duplicated implementations for MLG objects
#' 
#' internal use
#' 
#' @rdname unique-methods
#' @aliases unique,MLG-method
#' @param incomparables vector of values that cannot be compared
#' @param ... options passed on to the base function \link{unique} or \link{duplicated}.
#' @export
#' @keywords internal
#' @docType methods
#==============================================================================#
setMethod(
  f = "unique",
  signature("MLG"),
  definition = function(x, incomparables = FALSE, ...){
    unique(x@mlg[[x@visible]], incomparables, ...)
  }
)

setGeneric("duplicated")
#==============================================================================#
#' @rdname unique-methods
#' @aliases duplicated,MLG-method
#' @export
#' @docType methods
#==============================================================================#
setMethod(
  f = "duplicated",
  signature("MLG"),
  definition = function(x, incomparables = FALSE, ...){
    duplicated(x@mlg[[x@visible]], incomparables, ...)
  } 
)

setGeneric("levels")

#==============================================================================#
#' Unique and Duplicated implementations for MLG objects
#' 
#' internal use
#' 
#' @rdname levels-methods
#' @aliases levels,MLG-method
#' @param x an MLG object
#' @return a character vector showing the levels of custom MLGs or NULL if the
#' visible slot is not set to "custom"
#' @export
#' @keywords internal
#' @docType methods
#==============================================================================#
setMethod(
  f = "levels",
  signature("MLG"),
  definition = function(x){
    return(levels(x[]))  
  }
)

setGeneric("levels<-")

#==============================================================================#
#' @rdname levels-methods
#' @aliases levels<-,MLG-method
#' @export
#' @docType methods
#==============================================================================#
setMethod(
  f = "levels<-",
  signature("MLG"),
  definition = function(x, value){
    if (x@visible != "custom"){
      warning("Cannot assign levels unless you have custom MLGs.", .immediate = TRUE)
    } else {
      levels(x@mlg[[x@visible]]) <- value
    }
    return(x)
  }
)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ACCESSORS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#==============================================================================#
#' Accessors for the MLG object
#' 
#' This documentation is for future developers of poppr. The accessors here are 
#' preferred over accessing the elements via the @@ symbol.
#' 
#' @param x an MLG object
#' @param value see details
#' @return see details
#' @details These accessors are intended for internal use only. They only affect
#'   MLG objects, not genind objects. Only visible and MLG2df are general for
#'   all forms of MLG. The distargs and cutoff are specific for use in
#'   mlg.filter or any function that offers filtering as an option. The argument
#'   "value" will always take the type defined in the \code{\linkS4class{MLG}}
#'   class.
#'   
#' @rdname MLG-accessors
#' @aliases visible,MLG-method
#' @export
#' @keywords internal
#' @examples
#' 
#' \dontrun{
#' # These examples will simply show you what you can do with these
#' set.seed(5000)
#' (x <- sample(10, 20, replace = TRUE))
#' (m <- new("MLG", x))
#' 
#' # Visibility ------------------------------
#' visible(m) # original
#' visible(m) <- "contracted"
#' m          # shows contracted MLGS
#' 
#' # Conversion to data frame ----------------
#' MLG2df(m)  # Grab the internal data frame
#' 
#' # Distance function handling --------------
#' distname(m) # nei.dist
#' distargs(m) # list()
#' distalgo(m) # farthest
#' cutoff(m)
#' 
#' distname(m) <- substitute("diss.dist")
#' distargs(m) <- list(percent = TRUE)
#' distalgo(m) <- "average"
#' cutoff(m)["contracted"] <- 0.2
#' 
#' }
#==============================================================================#
setGeneric("visible", function(x) {
  standardGeneric("visible")
})

setMethod(
  f = "visible", 
  signature(x = "MLG"), 
  definition = function(x) {
    return(x@visible)
})

#==============================================================================#
#' @rdname MLG-accessors
#' @aliases visible<-,MLG-method
#' @export
#==============================================================================#
setGeneric("visible<-", function(x, value) {
  standardGeneric("visible<-")
})

setMethod(
  f = "visible<-", 
  signature(x = "MLG"), 
  definition = function(x, value) {
    x@visible <- value
    return(x)
})

#==============================================================================#
#' @rdname MLG-accessors
#' @aliases MLG2df,MLG-method
#' @export
#==============================================================================#
setGeneric("MLG2df", function(x) {
  standardGeneric("MLG2df")
})

setMethod(
  f = "MLG2df", 
  signature(x = "MLG"), 
  definition = function(x) {
    x@mlg
})


#==============================================================================#
#' @rdname MLG-accessors
#' @aliases distname,MLG-method
#' @export
#==============================================================================#
setGeneric("distname", function(x) {
  standardGeneric("distname")
})

setMethod(
  f = "distname", 
  signature(x = "MLG"), 
  definition = function(x) {
    x@distname
  })

#==============================================================================#
#' @rdname MLG-accessors
#' @aliases distname<-,MLG-method
#' @export
#==============================================================================#
setGeneric("distname<-", function(x, value) {
  standardGeneric("distname<-")
})

setMethod(
  f = "distname<-", 
  signature(x = "MLG"), 
  definition = function(x, value) {
    x@distname <- value
    return(x)
  })

#==============================================================================#
#' @rdname MLG-accessors
#' @aliases distargs,MLG-method
#' @export
#==============================================================================#
setGeneric("distargs", function(x) {
  standardGeneric("distargs")
})

setMethod(
  f = "distargs", 
  signature(x = "MLG"), 
  definition = function(x) {
    x@distargs
})

#==============================================================================#
#' @rdname MLG-accessors
#' @aliases distargs<-,MLG-method
#' @export
#==============================================================================#
setGeneric("distargs<-", function(x, value) {
  standardGeneric("distargs<-")
})

setMethod(
  f = "distargs<-", 
  signature(x = "MLG"), 
  definition = function(x, value) {
    x@distargs <- value
    return(x)
  })


#==============================================================================#
#' @rdname MLG-accessors
#' @aliases distalgo,MLG-method
#' @export
#==============================================================================#
setGeneric("distalgo", function(x) {
  standardGeneric("distalgo")
})

setMethod(
  f = "distalgo", 
  signature(x = "MLG"), 
  definition = function(x) {
    x@distalgo
  })

#==============================================================================#
#' @rdname MLG-accessors
#' @aliases distalgo<-,MLG-method
#' @export
#==============================================================================#
setGeneric("distalgo<-", function(x, value) {
  standardGeneric("distalgo<-")
})

setMethod(
  f = "distalgo<-", 
  signature(x = "MLG"), 
  definition = function(x, value) {
    x@distalgo <- value
    return(x)
  })


#==============================================================================#
#' @rdname MLG-accessors
#' @aliases cutoff,MLG-method
#' @export
#==============================================================================#
setGeneric("cutoff", function(x) {
  standardGeneric("cutoff")
})

setMethod(
  f = "cutoff", 
  signature(x = "MLG"), 
  definition = function(x) {
    x@cutoff
  })

#==============================================================================#
#' @rdname MLG-accessors
#' @aliases cutoff<-,MLG-method
#' @export
#==============================================================================#
setGeneric("cutoff<-", function(x, value) {
  standardGeneric("cutoff<-")
})

setMethod(
  f = "cutoff<-", 
  signature(x = "MLG"), 
  definition = function(x, value) {
    x@cutoff <- value
    return(x)
  })