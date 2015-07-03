



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
#' Default methods for acessing and subsetting MLG objects.
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









