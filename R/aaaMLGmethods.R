



#==============================================================================#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Below are methods for the MLG class. I will break up the sections into generic
# methods that are needed for them to make sense in the R world, custom generic
# methods that aren't generic in R, but need to be to make these work, and
# finally a new method specifically for this object.
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
    if (is.vector(mlg)){
      mlg <- list(expanded = mlg, original = mlg, 
                  contracted = mlg, custom = rep(NA, length(mlg)))
    }
    slot(.Object, "cutoff")   <- c(expanded = 0.0, contracted = 0.0) 
    slot(.Object, "mlg")      <- mlg
    slot(.Object, "visible")  <- "original"
    slot(.Object, "dist")     <- nei.dist
    slot(.Object, "distname") <- "nei.dist"
    return(.Object)
  }
)

#==============================================================================#
#' @rdname MLG-method
#' @param x an MLG object
#' @param i a vector of integers or logical values to index the MLG vector.
#==============================================================================#
setMethod(
  f = "[",
  signature = "MLG",
  definition = function(x, i){
    return(x@mlg[[x@visible]][i])
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
    
    if (j != "custom" & !is.numeric(value)){
      stop("can't set this manually")
    } else if (j == "custom"){
      if (is.null(i)){
        i <- 1:length(x@mlg[["original"]])
      }
      x@mlg[[j]][i] <- value
    } else {
      x@mlg[[j]] <- value
    }
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





#' @export
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
unique.MLG <- function(x, incomparables = FALSE, ...){
  unique(x@mlg[[x@visible]], incomparables, ...)
}

setMethod(
  f = "unique",
  signature("MLG"),
  definition = unique.MLG
)
#==============================================================================#
#' @rdname unique-methods
#' @aliases duplicated,MLG-method
#' @export
#' @docType methods
#==============================================================================#
duplicated.MLG <- function(x, incomparables = FALSE, ...){
  duplicated(x@mlg[[x@visible]], incomparables, ...)
}

#' @export
setGeneric("duplicated")

setMethod(
  f = "duplicated",
  signature("MLG"),
  definition = duplicated.MLG
)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# NEW METHODS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




#==============================================================================#
#' Access and manipluate multilocus genotypes or lineages
#' 
#' The MLG object in a \linkS4class{genclone} object contains 4 vectors to give
#' the user flexibility of MLG assignment. 
#' 
#' @rdname mll-methods
#' @aliases setmlg,MLG-method
#' @param x an MLG object
#' @param value one of the 4 possible MLG types.
#' 
#' @details
#' The four possible MLG types:
#' 
#' \itemize{
#' \item "original" - Strict definition of MLGs. Sensitive to small changes in 
#' data and missing data.
#' \item "contracted" - aka Multilocus Lineages. These MLGs are collapsed into 
#' groups given a certain genetic distance, cutoff value, and algorithm
#' \item "expanded" - Expanded multilocus genotype groups for when there are too
#' few loci. Psex is caluclated and mlgs are split into new MLGs based on a given
#' pvalue cutoff of Psex.
#' \item "custom" - User defined MLGs
#' }
#' 
#' @export
#' @docType methods
#==============================================================================#
setmlg <- function(x, value = NULL){
  standardGeneric("setmlg")
} 

#' @export
setGeneric("setmlg")

setMethod(
  f = "setmlg",
  signature(x = "MLG"),
  definition = function(x, value){
    if (!is.null(value)){
      ARGS  <- c("expanded", "original", "contracted", "custom")
      value <- match.arg(value, ARGS)
      x@visible <- value
    }
    return(x)
  }
)

#==============================================================================#
#' @rdname mll-methods
#' @aliases setmlg<-,MLG-method
#' @export
#' @docType methods
#==============================================================================#
"setmlg<-" <- function(x, value = NULL){
  standardGeneric("setmlg<-")
} 

#' @export
setGeneric("setmlg<-")

setMethod(
  f = "setmlg<-",
  signature(x = "MLG"),
  definition = function(x, value){
    return(setmlg(x, value))
  }
)





