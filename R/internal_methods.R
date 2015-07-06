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
# The point of this file serves as a way of parsing the internal functions that
# are specific to methods since there are now two classes that share common 
# slots. Having the methods specific to these internal slots ensures that the 
# code is a) DRY and b) debuggable (by calling the internal function).
#==============================================================================#
# gethierarchy.internal <- function(x, formula = NULL, combine = TRUE){
# 	if (is.null(formula)) return(x@hierarchy)
# 	vars <- all.vars(formula)
# 	if (any(!vars %in% names(x@hierarchy))){
# 	  stop(hier_incompatible_warning(vars, x@hierarchy))
# 	}
# 	if (combine){
# 	  hier <- make_hierarchy(formula, x@hierarchy)
# 	} else {
# 	  hier <- x@hierarchy[all.vars(formula)]
# 	}
# 	invisible(return(hier))
# }

# sethierarchy.internal <- function(x, value){
#   if (!inherits(value, "data.frame")){
#     stop(paste(substitute(value), "is not a data frame"))
#   }
#   if (nrow(value) != nInd(x)){
#     stop("Number of rows in data frame not equal to number of individuals in object.")
#   }
#   value <- data.frame(lapply(value, function(f) factor(f, unique(f))))
#   x@hierarchy <- value
#   return(x)
# }

# namehierarchy.internal <- function(x, value){
#   if (is.language(value)){
#     value <- all.vars(value)
#   }
#   if (!is.vector(value) | length(value) != length(x@hierarchy)){
#     stop(paste("Hierarchy, needs a vector argument of length", length(x@hierarchy)))
#   }
#   names(x@hierarchy) <- value
#   return(x)
# }

# splithierarchy.internal <- function(x, value, sep = "_"){
#   if (is.language(value)){
#     # valterms <- attr(terms(value), "term.labels")
#     # valterms <- valterms[length(valterms)]
#     # valterms <- gsub(":", sep, valterms)
#     value    <- all.vars(value)
#   } else {
#     stop("value must be a formula.")
#   }
#   if (length(value) < 1){
#     stop("value must have more than one hierarchical level.")
#   }
#   hierarchy  <- x@hierarchy
#   if (length(hierarchy) > 1){
#     warning("Hierarchy must be length 1. Taking the first column.")
#     hierarchy <- hierarchy[1]
#   }
#   seps     <- gregexpr(sep, hierarchy[[1]])
#   sepmatch <- vapply(seps, function(val) all(as.integer(val) > 0), logical(1))
#   seps     <- vapply(seps, length, numeric(1))
#   all_seps_match <- all(sepmatch)
#   given_seps     <- length(value) - 1 
#   if (!all_seps_match | all(seps != given_seps)){
#     seps <- ifelse(all_seps_match, seps[1], 0) + 1
#     msg1 <- paste("\n  Data has", seps, ifelse(seps == 1, "level", "levels"),
#                   "of hierarchy with the separator", sep, ".")
#     msg2 <- paste("Here is the fist column of the data:", hierarchy[1, ])
#     stop(paste(msg1, "\n ", msg2))
#   }
#   x@hierarchy <- reshape2::colsplit(as.character(hierarchy[[1]]), pattern = sep, value)
#   x@hierarchy <- data.frame(lapply(x@hierarchy, function(f) factor(f, levels = unique(f))))
#   # names(hierarchy) <- value
#   # x@hierarchy      <- hierarchy
#   return(x) 
# }

# addhierarchy.internal <- function(x, value, name = "NEW"){    
#   hierarchy  <- x@hierarchy
#   if ((is.vector(value) | is.factor(value)) & length(value) == nrow(hierarchy)){
#     value <- factor(value, levels = unique(value))
#     NEW <- data.frame(value)
#     names(NEW) <- name
#     hierarchy <- cbind(hierarchy, NEW)
#   } else if (is.data.frame(value) && nrow(value) == nrow(hierarchy)){
#     value <- data.frame(lapply(value, function(f) factor(f, unique(f))))
#     hierarchy <- cbind(hierarchy, value)
#   } else {
#     stop("value must be a vector or data frame.")
#   }
#   x@hierarchy <- hierarchy
#   return(x) 
# }

# setpop.internal <- function(x, formula = NULL){
#   if (is.null(formula) | !is.language(formula)){
#     stop(paste(substitute(formula), "must be a valid formula object."))
#   }
#   vars <- all.vars(formula)
#   if (!all(vars %in% names(x@hierarchy))){
#     stop(hier_incompatible_warning(vars, x@hierarchy))
#   }
#   pop(x) <- make_hierarchy(formula, x@hierarchy)[[length(vars)]]
#   return(x)
# }
#==============================================================================#
# Internal functions to deal with the MLG class.
#==============================================================================#

mll.reset.internal <- function(x, value){
  if (!is(x@mlg, "MLG")){
    x@mlg <- new("MLG", x@mlg)
    return(x)
  }
  if (is.logical(value) && length(value) == 1 && value == TRUE){
    x@mlg <- new("MLG", mlg.vector(x, reset = TRUE))
    return(x)
  }
  new_mlg <- x@mlg@mlg
  TYPES   <- c("original", "expanded", "contracted", "custom")
  types   <- match.arg(value, TYPES, several.ok = TRUE)
  
  if ("original" %in% types){
    newvec <- mlg.vector(x, reset = TRUE)
  } else {
    newvec <- mll(x, "original")
  }
  for (i in types){
    if (i == "custom"){
      new_mlg[i] <- as.factor(newvec)
    } else {
      new_mlg[i] <- newvec
    }
  }
  x@mlg@mlg <- new_mlg
  if ("contracted" %in% types){
    x@mlg@cutoff["contracted"] <- 0
    x@mlg@distname             <- "nei.dist"
    x@mlg@distargs             <- list()
  }
  return(x)
}


mll.custom.internal <- function(x, set = TRUE, value){
  if (!is(x@mlg, "MLG")){
    x@mlg <- new("MLG", x@mlg)
  }
  if (missing(value)){
    return(x@mlg[, "custom"])
  }
  mlgs <- x@mlg
  if (length(value) != length(mlgs)){
    stop("value must be the same length as the mlls")
  }
  if (!is.factor(value)){
    value <- factor(value)
  }
  if (set){
    mlgs@visible <- "custom"
  }
  mlgs@mlg[, "custom"] <- value
  x@mlg <- mlgs
  return(x)
}

mll.levels.internal <- function(x, set = TRUE, value){
  if (!is(x@mlg, "MLG")){
    x@mlg <- new("MLG", x)
  }
  mlgs <- x@mlg
  if (missing(value)){
    return(levels(mlgs))
  }
  if (length(value) != nlevels(mlgs@mlg[, "custom"])){
    stop("value length should match the number of values in mlg")
  }
  if (set){
    vis <- "custom"
  } else {
    vis <- mlgs@visible
  }
  mlgs@visible <- "custom"
  levels(mlgs) <- value
  x@mlg        <- mlgs
  mll(x)       <- vis
  return(x)
}
