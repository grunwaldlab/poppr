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
# This file will contain deprecated functions, namely hierarchy methods. They
# are here for backwards compatibility with poppr 1.1, but will strongly warn
# the user against them.
#==============================================================================#
#' DEPRECATED FUNCTIONS
#' 
#' These functions are deprecated as of \pkg{poppr} version 2.0 and have been 
#' moved to the \pkg{adegenet} package. See \code{\link[adegenet]{strata}} for 
#' details.
#' 
#' @param ... arguments to be passed on to adegenet 
#'   \code{\link[adegenet]{strata}} and \code{\link[adegenet]{setPop}} methods.
#' @param value argument passed on to adegenet methods.
#'   
#' @return a deprecation message and the same return as the adegenet methods.
#'   
#' @section Background:
#'   These functions originally appeared in \pkg{poppr} in July of 2014. In
#'   March 2015, a request was made for them to be ported to the \pkg{adegenet}
#'   package as part of the population genetics in R hackathon. The name of the
#'   functions were moved from "hierarchy" to "strata" because it better
#'   reflected the nature of the fields to contain any information possible.
#'   
#' @seealso \code{\link[adegenet]{strata}} \code{\link[adegenet]{setPop}}
#' @rdname deprecated
#' @export
#==============================================================================#
gethierarchy <- function(...){
	the_call <- match.call()
	the_call[[1]] <- quote(strata)
	dmsg <- paste("'gethierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'strata'.\n\n Please use:\n", 
								utils::capture.output(print(the_call)), "\n")
	.Deprecated("strata", package = "adegenet", dmsg)
	return(strata(...))
}

#' @rdname deprecated
#' @export
sethierarchy <- function(...){
	the_call <- match.call()
	the_call[[1]] <- quote(strata)
	# since sethierarchy and gethierarchy have merged, the arguments have also
	# merged. Since sethierarchy originally had only two arguments, this will make
	# sure that the arguments work. 
	names(the_call) <- c(NA, NA, "value")
	dots <- list(...)
	dmsg <- paste("'sethierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'strata'.\n\n Please use:\n", 
								utils::capture.output(print(the_call)), "\n")
	.Deprecated("strata", package = "adegenet", dmsg)
	return(strata(dots[[1]], value = dots[[2]]))
}

#' @rdname deprecated
#' @export
setpop <- function(...){
  the_call <- match.call()
  the_call[[1]] <- quote(setPop)
  dmsg <- paste("'setpop' has been deprecated, moved to the adegenet package, and",
                "renamed to 'setPop'.\n\n Please use:\n", 
                utils::capture.output(print(the_call)), "\n")
  .Deprecated("setPop", package = "adegenet", dmsg)
  return(setPop(...))
}

#' @rdname deprecated
#' @export
"setpop<-" <- function(..., value){
  the_call <- match.call()
  the_call[[1]] <- quote(setPop)
  dmsg <- paste("'setpop' has been deprecated, moved to the adegenet package, and",
                "renamed to 'setPop'.\n\n Please use:\n", 
                "setPop()", "\n")
  .Deprecated("setPop", package = "adegenet", dmsg)
  return(setPop(..., value))
}


#' @rdname deprecated
#' @export
"sethierarchy<-" <- function(..., value){
	the_call <- match.call()
	the_call[[1]] <- quote(strata)
	dmsg <- paste("'sethierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'strata'.\n\n Please use:\n", 
								"strata()", "\n")
	.Deprecated("strata", package = "adegenet", dmsg)
	return(strata(..., value = value))
}


#' @rdname deprecated
#' @export
addhierarchy <- function(...){
	the_call <- match.call()
	the_call[[1]] <- quote(addStrata)
	dmsg <- paste("'addhierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'addStrata'.\n\n Please use:\n", 
								utils::capture.output(print(the_call)), "\n")
	.Deprecated("addStrata", package = "adegenet", dmsg)
	return(addStrata(...))
}

#' @rdname deprecated
#' @export
"addhierarchy<-" <- function(..., value){
	the_call <- match.call()
	the_call[[1]] <- quote(addStrata)
	dmsg <- paste("'addhierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'addStrata'.\n\n Please use:\n", 
								"addStrata()", "\n")
	.Deprecated("addStrata", package = "adegenet", dmsg)
	return(addStrata(..., value = value))
}

#' @rdname deprecated
#' @export
splithierarchy <- function(...){
	the_call <- match.call()
	the_call[[1]] <- quote(splitStrata)
	dmsg <- paste("'splithierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'splitStrata'.\n\n Please use:\n", 
								utils::capture.output(print(the_call)), "\n")
	.Deprecated("splitStrata", package = "adegenet", dmsg)
	return(splitStrata(...))
}

#' @rdname deprecated
#' @export
"splithierarchy<-" <- function(..., value){
	the_call <- match.call()
	the_call[[1]] <- quote(splitStrata)
	dmsg <- paste("'splithierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'splitStrata'.\n\n Please use:\n", 
								"splitStrata()", "\n")
	.Deprecated("splitStrata", package = "adegenet", dmsg)
	return(splitStrata(..., value))
}

#' @rdname deprecated
#' @export
namehierarchy <- function(...){
	the_call <- match.call()
	the_call[[1]] <- quote(nameStrata)
	dmsg <- paste("'namehierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'nameStrata'.\n\n Please use:\n", 
								utils::capture.output(print(the_call)), "\n")
	.Deprecated("nameStrata", package = "adegenet", dmsg)
	return(nameStrata(...))
}

#' @rdname deprecated
#' @export
"namehierarchy<-" <- function(..., value){
	the_call <- match.call()
	the_call[[1]] <- quote(nameStrata)
	dmsg <- paste("'namehierarchy' has been deprecated, moved to the adegenet package, and",
								"renamed to 'nameStrata'.\n\n Please use:\n", 
								"nameStrata()", "\n")
	.Deprecated("nameStrata", package = "adegenet", dmsg)
	return(nameStrata(..., value))
}