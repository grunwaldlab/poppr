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
setClassUnion("charORLang", c("character", "language"))
#==============================================================================#
#' MLG class
#' 
#' A class to store multilocus genotypes in genclone objects. This is intended
#' for internal use only.
#' 
#' @name MLG-class
#' @rdname MLG-class
#' @aliases MLG
#' @export
#' @slot mlg a list containing four vectors, one for each type of MLG 
#'   manipulation.
#' @slot visible a character specifying which MLG type is to be displayed and 
#'   accessed.
#' @slot distname the name of the distance function or matrix used to collapse 
#'   mlgs.
#' @slot distenv the environment that contains the distance function or matrix
#' @slot distargs the arguments provided to compute the distance function.
#' @slot distalgo the algorithm used to contract multilocus genotypes.
#' @slot cutoff Two numbers specifying the cutoff value for expanding and 
#'   collapsing MLGs.
#' @author Zhian N. Kamvar
#' @seealso \code{\linkS4class{genclone}} \code{\linkS4class{snpclone}}
#'   \code{\link{mll}} For developers: \code{\link{visible}}
#' @keywords internal
#' @examples
#' 
#' # These examples will simply show you what you can do with these
#' set.seed(5000)
#' (x <- sample(10, 20, replace = TRUE))
#' (m <- new("MLG", x))
#' 
#'  visible(m) # original is always default
#'  
#'  m[]       # adding braces after the object will always return a vector of 
#'            # the same type as defined in "visible"
#'            
#'  m + 1     # You can do math on the numeric ones
#'  
#'  visible(m) <- "custom"
#'  m + 2     # This should throw a warning
#'  # The types are stored in a data frame. You can retrieve them easily:
#'  visible(m) <- "original"
#'  m
#'  m[, "custom"]
#'  
#'  # Important for subsetting, if you subset the object, normally, it will 
#'  # return a vector unless you specify all = TRUE
#'  m[1:10]             # original vector
#'  m[1:10, all = TRUE] # still class MLG
#' 
#==============================================================================#
setClass("MLG", 
         representation(visible = "character",
                        cutoff = "numeric",
                        distname = "charORLang",
                        distenv  = "environment",
                        distargs = "list",
                        distalgo = "character",
                        mlg = "data.frame"),
         prototype(visible = character(0), 
                   cutoff = numeric(0), 
                   distname = character(0), 
                   distenv = as.environment(.GlobalEnv),
                   distargs = list(),
                   distalgo = "farthest_neighbor", 
                   mlg = data.frame(expanded = numeric(0), 
                                    original = numeric(0), 
                                    contracted = numeric(0), 
                                    custom = factor(character(0))
                                    )
                   )
         )
