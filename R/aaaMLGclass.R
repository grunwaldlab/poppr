#' MLG class
#' 
#' A class to store multilocus genotypes in genclone objects. 
#' 
#' @name MLG-class
#' @rdname MLG-class
#' @aliases MLG
#' @export
#' @slot mlg a list containing four vectors, one for each type of MLG manipulation.
#' @slot visible a character specifying which MLG type is to be displayed and accessed.
#' @slot dist a function to calculate the genetic distance for collapsing MLGs.
#' @slot distname the name of the distance function
#' @slot Two numbers specifying the cutoff value for expanding and collapsing
#' MLGs.
#' @author Zhian N. Kamvar
#' @seealso \code{\link{as.genclone}} \code{\link{sethierarchy}} \code{\link{setpop}} 
#' \code{\linkS4class{genind}} 
setClass("MLG", 
         representation(visible = "character",
                        cutoff = "numeric",
                        dist = "function",
                        distname = "character",
                        mlg = "list"))