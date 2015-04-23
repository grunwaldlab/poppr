#' MLG class
#' 
#' A class to store multilocus genotypes in genclone objects. 
#' 
#' @name MLG-class
#' @rdname MLG-class
#' @export
#' @slot mlg a list containing four vectors, one for each type of MLG manipulation.
#' @slot visible a character specifying which MLG type is to be displayed and accessed.
#' @slot distname the name of the distance function used to collapse mlgs.
#' @slot distargs the arguments provided to compute the distance function.
#' @slot cutoff Two numbers specifying the cutoff value for expanding and collapsing
#' MLGs.
#' @author Zhian N. Kamvar
#' @seealso \code{\linkS4class{genclone}} \code{\link{mll}} \code{\linkS4class{snpclone}}
#' @keywords internal
setClassUnion("charORLang", c("character", "language"))
setClass("MLG", 
         representation(visible = "character",
                        cutoff = "numeric",
                        distname = "charORLang",
                        distargs = "list",
                        mlg = "data.frame"))