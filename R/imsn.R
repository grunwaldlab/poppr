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
#==============================================================================#
#' Create minimum spanning networks interactively
#' 
#' This function will launch an interactive interface that allows you to create,
#' plot, manipulate, and save minimum spanning networks. It runs using the
#' \pkg{shiny} R package.
#' 
#' @return NULL, invisibly
#' 
#' @details Creating and plotting MSNs requires three steps:
#' 
#' \enumerate{
#' \item Create a distance matrix from your data
#' \item Create a minimum spanning network with your data and the matrix
#' \item Visualize the minimum spanning network}
#' The function \code{\link{plot_poppr_msn}} is currently the most flexible way
#' of visualizing your minimum spanning network, but with 20 parameters, it can
#' become pretty intimidating trying to find the right display for your MSN. 
#' 
#' With this function, all three steps are combined into one interactive 
#' interface that will allow you to intuitively modify your minimum spanning 
#' network and even save the results to a pdf or png file. 
#' 
#' @section Interface:
#' \subsection{Buttons}{
#' In the left hand panel, there are three buttons to execute the functions.
#' These allow you to run the data set after you manipulate all of the
#' parameters.
#' 
#' \itemize{
#' 
#' \item \strong{GO!} - This button will start the application with the
#' specified parameters
#' \item \strong{reData} - Use this button when you have changed any parameters
#' under the section \strong{Data Parameters}. This involves recalculating the
#' distance matrix and msn.
#' \item \strong{reGraph} - Use this button when you have changed any parameters
#' under the section \strong{Graphical Parameters}. This involves superficial
#' changes to the display of the minimum spanning network.
#' 
#' }
#' }
#' 
#' \subsection{Tabs}{
#' 
#' The right hand panel contains different tabs related to your data set of
#' choice.
#' 
#' \itemize{
#' 
#' \item \strong{Plot} - The minimum spanning network itself
#' \item \strong{Data} - A display of your data set
#' \item \strong{Command} - The commands used to create the plot. You can copy
#' and paste this to an R file for reproducibility.
#' \item \strong{Save Plot} - This provides a tool for you to save the plot to a
#' PDF or PNG image.
#' \item \strong{Session Information} - displays the result of
#' \code{\link{sessionInfo}} for reproducibility.
#' 
#' }
#' }
#' 
#' @author Zhian N. Kamvar
#' 
#' @seealso \code{\link{plot_poppr_msn}} \code{\link{diss.dist}}
#'   \code{\link{bruvo.dist}} \code{\link{bruvo.msn}} \code{\link{poppr.msn}}
#'   \code{\link{nei.dist}} \code{\link{popsub}} \code{\link{missingno}}
#'   
#' @export
#' @examples 
#' \dontrun{
#' 
#' # Set up some data
#' library("poppr")
#' library("magrittr")
#' data(monpop)
#' splitStrata(monpop) <- ~Tree/Year/Symptom
#' summary(monpop)
#' monpop_ssr <- c(CHMFc4 = 7, CHMFc5 = 2, CHMFc12 = 4, 
#'                 SEA = 4, SED = 4, SEE = 2, SEG = 6, 
#'                 SEI = 3, SEL = 4, SEN = 2, SEP = 4, 
#'                 SEQ = 2, SER = 4)
#' t26 <- monpop %>% setPop(~Tree) %>% popsub("26") %>% setPop(~Year/Symptom)
#' t26
#' if (interactive()) {
#'   imsn() # select Bruvo's distance and enter "monpop_ssr" into the Repeat Length field.
#'   
#'   # It is also possible to run this from github if you are connected to the internet.
#'   # This allows you to access any bug fixes that may have been updated before a formal
#'   # release on CRAN
#' 
#'   shiny::runGitHub("grunwaldlab/poppr", subdir = "inst/shiny/msn_explorer")
#' 
#'   # You can also use your own distance matrices, but there's a small catch.
#'   # in order to do so, you must write a function that will subset the matrix
#'   # to whatever populations are in your data. Here's an example with the above
#  ' # data set:
#' 
#'   mondist <- bruvo.dist(monpop, replen = monpop_ssr)
#'   myDist <- function(x, d = mondist){
#'    dm <- as.matrix(d)          # Convert the dist object to a square matrix
#'    xi <- indNames(x)           # Grab the sample names that exist
#'    return(as.dist(dm[xi, xi])) # return only the elements that have the names
#'                                # in the data set
#'   }
#'   # After executing imsn, choose:
#'   # Distance: custom
#'   # myDist
#'   imsn() 
#' }
#' }
#' @import shiny
imsn <- function(){
  shiny::runApp(system.file("shiny", "msn_explorer", package = "poppr"))
  invisible(NULL)
}
