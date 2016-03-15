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
#' The \pkg{poppr} R package.
#' 
#' @description \pkg{Poppr} provides tools for population genetic analysis that
#' include genotypic diversity measures, genetic distances with bootstrap
#' support, native organization and handling of population hierarchies, and
#' clone correction.
#' 
#' To cite \pkg{poppr}, please use \code{citation("poppr")}. When referring to
#' \pkg{poppr} in your manuscript, please use lower case unless it occurs at the
#' beginning of a sentence.
#' 
#' @details This package relies on the \pkg{\link[adegenet]{adegenet}} package. 
#'   It is built around the \code{\linkS4class{genind}} and 
#'   \code{\linkS4class{genlight}} object. Genind objects store genetic
#'   information in a table of allele frequencies while genlight objects store
#'   SNP data efficiently by packing binary allele calls into single bits.
#'   \pkg{Poppr} has extended these object into new objects called 
#'   \code{\linkS4class{genclone}} and \code{\linkS4class{snpclone}},
#'   respectively. These objects are designed for analysis of clonal organisms
#'   as they add the \strong{@@mlg} slot for keeping track of multilocus
#'   genotypes and multilocus lineages.
#'   
#'   \subsection{Documentation}{ Documentation is available for any function by 
#'   typing \code{?function_name} in the R console. Detailed topic explanations
#'   live in the package vignettes:
#'   
#'   \tabular{ll}{
#'     \strong{Vignette}              \tab \strong{command}\cr
#'     Migration from poppr version 1 \tab \code{vignette("how_to_migrate", "poppr")}\cr
#'     Data import and manipulation   \tab \code{vignette("poppr_manual", "poppr")}\cr
#'     Algorithms and Equations      \tab \code{vignette("algo", "poppr")}\cr
#'     Multilocus Genotype Analysis   \tab \code{vignette("mlg", "poppr")}  
#'   }
#'   
#'   Important notes on the changes from poppr version 1 to poppr version 2 are
#'   in \emph{Migration from poppr version 1}, essential functions for importing
#'   and manipulating data are detailed within the \emph{Data import and 
#'   manipulation} vignette, details on algorithms used in \pkg{poppr} are
#'   within the \emph{Algorithms and equations} vignette, and details for
#'   working with multilocus genotypes are in \emph{Multilocus Genotype
#'   Analysis}. Examples of analyses are available in a primer written by 
#'   Niklaus J. Grünwald, Zhian N. Kamvar, and Sydney E. Everhart at 
#'   \url{http://grunwaldlab.github.io/Population_Genetics_in_R}.}
#'   
#'   \subsection{Getting help}{ If you have a specific question or issue with 
#'   \pkg{poppr}, feel free to contribute to the google group at 
#'   \url{https://groups.google.com/group/poppr}. If you find a bug and 
#'   are a github user, you can submit bug reports at 
#'   \url{https://github.com/grunwaldlab/poppr/issues}. Otherwise, leave a
#'   message on the groups. Personal emails are highly discouraged as they do 
#'   not allow others to learn.}
#'   
#' @section A note about version 2:
#'   \pkg{Poppr} has recently moved to version 2 partially because
#'   \pkg{adegenet} has also moved to version 2, and partially because we
#'   needed to change the way a few things were done. Major changes that have
#'   affected both packages: the main data structure in \pkg{adegenet}, the
#'   \code{\linkS4class{genind}} object, has changed drastically and reduced in
#'   size. You must convert your data if you want to use the new functionality.
#'   See the vignette "Migration from poppr version 1" by typing
#'   \code{vignette("how_to_migrate", "poppr")} for more details on
#'   major changes and how to deal with them.
#'   
#'   Many new features in \pkg{poppr} are described in Kamvar, Brooks, and
#'   Grünwald (2015), linked below.
#'   
#' @section Functions in \pkg{poppr}:
#' 
#'   Below are descriptions and links to functions found in \pkg{poppr}. Be
#'   aware that all functions in \pkg{\link[adegenet]{adegenet}} are also
#'   available. The functions are documented as:
#'   \itemize{
#'   \item \code{function} (data type) - Description
#'   }
#'   Where \sQuote{data type} refers to the type of data that can be used:
#'   \tabular{ll}{
#'     \strong{m} \tab a genclone or genind object \cr
#'     \strong{s} \tab a snpclone or genlight object \cr
#'     \strong{x} \tab a different data type (e.g. a matrix from \code{\link{mlg.table}})
#'   }
#'   
#' @section Data import/export:
#' \itemize{
#' \item \code{\link{getfile}} (x) - Provides a quick GUI to grab files for import
#' \item \code{\link{read.genalex}} (x) - Reads GenAlEx formatted csv files to a genind object
#' \item \code{\link{genind2genalex}} (m) - Converts genind objects to GenAlEx formatted csv files
#' \item \code{\link{genclone2genind}} (m) - Removes the @@mlg slot from genclone objects
#' }
#' @section Data manipulation:
#' \itemize{
#' \item \code{\link{as.genclone}} (m) - Converts genind objects to genclone objects
#' \item \code{\link{missingno}} (m) - Handles missing data
#' \item \code{\link{clonecorrect}} (m | s) - Clone-censors at a specified population hierarchy
#' \item \code{\link{informloci}} (m) - Detects and removes phylogenetically uninformative loci
#' \item \code{\link{popsub}} (m | s) - Subsets genind objects by population
#' \item \code{\link{shufflepop}} (m) - Shuffles genotypes at each locus using four different shuffling algorithms
#' \item \code{\link{recode_polyploids}} (m|x) - Recodes polyploid data sets with missing alleles imported as "0"
#' }
#' @section Genetic distances:
#' \itemize{
#' \item \code{\link{bruvo.dist}} (m) - Bruvo's distance (see also: \code{\link{fix_replen}})
#' \item \code{\link{diss.dist}} (m) - Absolute genetic distance (see prevosti.dist)
#' \item \code{\link{nei.dist}} (m | x) - Nei's 1978 genetic distance
#' \item \code{\link{rogers.dist}} (m | x) - Rogers' euclidean distance
#' \item \code{\link{reynolds.dist}} (m | x) - Reynolds' coancestry distance
#' \item \code{\link{edwards.dist}} (m | x) - Edwards' angular distance
#' \item \code{\link{prevosti.dist}} (m | x) - Prevosti's absolute genetic distance
#' \item \code{\link{bitwise.dist}} (s) - Calculates fast pairwise distances for genlight objects. 
#' }
#' @section Bootstrapping:
#' \itemize{
#' \item \code{\link{aboot}} (m | s | x) - Creates a bootstrapped dendrogram for any distance measure
#' \item \code{\link{bruvo.boot}} (m) - Produces dendrograms with bootstrap support based on Bruvo's distance
#' \item \code{\link{diversity_boot}} (x) - Generates boostrap distributions of diversity statistics for multilocus genotypes
#' \item \code{\link{diversity_ci}} (m | s | x) - Generates confidence intervals for multilocus genotype diversity.
#' }
#' @section Analysis:
#' \subsection{Multilocus Genotypes}{
#' \itemize{
#' \item \code{\link{mlg}} (m | s) - Calculates the number of multilocus genotypes
#' \item \code{\link{mll}} (m | s) - Displays the current multilocus lineages (genotypes) defined.
#' \item \code{\link{nmll}} (m | s) - Same as \code{\link{mlg}}.
#' \item \code{\link{mlg.crosspop}} (m | s) - Finds all multilocus genotypes that cross populations
#' \item \code{\link{mlg.table}} (m | s) - Returns a table of populations by multilocus genotypes
#' \item \code{\link{mlg.vector}} (m | s) - Returns a vector of a numeric multilocus genotype assignment for each individual
#' \item \code{\link{mlg.id}} (m | s) - Finds all individuals associated with a single multilocus genotype
#' \item \code{\link{mlg.filter}} (m | s) - Collapses MLGs by genetic distance
#' \item \code{\link{filter_stats}} (m | s) - Calculates mlg.filter for all algorithms and plots
#' \item \code{\link{cutoff_predictor}} (x) - Predicts cutoff threshold from mlg.filter. 
#' \item \code{\link{mll.custom}} (m | s) - Allows for the custom definition of multilocus lineages
#' \item \code{\link{mll.levels}} (m | s) - Allows the user to change levels of custom MLLs. 
#' \item \code{\link{mll.reset}} (m | s) - Reset multilocus lineages. 
#' \item \code{\link{diversity_stats}} (x) - Creates a table of diversity indices for multilocus genotypes. 
#' }
#' }
#' \subsection{Other}{
#' \itemize{
#' 
#' \item \code{\link{poppr.amova}} (m) - Analysis of Molecular Variance (as implemented in ade4)
#' \item \code{\link{ia}} (m) - Calculates the index of association
#' \item \code{\link{pair.ia}} (m) - Calculates the index of association for all loci pairs.
#' \item \code{\link{win.ia}} (s) - Index of association windows for genlight objects.
#' \item \code{\link{samp.ia}} (s) - Index of association on random subsets of loci for genlight objects.
#' \item \code{\link{poppr}} (m | x) - Returns a diversity table by population
#' \item \code{\link{poppr.all}} (m | x) - Returns a diversity table by population for all compatible files specified
#' \item \code{\link{private_alleles}} (m) - Tabulates the occurrences of alleles that only occur in one population.
#' \item \code{\link{locus_table}} (m) - Creates a table of summary statistics per locus.
#' \item \code{\link{rrmlg}} (m | x) - Round-robin multilocus genotype estimates.
#' \item \code{\link{rraf}} (m) - Round-robin allele frequency estimates.
#' \item \code{\link{pgen}} (m) - Probability of genotypes.
#' \item \code{\link{psex}} (m) - Probability of observing a genotype more than once.
#' }
#' }
#' @section Visualization:
#' \itemize{
#' \item \code{\link{imsn}} (m | s) - Interactive construction and visualization of minimum spanning networks
#' \item \code{\link{plot_poppr_msn}} (m | s | x) - Plots minimum spanning networks produced in poppr with scale bar and legend
#' \item \code{\link{greycurve}} (x) - Helper to determine the appropriate parameters for adjusting the grey level for msn functions
#' \item \code{\link{bruvo.msn}} (m) - Produces minimum spanning networks based off Bruvo's distance colored by population
#' \item \code{\link{poppr.msn}} (m | s | x) - Produces a minimum spanning network for any pairwise distance matrix related to the data
#' \item \code{\link{info_table}} (m) - Creates a heatmap representing missing data or observed ploidy
#' \item \code{\link{genotype_curve}} (m | x) - Creates a series of boxplots to demonstrate how many markers are needed to represent the diversity of your data. 
#' }
#' 
#' @section Datasets:
#' \itemize{
#' \item \code{\link{Aeut}} - (AFLP) Oomycete root rot pathogen
#' \emph{Aphanomyces euteiches} (Grünwald and Hoheisel, 2006)
#' \item \code{\link{monpop}} - (SSR) Peach brown rot pathogen \emph{Monilinia
#' fructicola} (Everhart and Scherm, 2015)
#' \item \code{\link{partial_clone}} - (SSR) partially-clonal data simulated via
#' simuPOP (Peng and Amos, 2008)
#' \item \code{\link{Pinf}} - (SSR) Potato late blight pathogen
#' \emph{Phytophthora infestans} (Goss et. al., 2014)
#' \item \code{\link{Pram}} - (SSR) Sudden Oak Death pathogen \emph{Phytophthora
#' ramorum} (Kamvar et. al., 2015; Goss et. al., 2009)
#' }
#' 
#' @author Zhian N. Kamvar, Jonah C. Brooks, Javier F. Tabima, Stacy
#'   Krueger-Hadfield, Erik Sotka, Niklaus J. Grünwald
#' 
#' Maintainer: Zhian N. Kamvar
#' 
#' @references --------- Papers announcing poppr ---------
#' 
#' Kamvar ZN, Tabima JF, Grünwald NJ. (2014) Poppr: an R package for genetic
#' analysis of populations with clonal, partially clonal, and/or sexual
#' reproduction. PeerJ 2:e281 \url{http://dx.doi.org/10.7717/peerj.281}
#' 
#' Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of 
#' genome-wide population genetic data with emphasis on clonality. Front. Genet.
#' 6:208. doi: 10.3389/fgene.2015.00208 
#' \url{http://dx.doi.org/10.3389/fgene.2015.00208}
#' 
#' --------- Papers referencing data sets ---------
#' 
#' Grunwald, NJ and Hoheisel, G.A. 2006. Hierarchical Analysis of Diversity, 
#' Selfing, and Genetic Differentiation in Populations of the Oomycete 
#' \emph{Aphanomyces euteiches}. Phytopathology 96:1134-1141 doi: 
#' \href{http://dx.doi.org/10.1094/PHYTO-96-1134}{10.1094/PHYTO-96-1134}
#' 
#' SE Everhart, H Scherm, (2015) Fine-scale genetic structure of \emph{Monilinia
#' fructicola} during brown rot epidemics within individual peach tree canopies.
#' Phytopathology 105:542-549 doi: 
#' \href{http://dx.doi.org/10.1094/PHYTO-03-14-0088-R}{10.1094/PHYTO-03-14-0088-R}
#' 
#' Bo Peng and Christopher Amos (2008) Forward-time simulations of nonrandom 
#' mating populations using simuPOP. \emph{bioinformatics}, 24 (11): 1408-1409.
#' 
#' Goss, Erica M., Javier F. Tabima, David EL Cooke, Silvia Restrepo, William E.
#' Fry, Gregory A. Forbes, Valerie J. Fieland, Martha Cardenas, and Niklaus J.
#' Grünwald. (2014) "The Irish potato famine pathogen \emph{Phytophthora 
#' infestans} originated in central Mexico rather than the Andes." Proceedings 
#' of the National Academy of Sciences 111:8791-8796. doi: 
#' \href{http://dx.doi.org/10.1073/pnas.1401884111}{10.1073/pnas.1401884111}
#' 
#' Kamvar, Z. N., Larsen, M. M., Kanaskie, A. M., Hansen, E. M., & Grünwald, N.
#' J. (2015). Spatial and temporal analysis of populations of the sudden oak
#' death pathogen in Oregon forests. Phytopathology 105:982-989. doi: 
#' \href{http://dx.doi.org/10.1094/PHYTO-12-14-0350-FI}{10.1094/PHYTO-12-14-0350-FI}
#' 
#' 
#' Goss, E. M., Larsen, M., Chastagner, G. A., Givens, D. R., and Grünwald, N. 
#' J. 2009. Population genetic analysis infers migration pathways of 
#' \emph{Phytophthora ramorum} in US nurseries. PLoS Pathog. 5:e1000583. doi: 
#' \href{http://dx.doi.org/10.1371/journal.ppat.1000583}{10.1371/journal.ppat.1000583}
#'   
#' 
#' @name poppr-package
#' @docType package
#==============================================================================#
NULL
