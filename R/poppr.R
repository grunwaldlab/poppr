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
#' The \pkg{poppr} R package
#' 
#' @md
#' @description \pkg{Poppr} provides tools for population genetic analysis that
#' include genotypic diversity measures, genetic distances with bootstrap
#' support, native organization and handling of population hierarchies, and
#' clone correction.
#' 
#' To cite \pkg{poppr}, please use `citation("poppr")`. When referring to
#' \pkg{poppr} in your manuscript, please use lower case unless it occurs at the
#' beginning of a sentence.
#' 
#' @details This package relies on the \pkg{\link[adegenet]{adegenet}} package. 
#'   It is built around the [genind-class] and 
#'   [genlight-class] object. Genind objects store genetic
#'   information in a table of allele frequencies while genlight objects store
#'   SNP data efficiently by packing binary allele calls into single bits.
#'   \pkg{Poppr} has extended these object into new objects called 
#'   [genclone-class] and [snpclone-class],
#'   respectively. These objects are designed for analysis of clonal organisms
#'   as they add the **@@mlg** slot for keeping track of multilocus
#'   genotypes and multilocus lineages.
#'   
#'   \subsection{Documentation}{ Documentation is available for any function by 
#'   typing `?function_name` in the R console. Detailed topic explanations
#'   live in the package vignettes:
#'   
#'   \tabular{ll}{
#'     \strong{Vignette}              \tab \strong{command}\cr
#'     Data import and manipulation   \tab `vignette("poppr_manual", "poppr")`\cr
#'     Algorithms and Equations       \tab `vignette("algo", "poppr")`\cr
#'     Multilocus Genotype Analysis   \tab `vignette("mlg", "poppr")`  
#'   }
#'   
#'   Essential functions for importing and manipulating data are detailed within
#'   the *Data import and manipulation* vignette, details on algorithms
#'   used in \pkg{poppr} are within the *Algorithms and equations*
#'   vignette, and details for working with multilocus genotypes are in
#'   *Multilocus Genotype Analysis*. 
#'
#'   Examples of analyses are available in a primer written by Niklaus J.
#'   Grünwald, Zhian N. Kamvar, and Sydney E.  Everhart at
#'   <http://grunwaldlab.github.io/Population_Genetics_in_R>.}
#'   
#'   \subsection{Getting help}{ If you have a specific question or issue with 
#'   \pkg{poppr}, feel free to contribute to the google group at 
#'   <https://groups.google.com/group/poppr>. If you find a bug and 
#'   are a github user, you can submit bug reports at 
#'   <https://github.com/grunwaldlab/poppr/issues>. Otherwise, leave a
#'   message on the groups. Personal emails are highly discouraged as they do 
#'   not allow others to learn.}
#'   
#' @section Functions in \pkg{poppr}:
#' 
#' Below are descriptions and links to functions found in \pkg{poppr}. Be
#' aware that all functions in \pkg{\link[adegenet]{adegenet}} are also
#' available. The functions are documented as: 
#' 
#' - `function_name()` (data type) - Description
#' 
#' Where \sQuote{data type} refers to the type of data that can be used:
#' \tabular{ll}{
#'   \strong{m} \tab a genclone or genind object \cr
#'   \strong{s} \tab a snpclone or genlight object \cr
#'   \strong{x} \tab a different data type (e.g. a matrix from [mlg.table()])
#' }
#'
#' @section Data import/export:
#' 
#' - [getfile()] (x) - Provides a quick GUI to grab files for import
#' - [read.genalex()] (x) - Reads GenAlEx formatted csv files to a genind object
#' - [genind2genalex()] (m) - Converts genind objects to GenAlEx formatted csv files
#' - [genclone2genind()] (m) - Removes the @@mlg slot from genclone objects
#' - [as.genambig()] (m) - Converts genind data to \pkg{polysat}'s [genambig][polysat::genambig-class] data structure.
#' - [bootgen2genind()] (x) - see [aboot()] for details)
#'
#' @section Data Structures:
#' 
#' Data structures "genclone" (based off of adegenet's [genind][adegenet::genind-class]) and 
#' "snpclone" (based off of adegenet's [genlight][adegenet::genlight-class] for large SNP data sets).
#' Both of these data structures are defined by the presence of an extra
#' MLG slot representing multilocus genotype assignments, which can be a numeric
#' vector or a MLG class object.
#' 
#' - [genclone-class] - Handles microsatellite, presence/absence, and small SNP data sets
#' - [snpclone-class] - Designed to handle larger binary SNP data sets. 
#' - [MLG-class] - An internal class holding a data frame of multilocus genotype
#'   assignments that acts like a vector, allowing the user to easily switch
#'   between different MLG definitions.
#' - [bootgen-class] - An internal class used explicitly for [aboot()] that
#'   inherits the [gen-class][adegenet::gen-class] virtual object. It is
#'   designed to allow for sampling loci with replacement.
#' - [bruvomat-class] - An internal class designed to handle bootstrapping for
#'   Bruvo's distance where blocks of integer loci can be shuffled.
#' 
#' @section Data manipulation:
#' 
#' - [as.genclone()] (m) - Converts genind objects to genclone objects
#' - [missingno()] (m) - Handles missing data
#' - [clonecorrect()] (m | s) - Clone-censors at a specified population hierarchy
#' - [informloci()] (m) - Detects and removes phylogenetically uninformative loci
#' - [popsub()] (m | s) - Subsets genind objects by population
#' - [shufflepop()] (m) - Shuffles genotypes at each locus using four different shuffling algorithms
#' - [recode_polyploids()] (m | x) - Recodes polyploid data sets with missing alleles imported as "0"
#' - [make_haplotypes()] (m | s) - Splits data into pseudo-haplotypes. This is mainly used in AMOVA.
#' - [test_replen()] (m) - Tests for inconsistent repeat lengths in microsatellite data. For use in [bruvo.dist()] functions.
#' - [fix_replen()] (m) - Fixes inconsistent repeat lengths. For use in [bruvo.dist()] functions.
#' 
#' @section Genetic distances:
#' 
#' - [bruvo.dist()] (m) - Bruvo's distance (see also: [fix_replen()])
#' - [diss.dist()] (m) - Absolute genetic distance (see [prevosti.dist()])
#' - [nei.dist()] (m | x) - Nei's 1978 genetic distance
#' - [rogers.dist()] (m | x) - Rogers' euclidean distance
#' - [reynolds.dist()] (m | x) - Reynolds' coancestry distance
#' - [edwards.dist()] (m | x) - Edwards' angular distance
#' - [prevosti.dist()] (m | x) - Prevosti's absolute genetic distance
#' - [bitwise.dist()] (s) - Calculates fast pairwise distances for genlight objects. 
#'  
#' @section Bootstrapping:
#' 
#' - [aboot()] (m | s | x) - Creates a bootstrapped dendrogram for any distance measure
#' - [bruvo.boot()] (m) - Produces dendrograms with bootstrap support based on Bruvo's distance
#' - [diversity_boot()] (x) - Generates boostrap distributions of diversity statistics for multilocus genotypes
#' - [diversity_ci()] (m | s | x) - Generates confidence intervals for multilocus genotype diversity.
#' - [resample.ia()] (m) - Calculates the index of association over subsets of data.
#'  
#' @section Multilocus Genotypes:
#' 
#' - [mlg()] (m | s) - Calculates the number of multilocus genotypes
#' - [mll()] (m | s) - Displays the current multilocus lineages (genotypes) defined.
#' - [nmll()] (m | s) - Same as [mlg()].
#' - [mlg.crosspop()] (m | s) - Finds all multilocus genotypes that cross populations
#' - [mlg.table()] (m | s) - Returns a table of populations by multilocus genotypes
#' - [mlg.vector()] (m | s) - Returns a vector of a numeric multilocus genotype assignment for each individual
#' - [mlg.id()] (m | s) - Finds all individuals associated with a single multilocus genotype
#' - [mlg.filter()] (m | s) - Collapses MLGs by genetic distance
#' - [filter_stats()] (m | s) - Calculates mlg.filter for all algorithms and plots
#' - [cutoff_predictor()] (x) - Predicts cutoff threshold from mlg.filter. 
#' - [mll.custom()] (m | s) - Allows for the custom definition of multilocus lineages
#' - [mll.levels()] (m | s) - Allows the user to change levels of custom MLLs. 
#' - [mll.reset()] (m | s) - Reset multilocus lineages. 
#' - [diversity_stats()] (x) - Creates a table of diversity indices for multilocus genotypes. 
#'    
#' 
#' @section Index of Association Analysis:
#'
#' Analysis of multilocus linkage disequilibrium.
#' 
#' - [ia()] (m) - Calculates the index of association
#' - [pair.ia()] (m) - Calculates the index of association for all loci pairs.
#' - [win.ia()] (s) - Index of association windows for genlight objects.
#' - [samp.ia()] (s) - Index of association on random subsets of loci for genlight objects.
#' 
#' @section Population Genetic Analysis:
#'
#' - [poppr.amova()] (m | s) - Analysis of Molecular Variance (as implemented in ade4)
#' - [poppr()] (m | x) - Returns a diversity table by population
#' - [poppr.all()] (m | x) - Returns a diversity table by population for all compatible files specified
#' - [private_alleles()] (m) - Tabulates the occurrences of alleles that only occur in one population.
#' - [locus_table()] (m) - Creates a table of summary statistics per locus.
#' - [rrmlg()] (m | x) - Round-robin multilocus genotype estimates.
#' - [rraf()] (m) - Round-robin allele frequency estimates.
#' - [pgen()] (m) - Probability of genotypes.
#' - [psex()] (m) - Probability of observing a genotype more than once.
#' - [rare_allele_correction][poppr::rare_allele_correction] (m) - rules for correcting rare alleles for round-robin estimates.
#' - [incomp()] (m) - Check data for incomparable samples.
#' 
#'
#' @section Visualization:
#' 
#' - [imsn()] (m | s) - Interactive construction and visualization of minimum spanning networks
#' - [plot_poppr_msn()] (m | s | x) - Plots minimum spanning networks produced in poppr with scale bar and legend
#' - [greycurve()] (x) - Helper to determine the appropriate parameters for adjusting the grey level for msn functions
#' - [bruvo.msn()] (m) - Produces minimum spanning networks based off Bruvo's distance colored by population
#' - [poppr.msn()] (m | s | x) - Produces a minimum spanning network for any pairwise distance matrix related to the data
#' - [info_table()] (m) - Creates a heatmap representing missing data or observed ploidy
#' - [genotype_curve()] (m | x) - Creates a series of boxplots to demonstrate how many markers are needed to represent the diversity of your data. 
#' 
#' @section Datasets:
#' 
#' - [Aeut()] - (AFLP) Oomycete root rot pathogen *Aphanomyces euteiches* (Grünwald and Hoheisel, 2006)
#' - [monpop()] - (SSR) Peach brown rot pathogen *Monilinia fructicola* (Everhart and Scherm, 2015)
#' - [partial_clone()] - (SSR) partially-clonal data simulated via simuPOP (Peng and Amos, 2008)
#' - [Pinf()] - (SSR) Potato late blight pathogen *Phytophthora infestans* (Goss et. al., 2014)
#' - [Pram()] - (SSR) Sudden Oak Death pathogen *Phytophthora ramorum* (Kamvar et. al., 2015; Goss et. al., 2009)
#' 
#' @author Zhian N. Kamvar, Jonah C. Brooks, Sydney E. Everhart, Javier F. 
#'   Tabima, Stacy Krueger-Hadfield, Erik Sotka, Niklaus J. Grünwald
#' 
#' Maintainer: Zhian N. Kamvar
#' 
#' @references --------- Papers announcing poppr ---------
#' 
#' Kamvar ZN, Tabima JF, Grünwald NJ. (2014) Poppr: an R package for genetic
#' analysis of populations with clonal, partially clonal, and/or sexual
#' reproduction. PeerJ 2:e281 <https://doi.org/10.7717/peerj.281>
#' 
#' Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of 
#' genome-wide population genetic data with emphasis on clonality. Front. Genet.
#' 6:208. doi: 10.3389/fgene.2015.00208 
#' <https://doi.org/10.3389/fgene.2015.00208>
#' 
#' --------- Papers referencing data sets ---------
#' 
#' Grünwald, NJ and Hoheisel, G.A. 2006. Hierarchical Analysis of Diversity, 
#' Selfing, and Genetic Differentiation in Populations of the Oomycete 
#' *Aphanomyces euteiches*. Phytopathology 96:1134-1141 doi: 
#' \href{https://doi.org/10.1094/PHYTO-96-1134}{10.1094/PHYTO-96-1134}
#' 
#' SE Everhart, H Scherm, (2015) Fine-scale genetic structure of *Monilinia
#' fructicola* during brown rot epidemics within individual peach tree canopies.
#' Phytopathology 105:542-549 doi: 
#' \href{https://doi.org/10.1094/PHYTO-03-14-0088-R}{10.1094/PHYTO-03-14-0088-R}
#' 
#' Bo Peng and Christopher Amos (2008) Forward-time simulations of nonrandom 
#' mating populations using simuPOP. *bioinformatics*, 24 (11): 1408-1409.
#' 
#' Goss, Erica M., Javier F. Tabima, David EL Cooke, Silvia Restrepo, William E.
#' Fry, Gregory A. Forbes, Valerie J. Fieland, Martha Cardenas, and Niklaus J.
#' Grünwald. (2014) "The Irish potato famine pathogen *Phytophthora 
#' infestans* originated in central Mexico rather than the Andes." Proceedings 
#' of the National Academy of Sciences 111:8791-8796. doi: 
#' \href{https://doi.org/10.1073/pnas.1401884111}{10.1073/pnas.1401884111}
#' 
#' Kamvar, Z. N., Larsen, M. M., Kanaskie, A. M., Hansen, E. M., & Grünwald, N.
#' J. (2015). Spatial and temporal analysis of populations of the sudden oak
#' death pathogen in Oregon forests. Phytopathology 105:982-989. doi: 
#' \href{https://doi.org/10.1094/PHYTO-12-14-0350-FI}{10.1094/PHYTO-12-14-0350-FI}
#' 
#' 
#' Goss, E. M., Larsen, M., Chastagner, G. A., Givens, D. R., and Grünwald, N. 
#' J. 2009. Population genetic analysis infers migration pathways of 
#' *Phytophthora ramorum* in US nurseries. PLoS Pathog. 5:e1000583. doi: 
#' \href{https://doi.org/10.1371/journal.ppat.1000583}{10.1371/journal.ppat.1000583}
#'   
#' 
#' @name poppr-package
#' @docType package
#==============================================================================#
NULL
#==============================================================================#
#' Oomycete root rot pathogen *Aphanomyces euteiches* AFLP data
#' 
#' @name Aeut
#' @docType data
#' @usage data(Aeut)
#' @description The Aeut dataset consists of 187 isolates of the Oomycete root 
#'   rot pathogen, *Aphanomyces euteiches* collected from two different 
#'   fields in NW Oregon and W Washington, USA.
#' @format a [genind()] object with two populations containing a data
#'   frame in the `other` slot called `population_hierarchy`. This
#'   data frame gives indices of the populations and subpopulations for the data
#'   set.
#' @references Grunwald, NJ and Hoheisel, G.A. 2006. Hierarchical Analysis of
#'   Diversity, Selfing, and Genetic Differentiation in Populations of the 
#'   Oomycete *Aphanomyces euteiches*. Phytopathology 96:1134-1141
#'   doi: \href{https://doi.org/10.1094/PHYTO-96-1134}{10.1094/PHYTO-96-1134}
#==============================================================================#
NULL
#==============================================================================#
#' Simulated data illustrating a Minimum Spanning Network based on Bruvo's
#' Distance
#' 
#' @name partial_clone
#' @rdname partial_clone
#' @docType data
#' @usage data(partial_clone)
#' @description These data were simulated using SimuPOP version 1.0.8 with 
#'   99.9\% clonal reproduction over 10,000 generations. Populations were
#'   assigned post-hoc and are simply present for the purposes of demonstrating
#'   a minimum spanning network with Bruvo's distance.
#' @format a [genind()] object with 50 individuals, 10 loci, and four 
#'   populations.
#' @references Bo Peng and Christopher Amos (2008) Forward-time simulations of 
#'   nonrandom mating populations using simuPOP. *bioinformatics*, 24 (11):
#'   1408-1409.
#==============================================================================#
NULL
#==============================================================================#
#' @name old_partial_clone
#' @rdname partial_clone
#==============================================================================#
NULL
#==============================================================================#
#' Phytophthora infestans data from Mexico and South America.
#' 
#' @name Pinf
#' @rdname Pinf
#' @docType data
#' @usage data(Pinf)
#' @description The Pinf data set contains 86 isolates genotyped over 11 
#'   microsatellite loci collected from Mexico, Peru, Columbia, and Ecuador. 
#'   This is a subset of the data used for the reference below.
#' @format a [genclone-class] object with 2 hierarchical levels 
#'   called "Continent" and "Country" that contain 2 and 4 populations, 
#'   respectively.
#' @references Goss, Erica M., Javier F. Tabima, David EL Cooke, Silvia 
#'   Restrepo, William E. Fry, Gregory A. Forbes, Valerie J. Fieland, Martha 
#'   Cardenas, and Niklaus J. Grünwald. "The Irish potato famine pathogen 
#'   *Phytophthora infestans* originated in central Mexico rather than the 
#'   Andes." Proceedings of the National Academy of Sciences 111:8791-8796. doi:
#'   \href{https://doi.org/10.1073/pnas.1401884111}{10.1073/pnas.1401884111}
#==============================================================================#
NULL
#==============================================================================#
#' @name old_Pinf
#' @rdname Pinf
#==============================================================================#
NULL
#==============================================================================#
#' Phytophthora ramorum data from OR Forests and Nurseries (OR and CA)
#' 
#' @name Pram
#' @rdname Pram
#' @docType data
#' @usage data(Pram)
#' @description This is the data set from 
#'   <https://doi.org/10.5281/zenodo.13007>. It has been converted to the 
#'   genclone object as of poppr version 2.0. It contains 729 samples of the 
#'   Sudden Oak Death pathogen *Phytophthora ramorum* genotyped over five 
#'   microsatellite loci (Kamvar et. al., 2015). 513 samples were collected from
#'   forests in Curry County, OR from 2001 to mid-2014 (labeled by watershed
#'   region). The other 216 samples represents genotypes collected from
#'   Nurseries in OR and CA from Goss et. al. (2009).
#'   
#' @format a [genclone-class] object with 3 hierarchical levels 
#'   called "SOURCE", "YEAR", and, "STATE". The \strong{other} slot contains a 
#'   named vector of repeat lengths called \strong{"REPLEN"}, a matrix of xy 
#'   coordinates for the forest samples called \strong{"xy"}, and a palette to 
#'   color the ~SOURCE/STATE stratification called \strong{"comparePal"}.
#'   
#' @references Kamvar, Z. N., Larsen, M. M., Kanaskie, A. M., Hansen, E. M., & 
#'   Grünwald, N. J. (2015). Spatial and temporal analysis of populations of the
#'   sudden oak death pathogen in Oregon forests. Phytopathology 105:982-989. 
#'   doi:
#'   \href{https://doi.org/10.1094/PHYTO-12-14-0350-FI}{10.1094/PHYTO-12-14-0350-FI}
#'   
#'   
#'   Zhian N. Kamvar, Meg M. Larsen, Alan M. Kanaskie, Everett M. Hansen, & 
#'   Niklaus J. Grünwald. 2014. Sudden_Oak_Death_in_Oregon_Forests: Spatial and 
#'   temporal population dynamics of the sudden oak death epidemic in Oregon 
#'   Forests. ZENODO, doi:
#'   \href{https://doi.org/10.5281/zenodo.13007}{10.5281/zenodo.13007}
#'   
#'   Goss, E. M., Larsen, M., Chastagner, G. A., Givens, D. R., and Grünwald, N.
#'   J. 2009. Population genetic analysis infers migration pathways of 
#'   *Phytophthora ramorum* in US nurseries. PLoS Pathog. 5:e1000583. doi:
#'   \href{https://doi.org/10.1371/journal.ppat.1000583}{10.1371/journal.ppat.1000583}
#'   
#' @examples
#' data(Pram)
#' 
#' # Repeat lengths (previously processed via fix_replen)
#' other(Pram)$REPLEN
#' 
#' # Color palette for source by state. Useful for minimum spanning networks
#' other(Pram)$comparePal
#==============================================================================#
NULL
#==============================================================================#
#' Peach brown rot pathogen *Monilinia fructicola*
#' 
#' @name monpop
#' @docType data
#' @usage data(monpop)
#' @description This is microsatellite data for a population of the haploid 
#'   plant pathogen *Monilinia fructicola* that causes disease within peach
#'   tree canopies (Everhart & Scherm, 2014). Entire populations within trees
#'   were sampled across 3 years (2009, 2010, and 2011) in a total of four
#'   trees, where one tree was sampled in all three years, for a total of 6
#'   within-tree populations. Within each year, samples in the spring were taken
#'   from affected blossoms (termed "BB" for blossom blight) and in late summer
#'   from affected fruits (termed "FR" for fruit rot). There are a total of 694 
#'   isolates with 65 to 173 isolates within each canopy population that were 
#'   characterized using a set of 13 microsatellite markers.
#' @format a [genclone-class] object with 3 hierarchical levels 
#'   coded into one population factor. These are named "Tree", "Year", and 
#'   "Symptom"
#' @references SE Everhart, H Scherm, (2015) Fine-scale genetic structure of 
#'   *Monilinia fructicola* during brown rot epidemics within individual
#'   peach tree canopies. Phytopathology 105:542-549 doi:
#'   \href{https://doi.org/10.1094/PHYTO-03-14-0088-R}{10.1094/PHYTO-03-14-0088-R}
#'   
#' @examples
#' data(monpop)
#' splitStrata(monpop) <- ~Tree/Year/Symptom
#' setPop(monpop) <- ~Symptom/Year
#' monpop
#==============================================================================#
NULL
#' @importFrom graphics abline axis frame hist layout legend lines par plot.new
#'   points polygon rasterImage text title
#' @importFrom grDevices as.raster gray grey topo.colors
#' @importFrom stats as.dist as.formula df dist median quantile rmultinom sd
#'   setNames terms update var dbinom printCoefmat
#' @importFrom utils combn head read.table setTxtProgressBar tail txtProgressBar
#'   write.table capture.output
#' @importFrom dplyr progress_estimated
#' @importFrom polysat Genotypes Ploidies PopInfo Missing PopNames
#' @useDynLib poppr, .registration = TRUE
NULL
