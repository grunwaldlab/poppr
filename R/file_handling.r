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
# This is simply a convenience function that serves as a wrapper for useful, but
# slightly obscure base R functions. 
#' Get a file name and path and store them in a list.
#'
#' getfile is a convenience function that serves as a wrapper for the functions
#' \code{\link{file.choose}, \link{file.path},} and \code{\link{list.files}}. 
#' If the user is working in a GUI environment, a window will pop up, allowing 
#' the user to choose a specified file regardless of path.
#'
#' @param multi this is an indicator to allow the user to store the names of
#' multiple files found in the directory. This is useful in conjunction with
#' \code{\link{poppr.all}}. 
#'
#' @param pattern a \code{\link{regex}} pattern for use while 
#' \code{multi == TRUE}. This will grab all files matching this pattern. 
#' 
#' @param combine \code{logical}. When this is set to \code{TRUE} (default), the
#' \code{$files} vector will have the path appended to them. When it is set to
#' \code{FALSE}, it will have the basename. 
#'
#' @return \item{path}{a character string of the absolute path to the
#' chosen file or files}
#' \item{files}{a character vector containing the chosen file
#' name or names.}
#' @author Zhian N. Kamvar
#'
#' @examples
#' \dontrun{
#'
#' x <- getfile()
#' poppr(x$files)
#'
#' y <- getfile(multi=TRUE, pattern="^.+?dat$") 
#' #useful for reading in multiple FSTAT formatted files.
#'
#' yfiles <- poppr.all(y$files)
#' 
#' # Write results to a file in that directory.
#' setwd(y$path)
#' write.csv(yfiles)
#' }  
#' @export
#==============================================================================#
getfile <- function(multi=FALSE, pattern=NULL, combine=TRUE){
  # the default option is to grab a single file in the directory. If multFile is 
  # set to TRUE, it will grab all the files in the directory corresponding to any
  # pattern that is set. If there is no pattern, all files will be grabbed.
  if (multi==TRUE){
    # this sets the path variable that the user can use to set the path
    # to the files with setwd(x$path), where x is the datastructure 
    # this function dumped into.
    pathandfile <- file.path(file.choose())
    path <- dirname(pathandfile)
    if (!is.null(pattern)){
      pat <- pattern
      x <- list.files(path, pattern=pat)
    }
    else {
      x <- list.files(path)
    }
  }
  else {
    # if the user chooses to analyze only one file, a pattern is not needed
    pathandfile <- file.path(file.choose())
    path <- dirname(pathandfile)
    x <- basename(pathandfile)
  }
  if(combine == TRUE){
    x <- paste(path, x, sep="/")
  }
  filepath <- list(files=x, path=path)
  return(filepath)
}

#==============================================================================#
#' Importing data from genalex formatted *.csv files.
#' 
#' read.genalex will read in a genalex-formatted file that has been exported in 
#' a comma separated format and will parse most types of genalex data. The 
#' output is a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} 
#' object.
#' 
#' @param genalex a *.csv file exported from genalex
#'   
#' @param ploidy indicate the ploidy of the dataset
#'   
#' @param geo indicates the presence of geographic data in the file. This data 
#'   will be included in a data frame labeled \code{xy} in the
#'   \code{\link{other}} slot.
#'   
#' @param region indicates the presence of regional data in the file.
#'   
#' @param genclone when \code{TRUE} (default), the output will be a
#'   \code{\linkS4class{genclone}} object. When \code{FALSE}, the output will be
#'   a \code{\linkS4class{genind}} object
#'   
#' @param sep A character specifying the column separator of the data. Defaults 
#'   to ",".
#'   
#' @param recode \strong{For polyploid data}: Do you want to recode your data to
#'   have varying ploidy? Default is \code{FALSE}, and the data will be returned
#'   with even ploidy where missing alleles are coded as "0". When \code{TRUE},
#'   the data is run through the function \code{\link{recode_polyploids}} before
#'   being returned. Note that this will prevent conversion to genpop objects in
#'   the future. See details.
#'   
#' @return A \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} 
#'   object.
#'   
#' @note This function cannot handle raw allele frequency data.
#'   
#'   In the case that there are duplicated names within the file, this function 
#'   will assume separate individuals and rename each one to a sequence of 
#'   integers from 1 to the number of individuals. A vector of the original 
#'   names will be saved in the \code{other} slot under \code{original_names}.
#'   
#'   
#' @details \subsection{if \code{genclone = FALSE}}{ The resulting genind object
#'   will have a data frame in the \code{other} slot called 
#'   \code{population_hierarchy}. This will contain a column for your population
#'   data and a column for your Regional data if you have set the flag.}
#'   
#'   \subsection{if \code{genclone = TRUE}}{ The resulting genclone object will 
#'   have a single strata defined in the strata slot. This will 
#'   be called "Pop" and will reflect the population factor defined in the 
#'   genalex input. If \code{region = TRUE}, a second column will be inserted 
#'   and labeled "Region". If you have more than two strata within 
#'   your data set, you should run the command \code{\link{splitStrata}} on 
#'   your data set to define the unique stratifications. }
#'   
#'   \subsection{FOR POLYPLOID (> 2n) DATA SETS}{ The genind object has 
#'   an all-or-none approach to missing data. If a sample has missing data at a 
#'   particular locus, then the entire locus is considered missing. This works 
#'   for diploids and haploids where allelic dosage is unambiguous. For 
#'   polyploids this poses a problem as much of the data set would be 
#'   transformed into missing data. With this function, I have created a 
#'   workaround.
#'   
#'   When importing polyploid data sets, missing data is scored as "0" and kept 
#'   within the genind object as an extra allele. This will break most analyses 
#'   relying on allele frequencies*. All of the functions in poppr will work 
#'   properly with these data sets as multilocus genotype analysis is agnostic 
#'   of ploidy and we have written both Bruvo's distance and the index of 
#'   association in such a way as to be able to handle polyploids presented in 
#'   this manner.
#'   
#'   * To restore functionality of analyses relying on allele frequencies, use
#'   the \code{\link{recode_polyploids}} function.}
#'   
#'   
#' @seealso \code{\link{clonecorrect}}, \code{\linkS4class{genclone}}, 
#'   \code{\linkS4class{genind}}, \code{\link{recode_polyploids}}
#'   
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' 
#' \dontrun{
#' Aeut <- read.genalex(system.file("files/rootrot.csv", package="poppr"))
#' 
#' genalex2 <- read.genalex("genalex2.csv", geo=TRUE)
#' # A genalex file with geographic coordinate data.
#' 
#' genalex3 <- read.genalex("genalex3.csv", region=TRUE) 
#' # A genalex file with regional information.
#' 
#' genalex4 <- read.genalex("genalex4.csv", region=TRUE, geo=TRUE) 
#' # A genalex file with both regional and geographic information.
#' }
#==============================================================================#

read.genalex <- function(genalex, ploidy = 2, geo = FALSE, region = FALSE, 
                         genclone = TRUE, sep = ",", recode = FALSE){
  # The first two lines from a genalex file contain all of the information about
  # the structure of the file (except for ploidy and geographic info)
  gencall  <- match.call()

  all.info <- strsplit(readLines(genalex, n = 2), sep)
  cskip    <- ifelse(inherits(genalex, "connection"), 0, 2)
  gena     <- utils::read.table(genalex, sep = sep, header = TRUE, skip = cskip, 
                                stringsAsFactors = FALSE, check.names = FALSE,
                                quote = "\"")
  num.info <- as.numeric(all.info[[1]])
  pop.info <- all.info[[2]][-c(1:3)]
  num.info <- num.info[!is.na(num.info)]
  pop.info <- pop.info[!pop.info %in% c("", NA)]
  nloci    <- num.info[1]
  ninds    <- num.info[2]
  npops    <- num.info[3]
  
  # Ensuring that all rows and columns have data
  data_rows <- apply(gena, 1, function(i) !all(is.na(i)))
  data_cols <- apply(gena, 2, function(i) !all(is.na(i)))
  gcols <- colnames(gena)
  dups  <- duplicated(gcols)
  if (any(dups) && any(gcols[dups] != "")){
    locations <- which(dups)
    locations <- locations[gcols[dups] != ""]
    colnames(gena) <- make.unique(gcols, sep = "_")
    the_loci  <- gcols[dups & gcols != ""]
    renamed   <- colnames(gena)[dups & gcols != ""]
    msg <- "\n I found duplicate column names in your data.\n They are being renamed:\n"
    renamed_data <- paste0("col ", locations, ": ", 
                           the_loci, " -> ", renamed, collapse = "\n ")
    msg <- paste(msg, renamed_data)
    warning(msg, immediate. = TRUE)
  }
  gena      <- gena[data_rows, data_cols]
  if (nrow(gena) != ninds) {
    theData <- if (inherits(genalex, "character")) genalex else deparse(substitute(genalex))
    msg <- paste0("\n The number of rows in your data do not match the number ",
                 "of individuals specified.",
                 "\n\t", ninds,      " individuals specified",
                 "\n\t", nrow(gena), " rows in data",
                 "\n Please inspect ", theData, " to ensure it's a properly ",
                 "formatted GenAlEx file.\n")
    stop(msg)
  }
  #----------------------------------------------------------------------------#
  # Checking for extra information such as Regions or XY coordinates
  #----------------------------------------------------------------------------#
  
  # Creating vectors that correspond to the different information fields. If the
  # regions are true, then the length of the pop.info should be equal to the
  # number of populations "npop"(npops) plus the number of regions which is the
  # npop+4th entry in the vector. Note that this strategy will only work if the
  # name of the first region does not match any of the populations.
  
  clm <- ncol(gena)

  if (region == TRUE && !is.na(num.info[npops + 4]) && length(pop.info) == npops + num.info[npops + 4]){
    # Info for the number of columns the loci can take on.
    loci.adj <- c(nloci, nloci*ploidy)
    
    # First question, do you have two or four extra columns? Two extra would 
    # indicate no geographic data. Four extra would indicate geographic data. 
    # Both of these indicate that, while a regional specification exists, a 
    # column indicating the regions was not specified, so it needs to be created
    if (((clm %in% (loci.adj + 4)) & (geo == TRUE)) | (clm %in% (loci.adj + 2))){
      
      pop.vec     <- gena[, 2]
      ind.vec     <- gena[, 1]
      xy          <- gena[, c((clm - 1), clm)]
      region.inds <- ((npops + 5):length(num.info)) # Indices for the regions
      reg.inds    <- num.info[region.inds] # Number of individuals per region
      reg.names   <- all.info[[2]][region.inds] # Names of the regions
      reg.vec     <- rep(reg.names, reg.inds) # Paste into a single vector
      names(reg.vec) <- ind.vec
      if (geo == TRUE){
        geoinds <- c((clm - 1), clm)
        xy      <- gena[, geoinds]
        gena    <- gena[, -geoinds, drop = FALSE]
      } else {
        xy <- NULL
      }
      gena <- gena[, c(-1, -2), drop = FALSE]
    } else {
      pop.vec      <- ifelse(any(gena[, 1] == pop.info[1]), 1, 2)
      reg.vec      <- ifelse(pop.vec == 2, 1, 2)
      orig.ind.vec <- NULL
      reg.vec      <- gena[, reg.vec] # Regional Vector 
      pop.vec      <- gena[, pop.vec] # Population Vector
      if (geo == TRUE){
        geoinds <- c((clm - 1), clm)
        xy      <- gena[, geoinds]
        gena    <- gena[, -geoinds, drop = FALSE]
      } else {
        xy <- NULL
      }
      ind.vec <- gena[, clm] # Individual Vector
      gena    <- gena[, -c(1, 2, clm), drop = FALSE] # removing the non-genotypic columns
    }
  } else if (geo == TRUE & length(pop.info) == npops){
    # There are no Regions specified, but there are geographic coordinates
    reg.vec <- NULL
    pop.vec <- gena[, 2]
    ind.vec <- gena[, 1]
    xy      <- gena[, c((clm - 1), clm)]
    gena    <- gena[, -c(1, 2, (clm - 1), clm), drop = FALSE]
  } else {
    # There are no Regions or geographic coordinates
    reg.vec <- NULL
    pop.vec <- gena[, 2]
    ind.vec <- gena[, 1]
    xy      <- NULL
    gena    <- gena[, -c(1, 2), drop = FALSE]
  }
  
  `%null%` <- function(a, b) if (is.null(a)) NULL else b
  
  xy_trail_na <- xy %null% (rev(cumsum(rev(rowSums(!is.na(xy))))) > 0)
  
  if (any(xy_trail_na)){
    xy         <- xy[xy_trail_na, , drop = FALSE]
    xy_nomatch <- TRUE
  } else {
    xy_nomatch <- is.null(xy)
  }
  # We need to set the names of the population vector to the original names.
  names(pop.vec) <- pop.vec %null% ind.vec
  names(reg.vec) <- reg.vec %null% ind.vec
  rownames(xy)   <- if (xy_nomatch) rownames(xy) else ind.vec
  
  #----------------------------------------------------------------------------#
  # The genotype matrix has been isolated at this point. Now this will
  # reconstruct the matrix in a way that adegenet likes it.
  #----------------------------------------------------------------------------#
  clm      <- ncol(gena)
  gena.mat <- as.matrix(gena)
  # Checking for greater than haploid data.
  if (nloci == clm/ploidy & ploidy > 1){
    # Missing data in genalex is coded as "0" for non-presence/absence data.
    # We will convert these data later, but this will take care of the 
    # data that's considered missing here
    gena[is.na(gena.mat)] <- "0"
    
    type  <- 'codom'
    # The remainder of the position when divided by ploidy is 1, which means
    # that this is the first column starting the locus
    loci  <- which(seq(clm) %% ploidy == 1)
    gena2 <- gena[, loci, drop = FALSE]

    # Collapsing all alleles over all loci into their respective loci.
    for (pos in loci){
      # Get the locus position in gena2. Example, in a triploid, the second
      # locus would be at gena[, 4:6], so the second element in pos would be 4
      # because position %% ploidy is 4 %% 3 == 1
      # We get our position out by:
      #  1. subtracting 1 to account for the modulo
      #  2. divide by ploidy
      #  3. add 1 to bring it up to R's 1-index
      pos_out <- ((pos - 1)/ploidy) + 1
      pos_in  <- pos:(pos + ploidy - 1) # allele positions in gena
      # paste columns of locus together
      donor   <- do.call("paste", c(gena[pos_in], list(sep = "/")))
      gena2[, pos_out] <- donor
    }

    # Treating missing data.
    gena2[gena2 == paste(rep("0", ploidy), collapse = "/")] <- NA_character_
    
    # Importing
    res.gid <- df2genind(gena2, sep="/", ind.names = ind.vec, pop = pop.vec,
                         ploidy = ploidy, type = type)
  } else if (nloci == clm & all(gena.mat %in% as.integer(-1:1))) {
    # Checking for AFLP data.
    # Missing data in genalex is coded as "-1" for presence/absence data.
    # this converts it to "NA" for adegenet.
    if(any(gena.mat == -1L)){
      gena[gena.mat == -1L] <- NA
    }
    type <- 'PA'
    res.gid <- df2genind(gena, ind.names = ind.vec, pop = pop.vec,
                         ploidy = ploidy, type = type, ncode = 1)
  } else if (nloci == clm & !all(gena.mat %in% as.integer(-1:1))) {
    # Checking for haploid microsatellite data or SNP data
    if (any(gena.mat == "0", na.rm = TRUE)){
      gena[gena.mat == "0"] <- NA
    }
    type    <- 'codom'
    res.gid <- df2genind(gena, ind.names = ind.vec, pop = pop.vec,
                         ploidy = 1, type = type)
  } else {
    weirdomsg <- paste("Something went wrong. Please ensure that your data is", 
                       "formatted correctly:\n\n",
                       "Is the ploidy flag set to the maximum ploidy of your data?\n",
                       ifelse(geo, 
                        "You set geo = TRUE. Do you have geographic coordinates in the right place?\n",
                        "If you have geographic coordinates: did you set the flag?\n"),
                       ifelse(region,  
                       "You set region = TRUE. Do you have regional data in the right place?\n\n",
                       "If you have regional data: did you set the flag?\n\n"),
                       "Otherwise, the problem may lie within the data structure itself.\n",
                       " - Inspect your data; if it looks fine then\n",
                       " - search the poppr forum for the error message; if you still can't find a solution:\n",
                       "   1. Make a minimum working example of your error (see http://stackoverflow.com/a/5963610)\n",
                       "   2. Use the 'reprex' package to reproduce the error (see https://github.com/jennybc/reprex#readme)\n",
                       "   3. Create a new issue on https://github.com/grunwaldlab/poppr/issues")
    stop(weirdomsg)
  }
  res.gid@call <- gencall
  # Checking for individual name duplications or removals -------------------
  
  same_names <- any(indNames(res.gid) %in% ind.vec)
  if (same_names){ # no duplications, only removals
    names(ind.vec) <- ind.vec
    ind.vec        <- ind.vec[indNames(res.gid)]
  } else {         # removals and/or duplciations
    ind.vec <- ind.vec[as.integer(indNames(res.gid))]
    other(res.gid)$original_names <- ind.vec
  }
  pop.vec <- pop.vec %null% pop.vec[ind.vec]
  reg.vec <- reg.vec %null% reg.vec[ind.vec]
  xy      <- if (xy_nomatch) xy else xy[ind.vec, , drop = FALSE]
  ind.vec <- indNames(res.gid)
  names(pop.vec) <- pop.vec %null% ind.vec
  names(reg.vec) <- reg.vec %null% ind.vec
  rownames(xy)   <- if (xy_nomatch) rownames(xy) else ind.vec
  
  # Keep the name if it's a URL ---------------------------------------------
  
  if (length(grep("://", genalex)) < 1 & !"connection" %in% class(genalex)){
    res.gid@call[2] <- basename(genalex)
  }
  if (region){
    strata(res.gid) <- data.frame(Pop = pop.vec, Region = reg.vec)
  } else {
    strata(res.gid) <- data.frame(Pop = pop.vec)
  }
  if (geo){
    if (!all(c("x", "y") %in% colnames(xy))){
      colnames(xy) <- c("x", "y")
    }
    other(res.gid)$xy <- xy
  }
  res.gid <- if (genclone) as.genclone(res.gid) else res.gid
  res.gid <- if (recode) recode_polyploids(res.gid, newploidy = TRUE) else res.gid

  return(res.gid)
}

#==============================================================================#
#' Exporting data from genind objects to genalex formatted *.csv files.
#' 
#' genind2genalex will export a genclone or genind object to a *.csv file
#' formatted for use in genalex.
#' 
#' @param gid a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}}
#'   object.
#'   
#' @param filename a string indicating the name and/or path of the file you wish
#'   to create.
#'   
#' @param quiet \code{logical} If \code{FALSE} a message will be printed to the 
#'   screen.
#'   
#' @param pop a character vector OR formula specifying the population factor. 
#'   This can be used to specify a specific subset of strata or custom 
#'   population factor for the output. Note that the \code{allstrata} command 
#'   has precedence over this unless the value of this is a new population
#'   factor.
#'   
#' @param allstrata if this is \code{TRUE}, the strata will be combined into a
#'   single population factor in the genalex file.
#' 
#' @param geo \code{logical} Default is \code{FALSE}. If it is set to 
#'   \code{TRUE}, the resulting file will have two columns for geographic data.
#'   
#' @param geodf \code{character} Since the \code{other} slot in the adegenet 
#'   object can contain many different items, you must specify the name of the 
#'   data frame in the \code{other} slot containing your geographic coordinates.
#'   It defaults to "xy".
#'   
#' @param sep a character specifying what character to use to separate columns. 
#'   Defaults to ",".
#'   
#' @param sequence when \code{TRUE}, sequence data will be converted to integers
#'   as per the GenAlEx specifications.
#'   
#' @note If you enter a file name that exists, that file will be overwritten. If
#'   your data set lacks a population structure, it will be coded in the new 
#'   file as a single population labeled "Pop". Likewise, if you don't have any
#'   labels for your individuals, they will be labeled as "ind1" through 
#'   "ind\emph{N}", with \emph{N} being the size of your population.
#'   
#' @seealso \code{\link{clonecorrect}}, \code{\linkS4class{genclone}} or
#'   \code{\linkS4class{genind}}
#'   
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' \dontrun{
#' data(nancycats)
#' genind2genalex(nancycats, "~/Documents/nancycats.csv", geo=TRUE)
#' }
#==============================================================================#
genind2genalex <- function(gid, filename = "genalex.csv", quiet = FALSE, pop = NULL, 
                           allstrata = TRUE, geo = FALSE, geodf = "xy", sep = ",",
                           sequence = FALSE){
  if (!is.genind(gid)) stop("A genind object is needed.")
  if (nchar(sep) != 1) stop("sep must be one byte/character (eg. \",\")")
  
  if (allstrata && !is.null(strata(gid))){
    if (is.null(pop) && length(pop) != nInd(gid)){
      allform <- paste0("~", paste(nameStrata(gid), collapse = "/"))
      allform <- stats::as.formula(allform)
      setPop(gid) <- allform
    } else {
      pop(gid) <- pop
    }
  } else if (!is.null(pop)){
    if ("formula" %in% class(pop)){
      setPop(gid) <- pop
    } else if (is.character(pop) && length(pop) == nInd(gid)){
      pop(gid) <- pop
    } else{
      warning("supplied population factor does not match the number of samples or strata.")
    }
  }
  
  if (is.null(pop(gid))){
    pop(gid) <- rep("Pop", nInd(gid))
  }
  popcall <- match.call()
  #topline is for the number of loci, individuals, and populations.
  topline <- c(nLoc(gid), nInd(gid), nPop(gid))
  popsizes <- table(pop(gid))
  # The sizes of the populations correspond to the second line, which is the_gid
  # names. 
  topline <- c(topline, popsizes)
  secondline <- c("", "", "", popNames(gid))
  ploid <- ploidy(gid)
  # Constructing the locus names. GenAlEx separates the alleles of the loci, so
  # There is one locus name for every p ploidy columns you have.
  if(all(ploid > 1) & gid@type == "codom"){
    # To intersperse spaces between the locus names, make a ploidy x loci
    # matrix, fill the first row with the loci names, fill the rest with
    # emptiness, and then convert it into a vector.
    locnames       <- matrix(character(nLoc(gid)*max(ploid)), nrow = max(ploid))
    locnames[1, ]  <- locNames(gid)
    locnames[-1, ] <- " "
    dim(locnames)  <- NULL
  } else {
    locnames <- locNames(gid)
  }
  thirdline <- c("Ind","Pop", locnames)
  
  # This makes sure that you don't get into a stacking error when stacking the
  # first three rows.
  if(length(thirdline) > length(topline)){
    lenfac <- length(thirdline) - length(topline)
    topline <- c(topline, rep("", lenfac))
    secondline <- c(secondline, rep("", lenfac))
  }
  else if(length(thirdline) < length(topline)){
    lenfac <- length(topline) - length(thirdline)
    thirdline <- c(thirdline, rep("", lenfac))
  }
  infolines <- rbind(topline, secondline, thirdline)

  if (sequence){
    alleles(gid) <- lapply(alleles(gid), function(i){
      i <- tolower(i)
      vapply(i, switch, integer(1), 
             a = 1L, c = 2L, g = 3L, t = 4L, 
             "-" = 5L, ":" = 5L, 0L)
    })
  }
  if(!quiet) cat("Extracting the table ... ")
  the_gid <- as.character(pop(gid))
  df      <- genind2df(gid, sep = "/", usepop = FALSE)
  if (any(ploid > 1)){
    df <- generate_bruvo_mat(df, maxploid = max(ploid), sep = "/", mat = TRUE)
  }
  df[is.na(df)] <- 0
  
  # making sure that the individual names are included.
  if(all(indNames(gid) == "") | is.null(indNames(gid))){
    indNames(gid) <- paste("ind", 1:nInd(gid), sep="")
  }
  df <- cbind(indNames(gid), the_gid, df)
  # setting the NA replacement. This doesn't work too well. 
  replacement <- ifelse(gid@type == "PA", "-1", "0")
  if(!quiet) cat("Writing the table to", filename, "... ")
  
  if(geo == TRUE & !is.null(gid$other[[geodf]])){
    replacemat <- matrix("", 3, 3)
    replacemat[3, 2:3] <- c("X", "Y")
    infolines <- cbind(infolines, replacemat)
    df2 <- data.frame(list("Space" = rep("", nInd(gid))))
    gdf <- as.matrix(gid@other[[geodf]])
    if (nrow(gdf) < nInd(gid)){
      gdf <- rbind(gdf, matrix("", nInd(gid) - nrow(gdf), 2))
    }
    df <- cbind(df, df2, gdf)
  } else if (geo == TRUE) {
    popcall <- popcall[2]
    warning(paste0("There is no data frame or matrix in ",
                  paste0(substitute(popcall)), "@other called ", geodf,
                  ".\nThe xy coordinates will not be represented in the",
                  " resulting file."))
  }
  
  df[df == "NA" | is.na(df)] <- replacement
  
  if (ncol(infolines) > ncol(df)){
    lendiff <- ncol(infolines) - ncol(df)
    padding <- matrix("", nrow = nInd(gid), ncol = lendiff)
    df      <- cbind(df, padding)
  }
  write.table(infolines, file = filename, quote = FALSE, row.names = FALSE, 
              col.names = FALSE, sep = sep)
  write.table(df, file = filename, quote = TRUE, na = replacement, append = TRUE, 
              row.names = FALSE, col.names = FALSE, sep = sep)
  if(!quiet) cat("Done.\n")
}
