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
#' \code{multFile==TRUE}.
#' 
#' @param combine \code{logical}. When this is set to \code{TRUE}, the
#' \code{\$files} vector will have the path appended to them. When it is set to
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
#' 
#'
#' y <- getfile(multFile=TRUE, pattern="^.+?dat$") 
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
# A way of dealing with the different types of data that adegenet can take in.
# This will detect whether or not one should drop all the non-informative loci
# from each separated population or not.
#==============================================================================#
.pop.divide <- function(x, drop=TRUE) {
  divcall <- match.call()
  if(!is.genind(x)){
    stop(c(as.character(divcall[2])," is not a valid genind object"))
  }
  if (is.null(pop(x))){
    pops <- NULL
  }
  #else if (!is.null(x@other$subpops)){
  #  pops <- lapply(seppop(x), function(y) seppop(y, pop=y@other$subpops))
  #}
  else if (x@type !="PA") {
    pops <- seppop(x, drop=drop)
  }
  else {
    pops <- seppop(x)
  }
	return(pops)
}
#==============================================================================#
#' Importing data from genalex formatted *.csv files.
#'
#' read.genalex will read in a genalex-formatted file that has been exported in
#' a comma separated format and will parse most types of genalex data. The
#' output is a \code{\link{genind}} object.
#'
#' @param genalex a *.csv file exported from genalex
#' 
#' @param ploidy indicate the ploidy of the dataset
#'
#' @param geo indicates the presence of geographic data in the file.
#' 
#' @param region indicates the presence of regional data in the file. 
#' 
#' @return A \code{\link{genind}} object. 
#'
#' @note This function cannot handle raw allele frequency data. 
#' 
#' The resulting genind object will have a data frame in the \code{other} slot called
#' population_hierarchy. This will contain a column for your population data and
#' a column for your Regional data if you have set the flag. 
#' 
#' If there is geographic data, it will be included in a data frame called xy in
#' the \code{other} slot. 
#' 
#' In the case that there are duplicated names within the file, this function 
#' will assume separate individuals and rename each one to a sequence of 
#' integers from 1 to the number of individuals. A vector of the original names
#' will be saved in the \code{other} slot under \code{original_names}.
#'
#' @seealso \code{\link{clonecorrect}}, \code{\link{genind}}
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

read.genalex <- function(genalex, ploidy=2, geo=FALSE, region=FALSE){
  # The first two lines from a genalex file contain all of the information about
  # the structure of the file (except for ploidy and geographic info)
  gencall <- match.call()
  all.info <- strsplit(readLines(genalex, n=2), ",")
  if (any(all.info[[1]]=="")){
    num.info <- as.numeric(all.info[[1]][-which(all.info[[1]]=="")])
  }
  else {
    num.info <- as.numeric(all.info[[1]])
  }
  if (any(all.info[[2]]=="")){
    pop.info <- all.info[[2]][c(-1,-2,-3,-which(all.info[[2]]==""))]  
  }
  else {
    pop.info <- all.info[[2]][c(-1,-2,-3)]
  }
  
  # glob.info is for the number of loci, individuals, and populations.
  glob.info <- num.info[1:3]
  if(any(is.na(glob.info))){
    stop("Something is wrong with your csv file. Perhaps it is not comma delimited?\n")
  }
  #   cat("Global Information:",glob.info,"\n")
  #   cat("Populations:",pop.info,"\n")
  #   cat("num.info:",num.info,"\n")
  
  gena <- read.csv(genalex, header=TRUE, skip=2)
  
  # Removing all null columns 
  if(!is.na(which(is.na(gena[1, ]))[1])){
    gena <- gena[, -which(is.na(gena[1, ]))]
  }
  
  #----------------------------------------------------------------------------#
  # Checking for extra information such as Regions or XY coordinates
  #----------------------------------------------------------------------------#
  
  # Creating vectors that correspond to the different information fields.
  # If the regions are true, then the length of the pop.info should be equal to
  # the number of populations "npop"(glob.info[3]) plus the number of regions
  # which is the npop+4th entry in the vector. 
  # Note that this strategy will only work if the name of the first region does
  # not match any of the populations. 
  
  if (region==TRUE & length(pop.info) == glob.info[3] + num.info[glob.info[3]+4]){
    # Info for the number of columns the loci can take on.
    loci.adj <- c(glob.info[1], glob.info[1]*ploidy)
    
    # First question, do you have two or four extra columns? Two extra would
    # indicate no geographic data. Four extra would indicate geographic data.
    # Both of these indicate that, while a regional specification exists, a 
    # column indicating the regions was not specified, so it needs to be created
    if(((ncol(gena) %in% (loci.adj + 4)) & (geo == TRUE)) | (ncol(gena) %in% (loci.adj + 2))){
      pop.vec <- gena[, 2]
      ind.vec <- gena[, 1]
      xy <- gena[, c((ncol(gena)-1), ncol(gena))]
      
      # Get the indices for the regions
      region.inds <- ((glob.info[3]+5):length(num.info))
      # Get the number of individuals per region
      reg.inds <- num.info[region.inds]
      # Get the names of the regions
      reg.names <- all.info[[2]][region.inds]
      # Paste them all into a single vector.
      reg.vec <- rep(reg.names, reg.inds)
      if(geo == TRUE){
        geoinds <- c((ncol(gena)-1), ncol(gena))
        xy <- gena[, geoinds]
        gena <- gena[, -geoinds]
      }
      else{
        xy <- NULL
      }
      gena <- gena[, c(-1,-2)]
    }
    #
    # The Regions are specified in one of the first two columns.
    #
    else{
      pop.vec <- ifelse(any(gena[, 1] == pop.info[1]), 1, 2)
      reg.vec <- ifelse(pop.vec == 2, 1, 2)
      orig.ind.vec <- NULL
      # Regional Vector    
      reg.vec <- gena[, reg.vec]
      # Population Vector
      pop.vec <- gena[, pop.vec]    
      if(geo == TRUE){
        geoinds <- c((ncol(gena)-1), ncol(gena))
        xy <- gena[, geoinds]
        gena <- gena[, -geoinds]
      }
      else{
        xy <- NULL
      }
      # Individual Vector
      ind.vec <- gena[, ncol(gena)]
      # removing the non-genotypic columns from the data frame
      gena <- gena[, c(-1,-2,-ncol(gena))]
    }
  }
  
  #
  # There are no Regions specified, but there are geographic coordinates
  #
  else if (geo == TRUE & length(pop.info) == glob.info[3]){
    reg.vec <- NULL
    pop.vec <- gena[, 2]
    ind.vec <- gena[, 1]
    xy <- gena[, c((ncol(gena)-1), ncol(gena))]
    gena <- gena[, c(-1,-2,-(ncol(gena)-1),-ncol(gena))]
  }
  #
  # There are no Regions or geographic coordinates
  #
  else{
    reg.vec <- NULL
    pop.vec <- gena[, 2]
    ind.vec <- gena[, 1]
    xy <- NULL
    gena <- gena[, c(-1,-2)]
  }
  
  #----------------------------------------------------------------------------#
  # The genotype matrix has been isolated at this point. Now this will
  # reconstruct the matrix in a way that adegenet likes it.
  #----------------------------------------------------------------------------#
  
  clm <- ncol(gena)
  gena.mat <- as.matrix(gena)
  # Checking for diploid data.
  if (glob.info[1] == clm/2){
    # Missing data in genalex is coded as "0" for non-presence/absence data.
    # this converts it to "NA" for adegenet.
    if(any(gena.mat =="0")){
      #gena.mat <- as.matrix(gena)
      gena[gena.mat == "0"] <- NA
      #gena <- as.data.frame(gena.mat)
    }
    type <- 'codom'
    loci <- which((1:clm)%%2==1)
    gena2 <- gena[, loci]
    lapply(loci, function(x) gena2[, ((x-1)/2)+1] <<-
             paste(gena[, x],"/",gena[, x+1], sep=""))
    #res <- list(Gena=gena2, Glob.info=glob.info, Ploid=ploidy)
    res.gid <- df2genind(gena2, sep="/", ind.names=ind.vec, pop=pop.vec,
                         ploidy=ploidy, type=type)
  }
  # Checking for AFLP data.
  else if (glob.info[1] == clm & all(gena.mat %in% as.integer(-1:1))) {
    # Missing data in genalex is coded as "-1" for presence/absence data.
    # this converts it to "NA" for adegenet.
    if(any(gena.mat == -1L)){
      #gena.mat <- as.matrix(gena)
      gena[gena.mat == -1L] <- NA
      #gena <- as.data.frame(gena.mat)
    }
    type <- 'PA'
    #res <- list(Gena=gena, Glob.info=glob.info, Ploid=ploidy)
    res.gid <- df2genind(gena, ind.names=ind.vec, pop=pop.vec,
                         ploidy=ploidy, type=type)
  }
  # Checking for haploid microsattellite data or SNP data
  else if (glob.info[1] == clm & !all(gena.mat %in% as.integer(-1:1))) {
    if(any(gena.mat == "0")){
      #gena.mat <- as.matrix(gena)
      gena[gena.mat == "0"] <- NA
      #gena <- as.data.frame(gena.mat)
    }
    type <- 'codom'
    #res <- list(Gena=gena, Glob.info=glob.info, Ploid=1)
    res.gid <- df2genind(gena, ind.names=ind.vec, pop=pop.vec,
                         ploidy= 1, type=type)
  }
  else {
    stop("Something went wrong. Check your geo and region flags to make sure they are set correctly. Otherwise, the problem may lie within the data structure itself.")
  }
  if (any(duplicated(ind.vec))){
    # preserving the names
    orig.ind.vec <- ind.vec
    # ensuring that all names are unique
    ind.vec <- 1:length(ind.vec)
    res.gid@other[["original_names"]] <- orig.ind.vec
  }
  
  res.gid@other[["population_hierarchy"]] <- as.data.frame(list(Pop=pop.vec))
  res.gid@call <- gencall
  res.gid@call[2] <- basename(genalex)
  if(region==TRUE){
    res.gid@other[["population_hierarchy"]]$Region <- reg.vec
    #return(res.gid)
  }
  if(geo==TRUE){
    res.gid@other[["xy"]] <- xy
    #return(res.gid)
  }
  return(res.gid)
}

#==============================================================================#
#' Exporting data from genind objects to genalex formatted *.csv files.
#'
#' genind2genalex will export a genind object to a *.csv file formatted for use
#' in genalex. 
#'
#' @param pop a \code{\link{genind}} object.
#' 
#' @param filename a string indicating the name and/or path of the file you wish
#' to create.
#'
#' @param quiet \code{logical} If \code{FALSE} a message will be printed to the
#' screen.
#' 
#' @param geo \code{logical} Default is \code{FALSE}. If it is set to 
#' \code{TRUE}, the resulting file will have two columns for geographic data.
#' 
#' @param geodf \code{character} Since the \code{other} slot in the adegenet object
#' can contain many different items, you must specify the name of the data frame
#' in the \code{other} slot containing your geographic coordinates. It defaults to
#' "xy". 
#' 
#' @note If you enter a file name that exists, that file will be overwritten.
#' If your data set lacks a population structure, it will be coded in the new
#' file as a single population lableled "Pop". Likewise, if you don't have any
#' labels for your individuals, they will be labeled as "ind1" through 
#' "ind\emph{N}", with \emph{N} being the size of your population.
#'
#' @seealso \code{\link{clonecorrect}}, \code{\link{genind}}
#'
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' \dontrun{
#' data(nancycats)
#' genind2genalex(nancycats, "~/Documents/nancycats.csv", geo=TRUE)
#' }
#==============================================================================#
genind2genalex <- function(pop, filename="genalex.csv", quiet=FALSE, geo=FALSE, geodf="xy"){
  if(!is.genind(pop)) stop("A genind object is needed.")
  if(is.null(pop@pop)){
    pop(pop) <- rep("Pop", nInd(pop))
  }
  popcall <- match.call()
  #topline is for the number of loci, individuals, and populations.
  topline <- c(nLoc(pop), nInd(pop), length(pop@pop.names))
  popsizes <- table(pop@pop)
  # The sizes of the populations correspond to the second line, which is the pop
  # names. 
  topline <- c(topline, popsizes)
  secondline <- c("", "", "", pop@pop.names)
  ploid <- ploidy(pop)
  # Constructing the locus names. GenAlEx separates the alleles of the loci, so
  # There is one locus name for every p ploidy columns you have.
  if(ploid > 1 & pop@type == "codom"){
    locnames <- unlist(strsplit(paste(pop@loc.names, 
                                      paste(rep(" ", ploidy(pop)-1), 
                                            collapse="/"), sep="/"),"/"))
  }
  else{
    locnames <- pop@loc.names
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
  
  # converting to a data frame
  if(any(!pop@tab %in% c(0, ((1:ploid)/ploid), 1, NA))){
    pop@tab[!pop@tab %in% c(0, ((1:ploid)/ploid), 1, NA)] <- NA
  }
  if(!quiet) cat("Extracting the table ... ")
  df <- genind2df(pop, oneColPerAll=TRUE)
  
  # making sure that the individual names are included.
  if(all(pop@ind.names == "") | is.null(pop@ind.names)){
    pop@ind.names <- paste("ind", 1:nInd(pop), sep="")
  }
  df <- cbind(pop@ind.names, df)
  # setting the NA replacement. This doesn't work too well. 
  replacement <- ifelse(pop@type =="PA","-1","0")
  if(!quiet) cat("Writing the table to",filename,"... ")
  
  if(geo == TRUE & !is.null(pop$other[[geodf]])){
    replacemat <- matrix("", 3, 3)
    replacemat[3, 2:3] <- c("X", "Y")
    infolines <- cbind(infolines, replacemat)
    df2 <- data.frame(list("Space" = rep("", nInd(pop))))
    gdf <- as.matrix(pop@other[[geodf]])
    if(nrow(gdf) < nInd(pop)){
      gdf <- rbind(gdf, matrix("", nInd(pop) - nrow(gdf), 2))
    }
    df <- cbind(df, df2, gdf)
  }
  else if (geo == TRUE){
    popcall <- popcall[2]
    warning(paste("There is no data frame or matrix in ",
                  paste(substitute(popcall), collapse=""), 
                  "@other called ",geodf,
                  ".\nThe xy coordinates will not be represented in the resulting file.", sep=""))
  }
  
  df[df == "NA" | is.na(df)] <- replacement
  write.table(infolines, file=filename, quote=FALSE, row.names=FALSE, 
              col.names=FALSE, sep=",")
  write.table(df, file=filename, quote=TRUE, na=replacement, append=TRUE, 
              row.names=FALSE, col.names=FALSE, sep=",")
  if(!quiet) cat("Done.\n")
}
