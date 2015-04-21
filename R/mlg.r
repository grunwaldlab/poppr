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
#' Create counts, vectors, and matrices of multilocus genotypes.
#'
#' @name mlg
#'
#' @param pop a \code{\linkS4class{genind}} or \code{\linkS4class{genclone}} object.
#'
#' @param sublist a \code{vector} of population names or indices that the user
#' wishes to keep. Default to "ALL".
#'
#' @param blacklist a \code{vector} of population names or indices that the user
#' wishes to discard. Default to \code{NULL}.
#'
#' @param mlgsub a \code{vector} of multilocus genotype indices with which to
#' subset \code{mlg.table} and \code{mlg.crosspop}. NOTE: The resulting table
#' from \code{mlg.table} will only contain countries with those MLGs
#'
#' @param quiet \code{Logical}. If FALSE, progress of functions will be printed
#' to the screen. 
#'
#' @param bar \code{logical} If \code{TRUE}, a bar graph for each population
#' will be displayed showing the relative abundance of each MLG within the
#' population.
#'
#' @param indexreturn \code{logical} If \code{TRUE}, a vector will be returned
#' to index the columns of \code{mlg.table}.
#'
#' @param df \code{logical} If \code{TRUE}, return a data frame containing the
#' counts of the MLGs and what countries they are in. Useful for making graphs
#' with \code{\link{ggplot}}. 
#'
#' @param total \code{logical} If \code{TRUE}, a row containing the sum of all
#' represented MLGs is appended to the matrix produced by mlg.table.
#' 
#' @return 
#' \subsection{mlg}{
#' an integer describing the number of multilocus genotypes observed.
#' }
#' \subsection{mlg.table}{
#' a matrix with columns indicating unique multilocus genotypes and rows
#' indicating populations.
#' }
#' \subsection{mlg.vector}{
#' a numeric vector naming the multilocus genotype of each individual in
#'   the dataset.
#' }
#' \subsection{mlg.crosspop}{
#' \itemize{
#' \item{default}{ a \code{list} where each element contains a named integer vector representing the number of individuals represented from each population in that MLG}
#' \item{\code{indexreturn = TRUE}}{ a \code{vector} of integers defining the multilocus genotypes that have individuals crossing populations}
#' \item{\code{df = TRUE}}{ A long form data frame with the columns: MLG, Population, Count. Useful for graphing with ggplot2}
#' }
#' }
#' \subsection{mlg.id}{
#' a list of multilocus genotypes with the associated individual names per MLG.
#' }
#'
#' @seealso \code{\link[vegan]{diversity}} \code{\link{popsub}}
#' @author Zhian N. Kamvar
#' @examples
#'
#' # Load the data set
#' data(Aeut)
#' 
#' # Investigate the number of multilocus genotypes.
#' amlg <- mlg(Aeut)
#' amlg # 119
#' 
#' # show the multilocus genotype vector 
#' avec <- mlg.vector(Aeut)
#' avec 
#' 
#' # Get a table
#' atab <- mlg.table(Aeut, bar = FALSE)
#' atab
#' 
#' # See where multilocus genotypes cross populations
#' acrs <- mlg.crosspop(Aeut) # MLG.59: (2 inds) Athena Mt. Vernon
#' 
#' # See which individuals belong to each MLG
#' aid <- mlg.id(Aeut)
#' aid["59"] # individuals 159 and 57
#' 
#' \dontrun{
#' 
#' # A simple example. 10 individuals, 5 genotypes.
#' mat1 <- matrix(ncol=5, 25:1)
#' mat1 <- rbind(mat1, mat1)
#' mat <- matrix(nrow=10, ncol=5, paste(mat1,mat1,sep="/"))
#' mat.gid <- df2genind(mat, sep="/")
#' mlg(mat.gid)
#' mlg.vector(mat.gid)
#' mlg.table(mat.gid)
#' 
#' # Now for a more complicated example.
#' # Data set of 1903 samples of the H3N2 flu virus genotyped at 125 SNP loci.
#' data(H3N2)
#' mlg(H3N2, quiet=FALSE)
#' 
#' H.vec <- mlg.vector(H3N2)
#' 
#' # Changing the population vector to indicate the years of each epidemic.
#' pop(H3N2) <- other(H3N2)$x$country
#' H.tab <- mlg.table(H3N2, bar=FALSE, total=TRUE)
#'
#' # Show which genotypes exist accross populations in the entire dataset.
#' res <- mlg.crosspop(H3N2, quiet=FALSE)
#'
#' # Let's say we want to visualize the multilocus genotype distribution for the
#' # USA and Russia
#' mlg.table(H3N2, sublist=c("USA", "Russia"), bar=TRUE)
#' 
#' # An exercise in subsetting the output of mlg.table and mlg.vector.
#' # First, get the indices of each MLG duplicated across populations.
#' inds <- mlg.crosspop(H3N2, quiet=FALSE, indexreturn=TRUE)
#' 
#' # Since the columns of the table from mlg.table are equal to the number of
#' # MLGs, we can subset with just the columns.
#' H.sub <- H.tab[, inds]
#'
#' # We can also do the same by using the mlgsub flag.
#' H.sub <- mlg.table(H3N2, mlgsub=inds)
#'
#' # We can subset the original data set using the output of mlg.vector to
#' # analyze only the MLGs that are duplicated across populations. 
#' new.H <- H3N2[H.vec %in% inds, ]
#' 
#' }
NULL
#==============================================================================#
#' @rdname mlg
#'
#'
#' @export
#==============================================================================#

mlg <- function(pop, quiet=FALSE){
  if (!is(pop, "genlight") & !is(pop, "genind")){
    stop(paste(substitute(pop), "is not a genind or genlight object"))
  }
  if ((is.snpclone(pop) | is.genclone(pop)) && length(pop@mlg) == nrow(pop@tab)){
    out <- length(unique(pop@mlg[]))
  } else {
    if (is(pop, "genlight")) return(nInd(pop))
    if(nrow(pop@tab) == 1){
      out <- 1
    }
    else {
      out <- nrow(unique(pop@tab[, 1:ncol(pop@tab)]))
    } 
  } 
  if(quiet!=TRUE){
    cat("#############################\n")
    cat("# Number of Individuals: ", nInd(pop), "\n")
    cat("# Number of MLG: ", out, "\n")
    cat("#############################\n")
  }
  return(out)
}
#==============================================================================#
#' @rdname mlg 
#'
#' @note The resulting matrix of \code{mlg.table} can be used for analysis with 
#' the \code{\link{vegan}} package.
#' 
#' @export
#==============================================================================#
mlg.table <- function(pop, sublist="ALL", blacklist=NULL, mlgsub=NULL, bar=TRUE, 
                      total=FALSE, quiet=FALSE){  
  if (!is.genind(pop)){
    stop("This function requires a genind object.")
  }
  mlgtab <- mlg.matrix(pop)
  if (!is.null(mlgsub)){
    mlgsub <- paste("MLG", mlgsub, sep = ".")
    mlgtab <- mlgtab[, mlgsub, drop = FALSE]
    mlgtab <- mlgtab[which(rowSums(mlgtab) > 0L), , drop = FALSE]
    pop <- popsub(pop, sublist=rownames(mlgtab))
  }
  if (sublist[1] != "ALL" | !is.null(blacklist)){
    pop <- popsub(pop, sublist, blacklist)
    mlgtab <- mlgtab[unlist(vapply(popNames(pop), 
                function(x) which(rownames(mlgtab) == x), 1)), , drop=FALSE]
    rows <- rownames(mlgtab)
  }
  if (total==TRUE & (nrow(mlgtab) > 1 | !is.null(nrow(mlgtab)) )){
    mlgtab <- rbind(mlgtab, colSums(mlgtab))
    rownames(mlgtab)[nrow(mlgtab)] <- "Total"
  }

  # Dealing with the visualizations.
  if (bar){
    # If there is a population structure
    if(!is.null(popNames(pop))){
      popnames <- popNames(pop)
      if(total & nrow(mlgtab) > 1){
        popnames[length(popnames) + 1] <- "Total"
      }
      # Apply this over all populations. 
      invisible(lapply(popnames, print_mlg_barplot, mlgtab, quiet=quiet))
    } else {
      print(mlg_barplot(mlgtab) + 
        theme_classic() %+replace%
        theme(axis.text.x=element_text(size=10, angle=-45, hjust=0, vjust=1)) +
        labs(title = paste("File:", as.character(pop@call[2]), "\nN =",
                           sum(mlgtab), "MLG =", length(mlgtab))))
    }
  }
  mlgtab <- mlgtab[, which(colSums(mlgtab) > 0)]
  return(mlgtab)
}

#==============================================================================#
#' @rdname mlg
#' 
#' @param reset logical. For genclone objects, the MLGs are defined by the input
#' data, but they do not change if more or less information is added (i.e. loci are dropped). Setting \code{reset = TRUE} will recalculate MLGs. Default is \code{FALSE}, returning the MLGs defined in the @@mlg slot.
#'   
#' @note mlg.vector will recalculate the mlg vector for
#'   \code{\linkS4class{genind}} objects and will return the contents of the mlg
#'   slot in \code{\linkS4class{genclone}} objects. This means that MLGs will be
#'   different for subsetted \code{\linkS4class{genind}} objects.
#'   
#' @export
#==============================================================================#

mlg.vector <- function(pop, reset = FALSE){

  # This will return a vector indicating the multilocus genotypes.
  # note that the genotype numbers will not match up with the original numbers,
  # but will be scattered as a byproduct of the sorting. This is inconsequential
  # as the naming of the MLGs is arbitrary.
  
  # Step 1: collapse genotypes into strings
  # Step 2: sort strings and keep indices
  # Step 3: create new vector in which to place sorted index number.
  # Step 4: evaluate strings in sorted vector and increment to the respective 
  # # index vector each time a unique string occurs.
  # Step 4: Rearrange index vector with the indices from the original vector.
  if (!reset && (is.snpclone(pop) || is.genclone(pop)) && length(pop@mlg) == nrow(pop@tab)){
    return(pop@mlg[])
  }
  if (is(pop, "genlight")){
    return(seq_len(nInd(pop)))
  }
  xtab <- pop@tab
  # concatenating each genotype into one long string.
  xsort <- vapply(seq(nrow(xtab)),function(x) paste(xtab[x, ], collapse = ""), "string")
  # creating a new vector to store the counts of unique genotypes.
  countvec <- vector(length = length(xsort), mode = "integer")
  # sorting the genotypes ($x) and preserving the index ($xi). 
  xsorted <- sort(xsort, index.return=TRUE)
  
  # simple function to count number of genotypes. Num is the index number, comp
  # is the vector of genotypes.
  # This works by the "<<-" assignment operator, which is a global assignment.
  # This will search in higher environments until it finds a corresponding
  # object to modify. In this case it is countvec, which was declared above.
  
  f1 <- function(num, comp){
    if (num - 1 == 0){
      countvec[num] <<- 1L
    }
    # These have the exact same strings, thus they are the same MLG. Perpetuate
    # The MLG index.
    else if (comp[num] == comp[num - 1]){
      countvec[num] <<- countvec[num - 1]
    } else {
    # These have differnt strings, increment the MLG index by one.
      countvec[num] <<- countvec[num - 1] + 1L
    }
  }

  # applying this over all genotypes.
  lapply(1:length(xsorted$x), f1, xsorted$x)
  
  # a new vector to take in the genotype indicators
  countvec2 <- 1:length(xsort)
  
  # replacing the numbers in the vector with the genotype indicators.
  countvec2[xsorted$ix] <- countvec
  return(countvec2)
}

#==============================================================================#
#' @rdname mlg 
#' 
#' @export
#==============================================================================#

mlg.crosspop <- function(pop, sublist="ALL", blacklist=NULL, mlgsub=NULL, indexreturn=FALSE, df=FALSE, quiet=FALSE){
  if (length(sublist) == 1 & sublist[1] != "ALL" | is.null(pop(pop))){
    cat("Multiple populations are needed for this analysis.\n")
    return(0)
  }
  visible <- "original"
  if (is.genclone(pop)){
    vec <- pop@mlg[]
    if (is(pop@mlg, "MLG")){
      visible <- pop@mlg@visible
    }
  } else {
    vec <- mlg.vector(pop) 
  }
  subind <- sub_index(pop, sublist, blacklist)
  vec    <- vec[subind]
  mlgtab <- mlg.matrix(pop)

  if (!is.null(mlgsub)){
    if (visible == "custom"){
      mlgsubnames <- mlgsub
    } else {
      mlgsubnames <- paste("MLG", mlgsub, sep = ".")
    }
    matches <- mlgsubnames %in% colnames(mlgtab)

    if (!all(matches)){
      rejects <- mlgsub[!matches]
      mlgsubnames  <- mlgsubnames[matches]
      warning(mlg_sub_warning(rejects))
    }
    
    mlgtab <- mlgtab[, mlgsubnames, drop = FALSE]
    mlgs   <- setNames(1:ncol(mlgtab), colnames(mlgtab))
  
  } else {
    if (sublist[1] != "ALL" | !is.null(blacklist)){
      pop    <- popsub(pop, sublist, blacklist)
      mlgtab <- mlgtab[pop@pop.names, , drop = FALSE]
    }

    mlgs <- colSums(ifelse(mlgtab == 0L, 0L, 1L)) > 1
    if (sum(mlgs) == 0){
      cat("No multilocus genotypes were detected across populations\n")
      return(0)
    }
  }
  if (indexreturn){
    if (visible == "custom"){
      mlgout <- names(mlgs[mlgs])
    } else {
      mlgout <- unlist(strsplit(names(mlgs[mlgs]), "\\."))
      mlgout <- as.numeric(mlgout[!mlgout %in% "MLG"])        
    }
    return(mlgout)
  }
  popop <- function(x, quiet=TRUE){
    popnames <- mlgtab[mlgtab[, x] > 0L, x]
    if (length(popnames) == 1){
      names(popnames) <- rownames(mlgtab[mlgtab[, x] > 0L, x, drop=FALSE])
    }
    if (!quiet)
      cat(paste(x, ":", sep=""),paste("(",sum(popnames)," inds)", sep=""),
          names(popnames), fill=80)
    return(popnames)
  }
  # Removing any populations that are not represented by the MLGs.
  mlgtab <- mlgtab[rowSums(mlgtab[, mlgs, drop=FALSE]) > 0L, mlgs, drop=FALSE]
  # Compiling the list.
  mlg.dup <- lapply(colnames(mlgtab), popop, quiet=quiet)
  names(mlg.dup) <- colnames(mlgtab)
  if (df == TRUE){
    mlg.dup <- as.data.frame(list(MLG = rep(names(mlg.dup), sapply(mlg.dup, length)), 
                             Population = unlist(lapply(mlg.dup, names)), 
                             Count = unlist(mlg.dup)))
    rownames(mlg.dup) <- NULL
  }
  return(mlg.dup)
}



#==============================================================================#
#' @rdname mlg
#' @export
#==============================================================================#


mlg.id <- function (pop){
  if (!is.genind(pop)){
    stop(paste(substitute(pop), "is not a genind or genclone object"))
  }
  return(split(indNames(pop), mlg.vector(pop)))
}
