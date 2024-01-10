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

#' Create counts, vectors, and matrices of multilocus genotypes.
#' 
#' @name mlg
#' @md
#'   
#' @param gid a [adegenet::genind], [genclone][genclone-class], [adegenet::genlight], or [snpclone][snpclone-class] object.
#' @param strata a formula specifying the strata at which computation is to be
#'   performed.
#' @param sublist a `vector` of population names or indices that the user 
#'   wishes to keep. Default to "ALL".
#' @param exclude a `vector` of population names or indexes that the user
#' wishes to discard. Default to `NULL`.
#' @param blacklist DEPRECATED, use exclude.
#' @param mlgsub a `vector` of multilocus genotype indices with which to 
#'   subset `mlg.table` and `mlg.crosspop`. NOTE: The resulting table 
#'   from `mlg.table` will only contain countries with those MLGs
#' @param quiet `Logical`. If FALSE, progress of functions will be printed 
#'   to the screen.
#' @param bar deprecated. Same as `plot`. Retained for compatibility.
#' @param plot `logical` If `TRUE`, a bar graph for each population 
#'   will be displayed showing the relative abundance of each MLG within the 
#'   population.
#' @param indexreturn `logical` If `TRUE`, a vector will be returned 
#'   to index the columns of `mlg.table`.
#' @param df `logical` If `TRUE`, return a data frame containing the 
#'   counts of the MLGs and what countries they are in. Useful for making graphs
#'   with [ggplot].
#' @param total `logical` If `TRUE`, a row containing the sum of all 
#'   represented MLGs is appended to the matrix produced by mlg.table.
#'   
#' @return 
#' 
#' ## mlg
#'
#' an integer describing the number of multilocus genotypes observed. 
#'
#' ## mlg.table
#'
#' a matrix with columns indicating unique multilocus genotypes and rows
#' indicating populations. This table can be used with the funciton
#' [diversity_stats] to calculate the Shannon-Weaver index (H), Stoddart and
#' Taylor's index (aka inverse Simpson's index; G), Simpson's index (lambda),
#' and evenness (E5).
#' 
#' ## mlg.vector
#'
#' a numeric vector naming the multilocus genotype of each individual in the
#' dataset. 
#'
#' ## mlg.crosspop
#'
#' - **default** a `list` where each element contains a named integer
#'   vector representing the number of individuals represented from each
#'   population in that MLG 
#' - `indexreturn = TRUE` a `vector` of integers defining the multilocus
#'   genotypes that have individuals crossing populations 
#' - `df = TRUE` A long form data frame with the columns: MLG, Population,
#'   Count. Useful for graphing with ggplot2   
#'
#' ## mlg.id
#'
#' a list of multilocus genotypes with the associated individual names per MLG. 
#' 
#' @details Multilocus genotypes are the unique combination of alleles across
#'   all loci. For details of how these are calculated see `vignette("mlg",
#'   package = "poppr")`. In short, for genind and genclone objects, they are
#'   calculated by using a rank function on strings of alleles, which is
#'   sensitive to missing data. For genlight and snpclone objects, they are
#'   calculated with distance methods via [bitwise.dist] and
#'   [mlg.filter], which means that these are insensitive to missing
#'   data. Three different types of MLGs can be defined in \pkg{poppr}: 
#'
#' - **original** the default definition of multilocus genotypes as
#'   detailed above 
#' - **contracted** these are multilocus genotypes collapsed into multilocus
#'   lineages ([mll]) with genetic distance via [mlg.filter] 
#' - **custom** user-defined multilocus genotypes. These are useful for
#'   information such as mycelial compatibility groups
#'
#' **All of the functions documented here will work on any of the MLG types
#' defined in \pkg{poppr}**
#' 
#' @seealso `\link[vegan]{diversity`} 
#'  [diversity_stats]
#'  [popsub]
#'  [mll]
#'  [mlg.filter]
#'  [mll.custom]
#'  
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
#' atab <- mlg.table(Aeut, color = TRUE)
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
#' # For the mlg.table, you can also choose to display the number of MLGs across
#' # populations in the background
#' 
#' mlg.table(Aeut, background = TRUE)
#' mlg.table(Aeut, background = TRUE, color = TRUE)
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
#' mlg(H3N2, quiet = FALSE)
#' 
#' H.vec <- mlg.vector(H3N2)
#' 
#' # Changing the population vector to indicate the years of each epidemic.
#' pop(H3N2) <- other(H3N2)$x$country
#' H.tab <- mlg.table(H3N2, plot = FALSE, total = TRUE)
#' 
#' # Show which genotypes exist accross populations in the entire dataset.
#' res <- mlg.crosspop(H3N2, quiet = FALSE)
#' 
#' # Let's say we want to visualize the multilocus genotype distribution for the
#' # USA and Russia
#' mlg.table(H3N2, sublist = c("USA", "Russia"), bar=TRUE)
#' 
#' # An exercise in subsetting the output of mlg.table and mlg.vector.
#' # First, get the indices of each MLG duplicated across populations.
#' inds <- mlg.crosspop(H3N2, quiet = FALSE, indexreturn = TRUE)
#' 
#' # Since the columns of the table from mlg.table are equal to the number of
#' # MLGs, we can subset with just the columns.
#' H.sub <- H.tab[, inds]
#' 
#' # We can also do the same by using the mlgsub flag.
#' H.sub <- mlg.table(H3N2, mlgsub = inds)
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

mlg <- function(gid, quiet=FALSE){
  if (!inherits(gid, c("genlight", "genind"))) {
    stop(paste(substitute(gid), "is not a genind, or genlight object"))
  }
  if (is.clone(gid) && length(gid@mlg) == nInd(gid)) {
    out <- length(unique(gid@mlg[]))
  } else {
    if (inherits(gid, "genlight"))
      return(nInd(gid))
    if (nrow(tab(gid)) == 1) {
      out <- 1
    } else {
      out <- nrow(unique(tab(gid)[, 1:ncol(tab(gid))]))
    }
  }
  if (quiet != TRUE) {
    cat("#############################\n")
    cat("# Number of Individuals: ", nInd(gid), "\n")
    cat("# Number of MLG: ", out, "\n")
    cat("#############################\n")
  }
  return(out)
}
#==============================================================================#
#' @rdname mlg 
#' 
#' @param color an option to display a single barchart for mlg.table, colored by
#'   population (note, this becomes facetted if `background = TRUE`).
#'   
#' @param background an option to display the the total number of MLGs across
#'   populations per facet in the background of the plot.
#'
#' @note The resulting matrix of `mlg.table` can be used for analysis with 
#' the \pkg{vegan} package.
#' 
#' @export
#==============================================================================#
mlg.table <- function(gid, strata = NULL, sublist = "ALL", exclude = NULL, blacklist = NULL, 
                      mlgsub = NULL, bar = TRUE, plot = TRUE, total = FALSE, 
                      color = FALSE, background = FALSE, quiet = FALSE){  
  if (!is.genind(gid) & !is(gid, "snpclone")){
    stop("This function requires a genind object.")
  }
  the_call <- match.call()
  if ("bar" %in% names(the_call)){
    plot <- bar
    warning("In poppr version 2.0, bar is deprecated. Please use plot.")
  }
  if (!is.null(blacklist)) {
    warning(
      option_deprecated(
        the_call, 
        "blacklist", 
        "exclude", 
        "2.8.7.", 
        "Please use `exclude` in the future"
       ), 
      immediate. = TRUE
    )
    exclude <- blacklist
  }
  the_data <- utils::capture.output(the_call[["gid"]])
  if (!is.null(strata)){
    setPop(gid) <- strata
  }
  mlgtab <- mlg.matrix(gid)
  if (!is.null(mlgsub)){
    if (is.numeric(mlgsub)){
      mlgsub <- paste("MLG", mlgsub, sep = ".")
    }
    mlgtab <- mlgtab[, mlgsub, drop = FALSE]
    mlgtab <- mlgtab[which(rowSums(mlgtab) > 0L), , drop = FALSE]
    gid <- popsub(gid, sublist = rownames(mlgtab))
  }
  if (sublist[1] != "ALL" | !is.null(exclude)){
    gid <- popsub(gid, sublist, exclude)
    mlgtab <- mlgtab[popNames(gid), , drop=FALSE]
    rows <- rownames(mlgtab)
  }
  if (total == TRUE && nrow(mlgtab) > 1){
    mlgtab <- rbind(mlgtab, colSums(mlgtab))
    rownames(mlgtab)[nrow(mlgtab)] <- "Total"
  }

  # Dealing with the visualizations.
  if (plot){
    # If there is a population structure
    if(!is.null(popNames(gid))){
      popnames <- popNames(gid)
      if(total & nrow(mlgtab) > 1){
        popnames[length(popnames) + 1] <- "Total"
      }
    }
    if (total && nrow(mlgtab) > 1 && (color | background)){
      ggmlg <- mlg_barplot(mlgtab[-nrow(mlgtab), , drop = FALSE], 
                           color = color, background = background)
    } else {
      ggmlg <- mlg_barplot(mlgtab, color = color, background = background)
    }
    print(ggmlg + 
          myTheme + 
          labs(title = paste("Data:", the_data, "\nN =",
                             sum(mlgtab), "MLG =", ncol(mlgtab))))
  }
  mlgtab <- mlgtab[, which(colSums(mlgtab) > 0), drop = FALSE]
  return(mlgtab)
}

#' @rdname mlg
#' 
#' @param reset logical. For genclone objects, the MLGs are defined by the input
#'   data, but they do not change if more or less information is added (i.e.
#'   loci are dropped). Setting `reset = TRUE` will recalculate MLGs.
#'   Default is `FALSE`, returning the MLGs defined in the @@mlg slot.
#'   
#' @note mlg.vector will recalculate the mlg vector for
#'   [adegenet::genind] objects and will return the contents of the mlg
#'   slot in [genclone][genclone-class] objects. This means that MLGs will be
#'   different for subsetted [adegenet::genind] objects.
#'   
#' @export
mlg.vector <- function(gid, reset = FALSE){

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
  if (!reset && is.clone(gid) && length(gid@mlg) == nInd(gid)){
    return(gid@mlg[])
  }
  if (inherits(gid, "genlight")){
    if (is.clone(gid)) gid@mlg <- seq_len(nInd(gid))
    return(mlg.filter(gid, 
                      threshold = .Machine$double.eps ^ 0.5,
                      distance = bitwise.dist,
                      missing_match = TRUE))
  } 

  xtab <- gid@tab

  # concatenating each genotype into one long string.
  xsort <- vapply(seq(nrow(xtab)),function(x) paste(xtab[x, ], collapse = ""), "string")

  # Interestingly enough, the methods below are only about 1% slower than using
  # rank. Both are effectively equivalent, except that rank gives the rank of
  # the order they were found. Since there is not much of a speedup and it would
  # change the names of things, I decided not to change this.
  # return(rank(xsort, ties.method = "min"))

  # creating a new vector to store the counts of unique genotypes.
  countvec <- countvec2 <- vector(length = length(xsort), mode = "integer")

  # sorting the genotypes ($x) and preserving the index ($xi). 
  xsorted <- sort(xsort, index.return=TRUE)
  
  # for loop to count number of genotypes. Num is the index number, comp
  # is the vector of genotypes.
  sorted_genos <- xsorted$x

  for (num in seq_along(sorted_genos)) {
    if (num - 1 == 0){
      countvec[num] <- 1L
    }
    # These have the exact same strings, thus they are the same MLG. Perpetuate
    # The MLG index.
    else if (sorted_genos[num] == sorted_genos[num - 1]){
      countvec[num] <- countvec[num - 1]
    } else {
    # These have differnt strings, increment the MLG index by one.
      countvec[num] <- countvec[num - 1] + 1L
    }
  }

  # replacing the numbers in the vector with the genotype indicators.
  countvec2[xsorted$ix] <- countvec
  return(countvec2)
}







#' @rdname mlg
#' @export
mlg.crosspop <- function(gid, strata = NULL, sublist = "ALL", exclude = NULL, blacklist = NULL, 
                         mlgsub = NULL,
                         indexreturn = FALSE, df = FALSE, quiet = FALSE){

  if (!is.null(blacklist)) {
    warning(
      option_deprecated(
        match.call(), 
        "blacklist", 
        "exclude", 
        "2.8.7.", 
        "Please use `exclude` in the future"
       ), 
      immediate. = TRUE
    )
    exclude <- blacklist
  }
  insufficient_subset <- length(sublist) == 1 && sublist[1] != "ALL"
  insufficient_pop    <- is.null(pop(gid)) || nPop(gid) == 1
  if (insufficient_subset | insufficient_pop){
    stop("Multiple populations are needed for this analysis.\n")
  }
  visible <- "original"
  if (is.genclone(gid) | is(gid, "snpclone")){
    vec <- gid@mlg[]
    if (is(gid@mlg, "MLG")){
      visible <- visible(gid@mlg)
    }
  } else {
    vec <- mlg.vector(gid) 
  }
  if (!is.null(strata)){
    setPop(gid) <- strata
  }
  subind <- sub_index(gid, sublist, exclude)
  vec    <- vec[subind]
  mlgtab <- mlg.matrix(gid)

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
    mlgs   <- stats::setNames(1:ncol(mlgtab), colnames(mlgtab))
  
  } else {
    if (sublist[1] != "ALL" | !is.null(exclude)){
      gid    <- popsub(gid, sublist, exclude)
      mlgtab <- mlgtab[popNames(gid), , drop = FALSE]
    }

    mlgs <- colSums(ifelse(mlgtab == 0L, 0L, 1L)) > 1
    if (sum(mlgs) == 0){
      message("No multilocus genotypes were detected across populations\n")
      return(NULL)
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



#' @rdname mlg
#' @export
mlg.id <- function (gid){
  if (!is.genind(gid) & !is(gid, "snpclone")){
    stop(paste(substitute(gid), "is not a genind or genclone object"))
  }
  return(split(indNames(gid), mlg.vector(gid)))
}
