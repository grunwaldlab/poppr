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



# Keeping this for compatability
new.read.genalex <- function(genalex, ploidy=2, geo=FALSE, region=FALSE){
  return(read.genalex(genalex, ploidy, geo, region))
}




#==============================================================================#
# hier = a nested formula such as ~ A/B/C where C is nested within B, which is
# nested within A.
#
# df = a data frame containing columns corresponding to the variables in hier.
#
# example:
# df <- data.frame(list(a = letters, b = LETTERS, c = 1:26))
# newdf <- make_hierarchy(~ a/b/c, df)
# df[names(newdf)] <- newdf # Add new columns.
#==============================================================================#
make_hierarchy <- function(hier, df, expand_label = FALSE){
  newlevs <- attr(terms(hier), "term.labels")
  levs <- all.vars(hier)
  if (length(levs) > 1){
    newlevs <- gsub(":", "_", newlevs)
  }
  if (!all(levs %in% names(df))){
    stop(hier_incompatible_warning(levs, df))
  }
  newdf <- df[levs[1]]
  if (!expand_label){
    newlevs <- levs
  }
  lapply(1:length(levs), function(x) newdf[[newlevs[x]]] <<- as.factor(pop_combiner(df, levs[1:x])))
  return(newdf)
}


#==============================================================================#
# Function for creating the structure data frame needed for ade4's AMOVA
# implementation.
#==============================================================================#
make_ade_df <- function(hier, df, expanded = FALSE){
  if (expanded){
    levs <- attr(terms(hier), "term.labels")
  } else {
    levs <- all.vars(hier)
  }
  if(length(levs) <= 1){
    # stop("Only one level present")
    return(NULL)
  }
  levs <- gsub(":", "_", levs)
  if(!all(levs %in% names(df))){
    stop(hier_incompatible_warning(levs, df))
  }
  smallest  <- df[[levs[length(levs)]]]
  smallinds <- !duplicated(smallest)
  newdf     <- df[smallinds, ]
  newdf     <- newdf[-length(levs)]
  if (length(newdf) > 1){
    factlist <- lapply(newdf, function(x) factor(x, unique(x)))
  } else {
    factlist        <- list(factor(newdf[[1]], unique(newdf[[1]])))
    names(factlist) <- names(newdf)
  }  
  return(rev(data.frame(factlist)))
}

#==============================================================================#
# Haplotype pooling. 
# The following functions are necessary to account for within sample variation. 
# They will separate the haplotypes of a genind object and repool them so that
# there are n*k individuals in the new data set where n is the number of
# individuals and k is the ploidy. 
#==============================================================================#
## Main Function. Lengthens the population hierarchy as well.
pool_haplotypes <- function(x, dfname = "population_hierarchy"){
  ploidy        <- ploidy(x)
  df            <- other(x)[[dfname]]
  df$Individual <- indNames(x)
  df            <- df[rep(1:nrow(df), ploidy), ]
  newx          <- repool(separate_haplotypes(x))
  pop(newx)     <- df$Individual
  other(newx)[[dfname]] <- df
  return(newx)
}
## Secondary function. Transforms the genind object into a loci object and
## collects information necessary to collect individual haplotypes.
separate_haplotypes <- function(x){
  ploidy        <- ploidy(x)
  allele_list   <- unlist(lapply(x@all.names, nchar))
  allele_length <- sum(allele_list)/length(allele_list)
  if (!all(allele_list == allele_length)){
    stop("not all alleles are of equal length.")
  }
  inds        <- indNames(x)
  x.loc       <- as.loci(x)
  loci_cols   <- attr(x.loc, "locicol")
  pop_col     <- which(names(x.loc) == "population")
  geno_length <- allele_length*ploidy + ploidy - 1
  sep         <- which(1:geno_length %% (allele_length + 1) == 0)
  sep         <- c(0, sep, geno_length + 1)
  hap_fac     <- rep(1:ploidy, each = allele_length + 1)
  allele_inds <- lapply(1:ploidy, function(i) c(sep[i] + 1, sep[i + 1] - 1))
  haplist     <- lapply(allele_inds, hap2genind, x.loc, loci_cols, inds, hap_fac)
  return(haplist)
}

## Obtains specific haplotypes from a genind object.
## 
## Arguments:
##   inds      - indexes of the first and last characters defining the haplotype
##   x         - the loci object
##   loci_cols - the columns defining the loci
##   ind_names - the sample names
##   hap_fac   - a vector defining the haplotype. 
##
## Since the inds and hap_fac args are not intuitive, I should demonstrate.
## Let's say we have a two individuals with two loci:
##      1      2
##   A: 10/20  30/40
##   B: 15/15  25/10
##
## They have the haplotypes of 
##   A1: 10 30 
##   A2: 20 40 
##   B1: 15 25
##   B2: 15 10
##
## Each genotype at a locus is 5 characters long. To get the haplotypes, we need
## to avoid the separator.
## inds in this case will be c(1, 2) for the first haplotype and c(4, 5) for the
## second haplotype and hap_fac will be c(1,1,1,2,2,2) for both.
##
hap2genind <- function(inds, x, loci_cols, ind_names, hap_fac){
  new_ind_names <- paste(ind_names, hap_fac[inds[2]], sep = ".")
  new_pop       <- rep(hap_fac[inds[2]], nrow(x))
  x2 <- lapply(x[loci_cols], substr, inds[1], inds[2])
  x2 <- data.frame(x2, stringsAsFactors = F)
  x2 <- df2genind(x2, ploidy = 1, pop=new_pop, ind.names=new_ind_names)
  return(x2)
}

## Function for determining if a genind object has any heterozygous sites.
check_Hs <- function(x){
  res <- any(x@tab > 0 & x@tab < 1, na.rm = TRUE)
  return(res)
}
#==============================================================================#
# Implementation of ade4's AMOVA function. Note that this cannot be used at the
# moment to calculate within individual variances. It will either compute a
# distance matrix or take in a distance matrix. Missing data must be treated
# here and is currently treated as extra alleles, but this can be modified by
# the user. Since ade4 needs a euclidean matrix, this function will, by default,
# correct the distance via the cailliez correction, which will add a number to 
# all the distances to satisfy euclidean nature. 
# Clone correction at the lowest level of the hierarchy is possible. 
#
# Note: This takes a nested formula argument. If you want to analyze the
# hierarchy of Year to Population to Subpopulation, you should make sure you
# have the right data frame in your "other" slot and then write the formula
# thusly: ~ Year/Population/Subpopulation
#
# arguments:
#
# x            = a genind object
# hier         = a formula such as ~Pop/Subpop
# clonecorrect = This refers to clone correction of the final output relative to
#                The lowest hierararchical level. FALSE
# within       = should witin individual variation be calculated? TRUE
# dist         = A user provided distance matrix NULL
# squared      = Is the distance matrix squared? TRUE
# correction   = A correction for non-euclidean distances provided by ade4.
#                The default, "quasieuclid", seems to give the best results.
# dfname       = the data frame containing the population hierarchy.
#                "population_hierarchy"
# sep          = the separator for the population hierarchy levels. "_"
# missing      = how to deal with missing data. Default is "loci".
# cutoff       = a cutoff for percent missing data to tolerate. 0.05
# quiet        = Should messages be printed? TRUE
#==============================================================================#
#' @importFrom ade4 amova is.euclid cailliez quasieuclid lingoes
ade4_amova <- function(x, hier = NULL, clonecorrect = FALSE, within = TRUE, dist = NULL, 
                       squared = TRUE, correction = "quasieuclid", 
                       dfname = "population_hierarchy", sep = "_", missing = "loci", 
                       cutoff = 0.05, quiet = FALSE){
  if (!is.genind(x)) stop(paste(substitute(x), "must be a genind object."))
  if (is.null(hier)) stop("A population hierarchy must be specified")
  parsed_hier <- gsub(":", sep, attr(terms(hier), "term.labels"))
  full_hier <- parsed_hier[length(parsed_hier)]
  
  if (is.genclone(x)){
    setpop(x) <- hier
    other(x)[[dfname]] <- gethierarchy(x, hier, combine = FALSE)
  } else {
    if (!dfname %in% names(other(x))){
      stop(paste(dfname, "is not present in the 'other' slot"))
    }
    if (!full_hier %in% names(other(x)[[dfname]])){
      hiers <- all.vars(hier)
      if (!all(hiers %in% names(other(x)[[dfname]]))){
        hier_incompatible_warning(hiers, df)
      }
      x <- splitcombine(x, hier = hiers, dfname = dfname, method = 2)
    } else {
      pop(x) <- other(x)[[dfname]][[full_hier]]
    }
  }
  # Treat missing data. This is a part I do not particularly like. The distance
  # matrix must be euclidean, but the dissimilarity distance will not allow
  # missing data to contribute to the distance. 
  #
  # This is corrected by setting missing to zero: more diversity
  # missing to mean of the columns: indiscreet distances.
  # remove loci at cutoff
  # remove individuals at cutoff
  if (clonecorrect){
    x <- clonecorrect(x, hier = hier, keep = 1:length(all.vars(hier)))
  }
  if (within & ploidy(x) > 1 & check_Hs(x)){
    hier <- update(hier, ~./Individual)
    x    <- pool_haplotypes(x, dfname = dfname)
  }
  x       <- missingno(x, type = missing, cutoff = cutoff, quiet = quiet)
  hierdf  <- make_hierarchy(hier, other(x)[[dfname]])
  xstruct <- make_ade_df(hier, hierdf)
  if (is.null(dist)){
    xdist <- sqrt(diss.dist(clonecorrect(x, hier = NA), frac = FALSE))
  } else {
    datalength <- choose(nInd(x), 2)
    mlgs       <- mlg(x, quiet = TRUE)
    mlglength  <- choose(mlgs, 2)
    if (length(dist) > mlglength & length(dist) == datalength){
      corrected <- .clonecorrector(x)
      xdist     <- as.dist(as.matrix(dist)[corrected, corrected])
    } else if(length(dist) == mlglength){
      xdist <- dist
    } else {
      msg <- paste("\nDistance matrix does not match the data.",
      "\nUncorrected observations expected..........", nInd(x),
      "\nClone corrected observations expected......", mlgs,
      "\nObservations in provided distance matrix...", ceiling(sqrt(length(dist)*2)),
      ifelse(within == TRUE, "\nTry setting within = FALSE.", "\n"))
      stop(msg)
    }
    if (squared){
      xdist <- sqrt(xdist)
    }
  }
  if (!is.euclid(xdist)){
    CORRECTIONS <- c("cailliez", "quasieuclid", "lingoes")
    try(correct <- match.arg(correction, CORRECTIONS))
    if (!exists("correct")){
      stop(not_euclid_msg(correction))
    } else {
      correct_fun <- match.fun(correct)
      if (correct == CORRECTIONS[2]){
        cat("Distance matrix is non-euclidean.\n")
        cat("Utilizing quasieuclid correction method. See ?quasieuclid for details.\n")
        xdist <- correct_fun(xdist)        
      } else {
        xdist <- correct_fun(xdist, print = TRUE, cor.zero = FALSE)        
      }
    }
  }
  xtab  <- t(mlg.matrix(x))
  xtab    <- as.data.frame(xtab[unique(mlg.vector(x)), ])
  return(ade4::amova(samples = xtab, distances = xdist, structures = xstruct))
}


#==============================================================================#
# Create a population hierarchy in a genind object if one is not there.
#==============================================================================#
genind_hierarchy <- function(x, df = NULL, dfname = "population_hierarchy"){
  if (is.null(df)){
    df <- data.frame(Pop = pop(x))
  } 
  other(x)[[dfname]] <- df
  return(x)
}


#==============================================================================#
# Idea, create a table representing missing data in the data set. This will 
# show missing data across loci and populations. 
# Args:
# 
#   x - a genind object
#   percent - should the data be presented as a percent or raw value? Note that 
#             this represents percentage with respect to the population. 
#   plot - should this be plotted
#   df - return a data frame in long format
#   returnplot - return the plot as the second element of a list
#   low - color representing no missing data
#   high - color representing missing data
#   plotlab - print the missing values on the plot.
#   scaled - if percent = TRUE, then this will make the scale from 0 to 100, if 
#   this is not turned on (default), then the highest missing value will be the max.
#==============================================================================#

# tabulate the amount of missing data per locus. 
number_missing_locus <- function(x, divisor){
  missing_result <- colSums(1 - propTyped(x, by = "both"))
  return(missing_result/divisor)
}

# tabulate the amount of missing data per genotype. 
number_missing_geno <- function(x, divisor){
  missing_result <- rowSums(1 - propTyped(x, by = "both"))
  return(missing_result/divisor)
}

# main function
missing_table <- function(x, percent = TRUE, plot = FALSE, df = FALSE, 
                          returnplot = FALSE, low = "blue", high = "red", 
                          plotlab = TRUE, scaled = TRUE){
  pops <- seppop(x, drop = FALSE)
  pops$Total <- x
  inds <- 1
  if (percent){
    inds <- c(table(pop(x)), nInd(x))
  }
  misstab <- matrix(0, nrow = nLoc(x) + 1, ncol = length(pops))
  misstab[1:nLoc(x), ] <- vapply(pops, number_missing_locus, numeric(nLoc(x)), 1)
  misstab[-nrow(misstab), ] <- t(apply(misstab[-nrow(misstab), ], 1, "/", inds))
  misstab[nrow(misstab), ] <- colMeans(misstab[-nrow(misstab), ])
  rownames(misstab) <- c(x@loc.names, "Mean")
  colnames(misstab) <- names(pops)
  if (all(misstab == 0)){
    cat("No Missing Data Found!")
    return(NULL)
  }
  if (plot){
    
    missdf <- melt(misstab, varnames = c("Locus", "Population"), value.name = "Missing")
    leg_title <- "Missing"
    missdf[1:2] <- data.frame(lapply(missdf[1:2], function(x) factor(x, levels = unique(x))))
    # levels(missdf[[1]]) <- rev(levels(missdf[[1]]))
    plotdf <- textdf <- missdf
    if (percent) {
      plotdf$Missing <- round(plotdf$Missing*100, 3)
      textdf$Missing <- paste(plotdf$Missing, "%")
      miss <- "0 %"
      title <- paste("Percent missing data per locus and population of", 
                     as.character(substitute(x)))
      leg_title <- paste("Percent", leg_title)
      if(!scaled) lims <- c(0, 100)
    } else {
      textdf$Missing <- round(textdf$Missing, 3)
      miss <- 0
      title <- paste("Missing data per locus and population of", 
                     as.character(substitute(x)))
    }
    if (scaled | !percent){
      lims <- c(0, max(plotdf$Missing))
    }
    linedata <- data.frame(list(yint = ncol(misstab) - 0.5, 
                                xint = nrow(misstab) - 0.5))
    textdf$Missing <- ifelse(textdf$Missing == miss, "", textdf$Missing)
    plotdf$Missing[plotdf$Locus == "Mean" & plotdf$Population == "Total"] <- NA
    outplot <- ggplot(plotdf, aes_string(x = "Locus", y = "Population")) + 
               geom_tile(aes_string(fill = "Missing")) +
               labs(list(title = title, x = "Locus", y = "Population")) +
               labs(fill = leg_title) + 
               scale_fill_gradient(low = low, high = high, na.value = "white", 
                                   limits = lims) +
               geom_hline(aes_string(yintercept = "yint"), data = linedata) + 
               geom_vline(aes_string(xintercept = "xint"), data = linedata) 
    if (plotlab){
      outplot <- outplot + geom_text(aes_string(label = "Missing"), 
                                     data = textdf)
    }
    outplot <- outplot +
               theme_classic() + 
               scale_x_discrete(expand = c(0, -1)) + 
               scale_y_discrete(expand = c(0, -1)) + 
               theme(axis.text.x = element_text(size = 10, angle = -45, 
                                                hjust = 0, vjust = 1))
    print(outplot)
  }
  if (df){
    if(!exists("missdf")){
      missdf <- melt(misstab, varnames = c("Locus", "Population"), 
                     value.name = "Missing")
    }
    misstab <- missdf
  } else {
    attr(misstab, "dimnames") <- list(Locus = rownames(misstab), 
                                      Population = colnames(misstab))
    misstab <- t(misstab)
  }
  if (returnplot){
    misstab <- list(table = misstab, plot = outplot)
  }
  return(misstab)
}


#==============================================================================#
# This function will simply return a table denoting private alleles within
# populations. The rows will indicate populations with private alleles and the
# columns will be the specific alleles. 
#==============================================================================#
private_alleles <- function(gid){
  if (!is.genind(gid) & !is.genpop(gid)) stop(paste(gid, "is not a genind or genpop object."))
    if (is.genind(gid) & !is.null(pop(gid)) | is.genpop(gid) & nrow(gid@tab) > 1){
      if (is.genind(gid)){
          gid.pop <- truenames(genind2genpop(gid, quiet = TRUE))
        } else {
          gid.pop <- truenames(gid)
        }
        privates <- gid.pop[, colSums(ifelse(gid.pop > 0, 1, 0), na.rm = TRUE) < 2]
        privates <- privates[rowSums(privates) > 0, ]
        return(privates)
      } else {
    stop("There are no populations detected")
  }
}

#==============================================================================#
# The following two functions serve as summary functions across loci. It was
# noted that there are no functions that really give explicit summaries
# accross loci, so these were written to do so. 
#
# The function locus_table_pegas is the internal workhorse. It will process a 
# summary.loci object into a nice table utilizing the various diversity indices 
# provided by vegan. Note that it has a catch which will remove all allele names
# that are any amount of zeroes and nothing else. The reason for this being that
# these alleles actually represent missing data. 
# The wrapper for this is locus_table, which will take in a genind object and 
# send it into locus_table_pegas. 
# Note: lev argument has only the options of "allele" or "genotype"
#==============================================================================#

locus_table_pegas <- function(x, index = "simpson", lev = "allele", type = "codom"){
  unique_types <- x[[lev]]
  # Removing any zero-typed alleles that would be present with polyploids.
  zero_names   <- grep("^0+?$", names(unique_types))
  
  if (length(zero_names) > 0 & type == "codom"){
    unique_types <- unique_types[-zero_names]
  }
  
  N       <- length(unique_types)
  H       <- vegan::diversity(unique_types)
  G       <- vegan::diversity(unique_types, "inv")
  Simp    <- vegan::diversity(unique_types, "simp")
  nei     <- (N/(N-1)) * Simp
  
  if (index == "simpson"){
    idx        <- Simp
    names(idx) <- "1-D"
  } else if (index == "shannon"){
    idx        <- H
    names(idx) <- "H"
  } else {
    idx        <- G
    names(idx) <- "G"
  }
  
  E.5        <- (G - 1)/(exp(H) - 1)
  names(N)   <- lev
  return(c(N, idx, Hexp = nei, Evenness = E.5))
}

locus_table <- function(x, index ="simpson", lev = "allele", population = "ALL", information = TRUE){
  INDICES <- c("shannon", "simpson", "invsimpson")
  index   <- match.arg(index, INDICES)
  x       <- popsub(x, population, drop = FALSE)
  x.loc   <- summary(as.loci(x))
  outmat  <- vapply(x.loc, locus_table_pegas, numeric(4), index, lev, x@type)
  loci    <- colnames(outmat)
  divs    <- rownames(outmat)
  outmat  <- t(outmat)
  dimlist <- list(`locus` = loci, `summary` = divs)
  attr(outmat, "dimnames") <- dimlist
  if (information){
    if (index == "simpson"){
      msg <- "Simpson index"
    } else if (index == "shannon"){
      msg <- "Shannon-Wiener index"
    } else {
      msg <- "Stoddard and Taylor index"
    }
    cat("\n", divs[1], "= Number of observed", paste0(divs[1], "s"))
    cat("\n", divs[2], "=", msg)
    cat("\n", divs[3], "= Nei's 1978 expected heterozygosity\n")
    cat("------------------------------------------\n")
  }
  return(outmat)
}


pair_ia <- function(pop){

  if(pop@type == "codom"){
    pop_loci <- seploc(pop)
    loci_pairs <- combn(pop@loc.names, 2)
    pair_ia_vector <- apply(loci_pairs, 2, function(x) .Ia.Rd(pop_loci[x]))
    colnames(pair_ia_vector) <- apply(loci_pairs, 2, paste, collapse = ":")
  } else {
    loci_pairs <- combn(1:nLoc(pop), 2)
    pair_ia_vector <- apply(loci_pairs, 2, function(x) .PA.Ia.Rd(pop[, x], missing = "ignore"))
    colnames(pair_ia_vector) <- apply(combn(pop@loc.names, 2), 2, paste, collapse = ":")
  }
  rownames(pair_ia_vector) <- c("Ia", "rbarD")    
  return(pair_ia_vector)
}


poppr_pair_ia <- function(pop){
  if(is.null(pop(pop))){
    return(pair_ia(pop))
  }
  pops <- seppop(pop, drop = FALSE)
  loci_pairs <- choose(nLoc(pop), 2)
  res_mat <- matrix(0.5, 2, loci_pairs)
  pops_array <- vapply(pops, pair_ia, res_mat)
  return(pops_array)
}


testing_funk <- function(){
  cat("This test worked...maybe.\n")
}





new.poppr <- function(pop,total=TRUE,sublist=c("ALL"),blacklist=c(NULL), sample=0,
  method=1,missing="ignore", quiet=FALSE,clonecorrect=FALSE,hier=c(1),dfname="hier",
  hist=TRUE, minsamp=10){
  METHODS = c("multilocus", "permute alleles", "parametric bootstrap",
      "non-parametric bootstrap")
	x <- .file.type(pop, missing=missing, clonecorrect=clonecorrect, hier=hier, 
                  dfname=dfname, quiet=TRUE)	
  # The namelist will contain information such as the filename and population
  # names so that they can easily be ported around.
  namelist <- NULL
  callpop <- match.call()
  if(!is.na(grep("system.file", callpop)[1])){
    popsplt <- unlist(strsplit(pop, "/"))
    namelist$File <- popsplt[length(popsplt)]
  }
  else if(is.genind(pop)){
    namelist$File <- x$X
  }
  else{
    namelist$File <- basename(x$X)
  }
  #poplist <- x$POPLIST
  pop <- popsub(x$GENIND, sublist=sublist, blacklist=blacklist)
  poplist <- .pop.divide(pop)
  # Creating the genotype matrix for vegan's diversity analysis.
  pop.mat <- mlg.matrix(pop)
  if (total==TRUE & !is.null(poplist)){
    poplist$Total <- pop
    pop.mat <- rbind(pop.mat, colSums(pop.mat))
  }
  sublist <- names(poplist)
  Iout <- NULL
  result <- NULL
  origpop <- x$GENIND
  rm(x)
  total <- toupper(total)
  missing <- toupper(missing)
  type <- pop@type
  # For presence/absences markers, a different algorithm is applied. 
  if(type=="PA"){
    .Ia.Rd <- .PA.Ia.Rd
  }
  #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
  # Creating an indicator for multiple subpopulations.
  # MPI = Multiple Population Indicator
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
  if (is.null(poplist)){
    MPI <- NULL
  }
  else{
    MPI <- 1
  }

  #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
  # ANALYSIS OF MULTIPLE POPULATIONS.
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,# 
 
	if (!is.null(MPI)){
    
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
    # Calculations start here.
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#

    MLG.vec <- vapply(sublist, function(x) mlg(poplist[[x]], quiet=TRUE), 1)
    N.vec <- vapply(sublist, function(x) length(poplist[[x]]@ind.names), 1)
    # Shannon-Weiner vegan::diversity index.
    H <- vegan::diversity(pop.mat)
    # E_1, Pielou's evenness.
    # J <- H / log(rowSums(pop.mat > 0))
    # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
    G <- vegan::diversity(pop.mat, "inv")
    Hexp <- (N.vec/(N.vec-1))*vegan::diversity(pop.mat, "simp")
    # E_5
    E.5 <- (G-1)/(exp(H)-1)
    # rarefaction giving the standard errors. This will use the minimum pop size
    # above a user-defined threshold.
    raremax <- ifelse(is.null(nrow(pop.mat)), sum(pop.mat), 
                      ifelse(min(rowSums(pop.mat)) > minsamp, 
                             min(rowSums(pop.mat)), minsamp))
    
    N.rare <- rarefy(pop.mat, raremax, se=TRUE)
    IaList <- NULL
    invisible(lapply(sublist, function(x) 
                              IaList <<- rbind(IaList, 
                                               .new.ia(poplist[[x]], 
                                                       sample=sample, 
                                                       method=method, 
                                                       quiet=quiet, 
                                                       missing=missing, 
                                                       namelist=list(File=namelist$File, population = x),
                                                       hist=hist
              ))))

    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
    # Making the data look pretty.
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
    Iout <- as.data.frame(list(Pop=sublist, N=N.vec, MLG=MLG.vec, 
                                eMLG=round(N.rare[1, ], 3), 
                                SE=round(N.rare[2, ], 3), 
                                H=round(H, 3), 
                                G=round(G,3),
                                Hexp=round(Hexp, 3),
                                E.5=round(E.5,3),
                                IaList,
                                File=namelist$File))
    rownames(Iout) <- NULL
    return(final(Iout, result))
	}
  #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
  # ANALYSIS OF SINGLE POPULATION. This is for if there are no subpopulations to
  # be analyzed. For details of the functions utilized, see the notes above.
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
  else { 
    MLG.vec <- mlg(pop, quiet=TRUE)
    N.vec <- length(pop@ind.names)
    # Shannon-Weiner vegan::diversity index.
    H <- vegan::diversity(pop.mat)
    # E_1, Pielou's evenness.
    # J <- H / log(rowSums(pop.mat > 0))
    # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
    G <- vegan::diversity(pop.mat, "inv")
    Hexp <- (N.vec/(N.vec-1))*vegan::diversity(pop.mat, "simp")
    # E_5
    E.5 <- (G-1)/(exp(H)-1)
    # rarefaction giving the standard errors. No population structure means that
    # the sample is equal to the number of individuals.
    N.rare <- rarefy(pop.mat, sum(pop.mat), se=TRUE)
    IaList <- .new.ia(pop, sample=sample, method=method, quiet=quiet, missing=missing,
                      namelist=(list(File=namelist$File, population="Total")),
                      hist=hist)
    Iout <- as.data.frame(list(Pop="Total", N=N.vec, MLG=MLG.vec, 
                          eMLG=round(N.rare[1, ], 3), 
                          SE=round(N.rare[2, ], 3),
                          H=round(H, 3), 
                          G=round(G,3), 
                          Hexp=round(Hexp, 3), 
                          E.5=round(E.5,3), 
                          as.data.frame(t(IaList)),
                          File=namelist$File))
    rownames(Iout) <- NULL
    return(final(Iout, result))
  }
}


new.poppr.all <- function(filelist, ...) {
	result <- NULL
	for(a in filelist){
    cat("| File: ",basename(a),"\n")
		result <- rbind(result, new.poppr(a, ...))
	}
	result <- extract.info(result)
	return(result)
}

.new.ia <- function(pop,sample=0,method=1,quiet=FALSE,namelist=NULL,missing="ignore",
                    hist=TRUE){
  METHODS = c("multilocus", "permute alleles", "parametric bootstrap",
      "non-parametric bootstrap")
  if(pop@type!="PA"){
    type <- pop@type
    popx <- seploc(pop)
  }
  else {
    type <- pop@type
    popx <- pop
    .Ia.Rd <- .PA.Ia.Rd
  }

  # if there are less than three individuals in the population, the calculation
  # does not proceed. 
  if (nInd(pop) < 3){
    IarD <- as.numeric(c(NA,NA))
    names(IarD) <- c("Ia", "rbarD")
    if(sample==0){
      return(IarD)
    }
    else{
      IarD <- as.numeric(rep(NA,4))
      names(IarD) <- c("Ia","p.Ia","rbarD","p.rD")
      return(IarD)
    }
  }
	IarD <- .Ia.Rd(popx, missing)
  # data vomit options.
  if (!quiet){
    cat("|", namelist$population ,"\n")
  }
  names(IarD) <- c("Ia", "rbarD")
  # no sampling, it will simply return two named numbers.
  if (sample==0){
    Iout <- IarD
    result <- NULL
  }
  # sampling will perform the iterations and then return a data frame indicating
  # the population, index, observed value, and p-value. It will also produce a 
  # histogram.
  else{
    Iout <- NULL 
    idx <- as.data.frame(list(Index=names(IarD)))
    samp <- .sampling(popx, sample, missing, quiet=quiet, type=type, method=method)
    sampstatIa <- c(range(samp[[1]]), median(samp[[1]]), mean(samp[[1]]), var(samp[[1]]), sd(samp[[1]]))
    names(sampstatIa) <- c("Ia.min", "Ia.max", "Ia.med", "Ia.mean", "Ia.var", "Ia.sd")
    sampstatRd <- c(range(samp[[2]]), median(samp[[2]]), mean(samp[[2]]), var(samp[[2]]), sd(samp[[2]]))
    names(sampstatRd) <- c("Rd.min", "Rd.max", "Rd.med", "Rd.mean", "Rd.var", "Rd.sd")
    samp2 <- rbind(samp, IarD)
    p.val <- ia.pval(index="Ia", samp2, IarD[1])
    p.val[2] <- ia.pval(index="rbarD", samp2, IarD[2])
    if(hist == TRUE){
      poppr.plot(samp, observed=IarD, pop=namelist$population,
                        file=namelist$File, pval=p.val, N=nrow(pop@tab))
    }
    result <- 1:4
    result[c(1,3)] <- IarD
    result[c(2,4)] <- p.val
    names(result) <- c("Ia","p.Ia","rbarD","p.rD")
    result <- c(result[1:2], sampstatIa, result[3:4], sampstatRd)
  } 
  return(final(Iout, result))
}

.new.Ia.Rd <- function (pop, missing = NULL) 
{
  vard.vector <- NULL
  numLoci <- length(pop)
  numIsolates <- length(pop[[1]]@ind.names)
  np <- choose(numIsolates, 2)
  if (np < 2) {
    return(as.numeric(c(NaN, NaN)))
  }
  V <- new.pair_diffs(pop, numLoci, np)
  varD <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector)
  sigVarj <- sum(vard.vector)
  rm(vard.vector)
  Ia <- (varD/sigVarj) - 1
  rbarD <- (varD - sigVarj)/(2 * sum(vardpair.vector))
  return(c(Ia, rbarD))
}

#==============================================================================#
# This creates a pairwise difference matrix via the C function pairdiffs in
# src/poppr_distance.c
# 
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # .Ia.Rd
#
#==============================================================================#
new.pair_diffs <- function(pop, numLoci, np)
{
  ploid <- ploidy(pop[[1]])
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab)*(ploid/2), numeric(np))
  return(temp.d.vector)
  d.vector <- colSums(temp.d.vector)
  d2.vector <- colSums(temp.d.vector^2)
  D.vector <- rowSums(temp.d.vector)
  return(list(d.vector = d.vector, d2.vector = d2.vector, D.vector = D.vector))
}

getinds <- function(x){
  # x is a vector indicating the size of loci in a data set (must be continuous).
  # The sum of x is equal to the number of observations.
  indices <- matrix(ncol = 2, nrow = length(x))
  to   <- cumsum(x)
  from <- c(1, to[-length(to)]+1)
  indices[] <- as.integer(c(from, to))
  return(indices)
}


brian.ia <- function(pop, sample=0, valuereturn = FALSE, method=1, 
                     quiet="minimal", missing="ignore",hist=TRUE){
  METHODS = c("multilocus", "permute alleles", "parametric bootstrap",
      "non-parametric bootstrap")
  namelist <- NULL
  namelist$population <- ifelse(length(levels(pop@pop)) > 1 | 
                                is.null(pop@pop), "Total", pop@pop.names)
  namelist$File <- as.character(pop@call[2])
  popx <- pop
  missing <- toupper(missing)
  type <- pop@type
  if(type=="PA"){
    .Ia.Rd <- .PA.Ia.Rd
  }
  else {
    popx <- seploc(popx)
  }

  # if there are less than three individuals in the population, the calculation
  # does not proceed. 
  if (nInd(pop) < 3){
    IarD <- as.numeric(c(NA,NA))
    names(IarD) <- c("Ia", "rbarD")
    if(sample==0){
      return(IarD)
    }
    else{
      IarD <- as.numeric(rep(NA,4))
      names(IarD) <- c("Ia","p.Ia","rbarD","p.rD")
      return(IarD)
    }
  }
  IarD <- .Ia.Rd(popx, missing)
  names(IarD) <- c("Ia", "rbarD")
  # no sampling, it will simply return two named numbers.
  if (sample==0){
    Iout <- IarD
    result <- NULL
  }
  # sampling will perform the iterations and then return a data frame indicating
  # the population, index, observed value, and p-value. It will also produce a 
  # histogram.
  else{
    Iout <- NULL 
    idx <- as.data.frame(list(Index=names(IarD)))
    samp <- .sampling(popx, sample, missing, quiet=quiet, type=type, method=method)
    samp2 <- rbind(samp, IarD)
    p.val <- ia.pval(index="Ia", samp2, IarD[1])
    p.val[2] <- ia.pval(index="rbarD", samp2, IarD[2])
    if(hist == TRUE){
      poppr.plot(samp, observed=IarD, pop=namelist$population,
                        file=namelist$File, pval=p.val, N=nrow(pop@tab))
    }
    result <- 1:4
    result[c(1,3)] <- IarD
    result[c(2,4)] <- p.val
    names(result) <- c("Ia","p.Ia","rbarD","p.rD")
    if (valuereturn == TRUE){
      return(list(index = final(Iout, result), samples = samp))
    }
  }  
  return(final(Iout, result))
}



jack.ia <- function(pop, iterations, type = type, divisor = 2){ 
  if(!is.list(pop)){
    inds <- nInd(pop)
  }
  else{
    inds <- nInd(pop[[1]])
  }
  half <- round(inds/divisor)
  np <- choose(half, 2)    
  V <- jack.Ia.Rd(pop)
	funrun <- function(V, inds = 10, half = 5, np = 10){
	
	  samp <- sample(inds, half)
    comb <- combn(inds, 2) # pairwise vector
    jackvec <- colSums(matrix(comb %in% samp, nrow = 2)) > 1

    V <- list(d.vector = colSums(V[jackvec, ]), 
              d2.vector = colSums(V[jackvec, ]^2), 
              D.vector = rowSums(V[jackvec, ])
              )
    IarD <- jack.calc(V, np)
	}
	sample.data <- vapply(1:iterations, function(x) funrun(V, inds = inds, half = half, np = np), c(pi, pi))
	rownames(sample.data) <- c("Ia", "rbarD")
	return(data.frame(t(sample.data)))
}

jack.Ia.Rd <- function (pop){

  vard.vector <- NULL
  if(!is.list(pop)){
    numLoci <- nLoc(pop)
    numIsolates <- nInd(pop)
  } else {
    numLoci <- length(pop)
    numIsolates <- length(pop[[1]]@ind.names)
  }
  np <- choose(numIsolates, 2)
  if (np < 2) {
    return(as.numeric(c(NaN, NaN)))
  }
  V <- jack.pair_diffs(pop, numLoci, np)
  return(V)
}

jack.calc <- function(V, np){
  varD <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector)
  sigVarj <- sum(vard.vector)
  rm(vard.vector)
  Ia <- (varD/sigVarj) - 1
  rbarD <- (varD - sigVarj)/(2 * sum(vardpair.vector))
  return(c(Ia, rbarD))
}

jack.pair_diffs <- function(pop, numLoci, np)
{
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  if(!is.list(pop)){
    ploid <- 2
    temp.d.vector <- vapply(seq(numLoci), 
                      function(x) as.vector(dist(pop@tab[,x])), 
                      temp.d.vector[,1])
    # checking for missing data and imputing the comparison to zero.
    if(any(is.na(temp.d.vector))){
      temp.d.vector[which(is.na(temp.d.vector))] <- 0
    }
  } 
  else {
    ploid <- ploidy(pop[[1]])
    temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab)*(ploid/2), 
                            temp.d.vector[, 1])
  }
  return(temp.d.vector)
}

# Function jackcomp is not fully thought out yet. Right now, it is unclear what 
# the application for this could be. The idea behind it is that if you randomly
# sampled half of the individuals from your population, how does that affect the
# index of association? This
#
#

jackcomp <- function(pop, sample = 999, quiet = TRUE, method = 1, divisor = 1.587302){
  inds <- nInd(pop)
  half <- round(inds/divisor)
  cat("Creating Null Distribution...\t")
  null <- brian.ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
    hist = FALSE, method = method)
  cat("Alternative Distribution...\t")
  if(pop@type == "codom"){
    popx <- seploc(pop)
  }
  else{
    popx <- pop
  }
  alt <- jack.ia(popx, sample, type = pop@type, divisor)
  # library(reshape)
  mnull <- melt(null$sample)
  malt <- melt(alt)
  mnull$Distribution <- "null"
  malt$Distribution <- "alternative"
  dat <- rbind(mnull, malt)
  cat("Creating Plots\n")
  obs.df <- data.frame(list(variable = names(null$samples), value = null$index[c(1,3)]))
  distplot <- jackbootplot(dat, obs.df) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals,", half, "Sampled for Alt.")
                     ) +
              theme(strip.text = element_text(size = rel(2)))
  #print(distplot)
  return(list(observed = null$index, null_samples = null$samples, alt_samples = alt, plot = distplot))
}

jackbootplot <- function(df, obs.df){
  Indexfac <- c("I[A]","bar(r)[d]")
  levels(df$variable) <- Indexfac
  obs.df$variable <- Indexfac
  dfIanull <- df[df$variable == Indexfac[1] & df$Distribution == "null", ]
  dfrbarDnull <- df[df$variable == Indexfac[2] & df$Distribution == "null", ]
  dfIaalt <- df[df$variable == Indexfac[1] & df$Distribution != "null", ]
  dfrbarDalt <- df[df$variable == Indexfac[2] & df$Distribution != "null", ]
  distplot <- ggplot(df, aes_string(x = "value", fill = "Distribution")) +
              geom_histogram(alpha = 0.5, position = "identity", 
                             data = dfIanull)+#,
                             # binwidth=diff(range(dfIanull$value))/30) +
              geom_histogram(alpha = 0.5, position = "identity", 
                             data = dfrbarDnull)+#,
                             # binwidth=diff(range(dfrbarDnull$value))/30) +
              geom_histogram(alpha = 0.5, position = "identity", 
                             data = dfIaalt)+#,
                             # binwidth=diff(range(dfIaalt$value))/30) +
              geom_histogram(alpha = 0.5, position = "identity", 
                             data = dfrbarDalt)+#,
                             # binwidth=diff(range(dfrbarDalt$value))/30) +

              geom_rug(alpha = 0.5, aes_string(color = "Distribution")) +
              facet_grid(" ~ variable", scales = "free_x", labeller = label_parsed) + theme_classic() +
              geom_vline(data = obs.df, aes_string(xintercept = "value", 
                                                   group = "variable"), 
                         linetype = "dashed") 
  return(distplot)      
}

#==============================================================================#
# This will perform a bootstrap analysis using the resulting distance matrix
# produced from pair_diffs(). Since the matrix is values of pairwise comparisons,
# the method of bootstrapping involves a few steps in order to resort the matrix.
#==============================================================================#
old.bootia <- function(pop, inds, iterations, comb){
  V <- jack.Ia.Rd(pop) # calculate matrix. rows = inds, cols = loci
  
  bootit <- function(V, inds, comb){
    #set.seed(9000)
    sam <- sample(inds, replace = TRUE) 
    #print(head(sam))
    
    
    # Indicating which columns contain ONLY values included in `sam`.
    uneak <- colSums(matrix(comb %in% sam, nrow = 2)) > 1
    samcomb <- comb[, uneak, drop = FALSE]
    
    comsam <- combn(sam, 2)
    compast <- vapply(1:ncol(comsam), function(x) paste(sort(comsam[, x]), collapse = ""), "a")
    sampast <- vapply(1:ncol(samcomb), function(x) paste(sort(samcomb[, x]), collapse = ""), "a")
    inds <- vapply(sampast, function(x) sum(compast %in% x), 1)
    samdump <- rep(as.numeric(names(inds)), inds)  
    
    # Now we have to subset the vector. Adding the duplicated comparisons will
    # not complete the entire data set, so the rest has to be filled in with 
    # zeroes, since that's the value one would obtain from comparisons to self.
    V[1:length(samdump), ] <- V[samdump, ]
    V[-(1:length(samdump)), ] <- 0
    
    # Calculate necessary vectors and send them to be calculated into Ia and rd.
    V <- list(d.vector = colSums(V), 
          d2.vector = colSums(V^2), 
          D.vector = rowSums(V)
          )
    return(jack.calc(V, ncol(comb)))
  }
  comb <- combn(inds, 2) # matrix of pairwise indexes.
  
  # Labelling the columns of this matrix provides tractability. When we subset
  # it, we will know where those subsetted columns came from. 
  colnames(comb) <- 1:ncol(comb)
	sample.data <- vapply(1:iterations, function(x) bootit(V, inds, comb), c(pi, pi))
	rownames(sample.data) <- c("Ia", "rbarD")
	return(data.frame(t(sample.data)))
}


bootia <- function(pop, inds, quiet = TRUE, method = 1){
  pop <- pop[sample(inds, replace = TRUE), ]
  pop@ind.names <- paste("ind", 1:inds)
  return(.ia(pop, quiet = quiet))
} 

bootcomp <- function(pop, sample = 999, quiet = TRUE, method = 1){
  inds <- nInd(pop)
  cat("Creating Null Distribution...\t")
  null <- brian.ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
    hist = FALSE, method = method)
  cat("Alternative Distribution...\t")
  alt <- vapply(1:sample, function(x) bootia(pop, inds, quiet, method), c(pi, pi))
  #if(pop@type == "codom"){
  #  popx <- seploc(pop)
  #}
  #else{
  #  popx <- pop
  #}
  #alt <- bootia(popx, inds, sample)
  alt <- data.frame(t(alt))
  # library(reshape)
  mnull <- melt(null$sample)
  malt <- melt(alt)
  mnull$Distribution <- "null"
  malt$Distribution <- "alternative"
  dat <- rbind(mnull, malt)
  cat("Creating Plots\n")
  obs.df <- data.frame(list(variable = names(null$samples), value = null$index[c(1,3)]))
  distplot <-distplot <- jackbootplot(dat, obs.df) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals, Method: Bootstrap")
                     ) +
              theme(strip.text = element_text(size = rel(2)))
  #print(distplot)
  return(list(observed = null$index, null_samples = null$samples, alt_samples = alt, plot = distplot))
}

#' @importFrom reshape2 melt

jackbootcomp <- function(pop, sample = 999, quiet = TRUE, method = 1, divisor = 1.587302){
  inds <- nInd(pop)
  cat("Creating Null Distribution...\t")
  null <- brian.ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
    hist = FALSE, method = method)
  if(pop@type == "codom"){
    popx <- seploc(pop)
  }
  else{
    popx <- pop
  }
  cat("Alternative Distribution (Bootstrap)...\t")
  altboot <- vapply(1:sample, function(x) bootia(pop, inds, quiet, method), c(pi, pi))
  altboot <- data.frame(t(altboot))
  #altboot <- bootia(popx, inds, sample)
  cat("Alternative Distribution (Jack Knife)...\t")
  half <- round(inds/divisor)
  altjack <- jack.ia(popx, sample, type = pop@type, divisor)  
  # library(reshape)
  mnull <- melt(null$sample)
  bmalt <- melt(altboot)
  jmalt <- melt(altjack)
  mnull$Distribution <- "null"
  bmalt$Distribution <- "bootstrap"
  jmalt$Distribution <- "jack knife"
  dat <- rbind(mnull, bmalt, jmalt)
  cat("Creating Plots\n")
  obs.df <- data.frame(list(variable = names(null$samples), value = null$index[c(1,3)]))
  distplot <-distplot <- jackbootplot(dat, obs.df) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals,", half, "Sampled for Jack knife")
                     ) +
              theme(strip.text = element_text(size = rel(2)))
  #print(distplot)
  return(list(observed = null$index, null_samples = null$samples, 
              alt_samples_boot = altboot, alt_samples_jack = altjack, 
              plot = distplot)
        )
} 






################################################################################
#################### Zhian's Functions above ###################################
################################################################################
################################################################################
################################################################################
#################### DO NOT CROSS THIS LINE!!!! ################################
################################################################################
################################################################################
################################################################################
#################### Javier's Functions below ##################################
################################################################################

genoid.bruvo.boot <- function(pop, replen = 1, add = TRUE, loss = TRUE, 
                              sample = 100,  tree = "upgma", showtree = TRUE, 
                              cutoff = NULL,  quiet = FALSE, ...){
  # This attempts to make sure the data is true microsatellite data. It will
  # reject snp and aflp data. 
  if(pop@type != "codom" | all(is.na(unlist(lapply(pop@all.names, as.numeric))))){
    stop("\nThis dataset does not appear to be microsatellite data. Bruvo's Distance can only be applied for true microsatellites.")
  }
  ploid <- ploidy(pop)
  # Bruvo's distance depends on the knowledge of the repeat length. If the user
  # does not provide the repeat length, it can be estimated by the smallest
  # repeat difference greater than 1. This is not a preferred method. 
  if (length(replen) != length(pop@loc.names)){
    replen <- vapply(pop@all.names, function(x) guesslengths(as.numeric(x)), 1)
    #    replen <- rep(replen[1], numLoci)
    warning("\n\nRepeat length vector for loci is not equal to the number of loci represented.\nEstimating repeat lengths from data:\n", immediate.=TRUE)
    cat(replen,"\n\n")
  }
  # This controlls for the user correcting missing data using "mean". 
  if(any(!pop@tab,10 %in% c(0,((1:ploid)/ploid), 1, NA))){
    pop@tab[!pop@tab %in% c(0,((1:ploid)/ploid), 1, NA)] <- NA
  }
  # Converting the genind object into a matrix with each allele separated by "/"
  bar <- as.matrix(genind2df(pop, sep="/", usepop = FALSE))
  # The bruvo algorithm will ignore missing data, coded as 0.
  bar[bar %in% c("", NA)] <- paste(rep(0, ploid), collapse="/")
  # stopifnot(require(phangorn))
  # Steps: Create initial tree and then use boot.phylo to perform bootstrap
  # analysis, and then place the support labels on the tree.
  if(tree == "upgma"){
    newfunk <- match.fun(upgma)
  }
  else if(tree == "nj"){
    newfunk <- match.fun(nj)
  }
  tre <- newfunk(phylo.bruvo.dist(bar, replen=replen, ploid=ploid, add = add, loss = loss))
  if (any (tre$edge.length < 0)){
    tre$edge.length[tre$edge.length < 0] <- 0
  }
  if(quiet == FALSE){
    cat("\nBootstrapping... (note: calculation of node labels can take a while even after the progress bar is full)\n\n")
  }
  bp <- boot.phylo(tre, bar, FUN = function (x) newfunk(phylo.bruvo.dist(x, replen = replen, ploid = ploid, add = add, loss = loss)), B = sample, quiet = quiet, ...)
  tre$node.labels <- round(((bp / sample)*100))
#  if (!is.null(cutoff)){
#    if (cutoff < 1 | cutoff > 100){
#      cat("Cutoff value must be between 0 and 100.\n")
#      cutoff<- as.numeric(readline(prompt = "Choose a new cutoff value between 0 and 100:\n"))
#    }
    tre$node.labels[tre$node.labels < cutoff] <- NA
#  }
  tre$tip.label <- pop@ind.names
#  if(showtree == TRUE){
#    plot(tre, show.node.label=TRUE)
#  }
#  if(tree=="upgma"){
#    axisPhylo(3)
#  }
  return(tre)
}

javier<-function(x){
  cat ("http://www.youtube.com/watch?v=1-ctsxVXvO0")
}

#######################################################
### Functions for Sequence Data ------------------- ###
### Tajima's D, Nuc. Diversity and Ka/Ks estimation ###
#######################################################
#
#
# using the packages pegas, ape and sequin, Im creating some wrappers to be able to do multiple calculations of basic pop. gen. statistical tests such as Tajimas D, Nucleotide diversity (pi) and dN/dS selection (as Ka/Ks)
#
#
####################################################################################################################################
# ## Tajimas D calculation
# Function td: Tajima-D wrapper to get a data frame out of the data with the Tajimas D value, the p-values for normality and poisson, the files and if the p-value (for normality) is significant or not.
# Usage: 
# td(fasta file location)
####################################################################################################################################
#' td
#' 
#' Tajima-D wrapper of \code{pegas} function \code{\link{tajima.test}} to get a data frame out of the data with the Tajima D value, the p-values for normality and poisson, the files and if the p-value (for normality) is significant or not.
#' 
#' @param x an alignment \code{fasta} file (preferably full path, fasta|fas extension).
#' 
#' @return a table with the Tajima's D calculation results: Tajimas D value, the p-values for normality and poisson, the files and if the p-value (for normality) is significant or not.
#' 
#' @seealso \code{\link{tajima.test}}
#' 
#' @note The input FASTA file MUST BE AN ALIGNMENT! Be sure to use the full path. The FASTA files are read by \code{\link{ape}} function \code{read.dna} and then processed and ran trough \code{pegas} function \code{\link{tajima.test}}.
#'
#' @export
#' 
#' @author Javier F. Tabima
#' 
#' @references
#' Tajima, F. (Nov 1989). "Statistical method for testing the neutral mutation hypothesis by DNA polymorphism.". 
#' Genetics 123 (3): 585-95. PMC 1203831. PMID 2513255.
#' 
#' @examples
#'
#' # Analyze one of the included FASTA alignment files
#' # td("data_1.fasta",format="fasta")
#' 
#' @importFrom ape read.dna
#' @importFrom pegas tajima.test
#' 
td<- function (x){
  seq <- read.dna(x,format = "fasta")
  tt <- tajima.test(seq)
  tt <- unlist(tt)
  return ((tt))
}
############################################################################
## Nucleotide Diversity
#Function n.diversity: Nucleotide diversity calculation for sequence data.
# Usage:
# n.diversity(fasta file location)
############################################################################
#' n.diversity
#' 
#' Nucleotide diversity calculation for sequence data. Wrapper of \code{pegas} function \code{nuc.div} to get a data frame out of the data with the nucleotide diversity estimate.
#' 
#' @param x an alignment \code{fasta} file (preferably full path, fasta|fas extension).
#' 
#' @return a table with the Nucleotide Diversity calculation results.
#' 
#' @seealso \code{\link{nuc.div}}
#' 
#' @note The input FASTA file MUST BE AN ALIGNMENT! Be sure to use the full path. The FASTA files are read by \code{\link{ape}} function \code{read.dna} and then processed and ran trough \code{pegas} function \code{\link{nuc.div}}.
#'
#' @export
#' 
#' @author Javier F. Tabima
#' 
#' @references
#' Nei, M. (1987) Molecular evolutionary genetics. New York: Columbia University Press.
#' Nei, M and Li, W. 1979. "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". PNAS 76 
#' 
#' @examples
#'
#' # Analyze one of the included FASTA alignment files
#' # n.diversity("data_1.fasta",format="fasta")
#'  
#' @importFrom ape read.dna
#' @importFrom pegas nuc.div
#' 
n.diversity <- function (x){
  seq <- read.dna(x, format="fasta")
  nuc.div <- nuc.div(seq)
  unlist(nuc.div)
  return(nuc.div)
}
############################################################################
## dN/dS Estimation
# Function dnds: Calculation of dN/dS
############################################################################

#' dn.ds
#' 
#' dN/dS (as Ka/Ks estimation) estimation of selection. Wrapper of \code{seqinr} function \code{\link{kaks}}
#' 
#' @param x an alignment \code{fasta} file (preferably full path, fasta|fas extension).
#' 
#' @return a table with the dN/dS (Ka/Ks) value and the estimate for selection (Positive, Negative or Neutral)
#' 
#' @seealso \code{\link{kaks}}
#' 
#' @note The input FASTA file MUST BE AN ALIGNMENT! Be sure to use the full path. The FASTA files are read by \code{\link{seqinr}} function \code{read.alignment}.
#'
#' @export
#' 
#' @author Javier F. Tabima
#' 
#' @references
#' Li WH, Wu CI, Luo CC. 1985. A new method for estimating synonymous and nonsynonymous rates of nucleotide substitution considering the relative likelihood of nucleotide and codon changes. Mol Biol Evol 2(2):150-174. 
#' Zhang, Z., Li, J., Zhao, X., Wang, J., Wong, G.K. and Yu, J. (2006) KaKs_Calculator: calculating Ka and Ks through model selection and model averaging, Genomics Proteomics Bioinformatics, 4(4): 259-263.
#' @examples
#'
#' # Analyze one of the included FASTA alignment files
#' # dn.ds("data_1.fasta",format="fasta")
#' @importFrom seqinr read.alignment kaks
#' 
dn.ds <- function (x){
  seq <- read.alignment(x, format="fasta") 
  if (!seq$nb == 2){
    stop("This alignment has more than 2 sequences. dN/dS compares between orthologous genes in a pair of species")
  } else {
    tab <- suppressWarnings(kaks(seq))
    tab <- unlist(tab)
    if (is.na(tab)){
      tab <- as.numeric(c(ka = NA, ks = NA, vka = NA, vks = NA))
    }
    return(tab)
  }
}

############################################################################
############################################################################
####
#### Wrappers for multiple sequences 
####
############################################################################
############################################################################

############################################################################
# Function multi.td: Tajima-D wrapper for multiple sets
# 
# Usage: 
# multi.td(list of fasta files)
# 
# To get the list of FASTA files use list.files()
############################################################################

#' multi.td
#' 
#' Tajima-D wrapper for multiple sets of FASTA files and the \code{\link{td}} function
#' 
#' @param x a \code{list} of multiple alignment \code{fasta} files (preferably full path, fasta|fas extension).
#' @param cutoff a \code{integer} of the cutoff value for the p-statistic.
#' @param quiet a \code{boolean} for verbosity.
#' 
#' @return a table with the Tajima's D calculation for all the files in the list.
#' 
#' @seealso \code{\link{td}}
#' 
#' @note To create the list of files use the \code{\link{list.files}} function. (See examples for help)
#'
#' @export
#' 
#' @author Javier F. Tabima
#' 
#' @examples
#'
#' # Create the list of FASTA files
#' #list.pairs <- list.files(path="fasta_file_folder/", pattern="*.fasta$", full.names=TRUE)
#' #
#' # Run the function
#' # multi.td(list.pairs)
#' 
multi.td <- function(x, cutoff = 0.05, quiet = FALSE){
  if (quiet == FALSE){
    cat("Tajimas D Calculation\n")
  }
  tt.df <- t(vapply(x, td,numeric(3))) 
  rownames(tt.df) <- sub("^([^.]*).*", "\\1", basename(x)) 
  tt.df <- as.data.frame(tt.df)
  pna <- is.na(tt.df$Pval.normal)  
  pcutoff <- tt.df$Pval.normal > cutoff
  tt.df$stat[pcutoff] <- "Neutral"
  tt.df$stat[!pcutoff&tt.df$D > 1] <- "Positive"
  tt.df$stat[!pcutoff&tt.df$D < 1] <- "Negative"
  tt.df$stat[pna] <- "Non_segregating"
  return(tt.df)  
}
############################################################################
#Function multi.nd: Nucleotide diversity wrapper for multiple sets
# 
# Usage: 
# multi.nd(list of fasta files)
# 
# To get the list of FASTA files use list.files() 
############################################################################
#' multi.nd
#' 
#' Nucleotide diversity wrapper for multiple sets of FASTA files and the \code{\link{n.diversity}} function
#' 
#' @param x a \code{list} of multiple alignment \code{fasta} files (preferably full path, fasta|fas extension).
#' @param quiet a \code{boolean} for verbosity.
#'
#' 
#' @return a table with the nucleotide diversity calculation for all the files in the list.
#' 
#' @seealso \code{\link{td}}
#' 
#' @note To create the list of files use the \code{\link{list.files}} function. (See examples for help)
#'
#' @export
#' 
#' @author Javier F. Tabima
#' 
#' @examples
#'
#' # Create the list of FASTA files
#' #list.pairs <- list.files(path="fasta_file_folder/", pattern="*.fasta$", full.names=TRUE)
#' #
#' # Run the function
#' # multi.nd(list.pairs)
#' 
multi.nd <- function(x, quiet = FALSE){
  if (quiet == FALSE){
    cat("Nucleotide Diversity\n")
  }
  nd.df <- (vapply(x,n.diversity,numeric(1))) # Add spaces after commas
  nd.df <- as.data.frame(nd.df)
  rownames(nd.df) <- sub("^([^.]*).*", "\\1", basename(x)) 
  colnames(nd.df) <- c("Nuc.Div")
  return(nd.df)  
}
##########################################################################

# Function multi.dnds: dN/dS wrapper for multiple pairs
# 
# Usage: 
# multi.dnds(list of fasta files)
# 
# To get the list of FASTA files use list.files()
############################################################################
#' multi.dnds
#' 
#' dN/dS (as Ka/Ks) wrapper for multiple sets of FASTA files and the \code{\link{dn.ds}} function
#' 
#' @param x a \code{list} of multiple alignment \code{fasta} files (preferably full path, fasta|fas extension).
#' @param quiet a \code{boolean} for verbosity.
#' 
#' @return a table with the dN/dS calculation for all the files in the list.
#' 
#' @seealso \code{\link{dn.ds}}
#' 
#' @note To create the list of files use the \code{\link{list.files}} function. (See examples for help)
#'
#' @export
#' 
#' @author Javier F. Tabima
#' 
#' @examples
#'
#' # Create the list of FASTA files
#' #list.pairs <- list.files(path="fasta_file_folder/", pattern="*.fasta$", full.names=TRUE)
#' #
#' # Run the function
#' # multi.dnds(list.pairs)

multi.dnds <- function(x, quiet=FALSE){
  if (quiet == FALSE){  
    cat("dN/dS calculation\n")
  }
  dnds.df <- (vapply(x, dn.ds,numeric(4)))
  dnds.df <-as.data.frame(t(dnds.df))
  rownames(dnds.df) <- sub("^([^.]*).*", "\\1", basename(x)) 
  dnds.df$dnds <- dnds.df$ka/dnds.df$ks
  dnds.df$stat[dnds.df$dnds>1] <- "Positive"
  dnds.df$stat[dnds.df$dnds<1] <- "Negative"
  return(dnds.df)  
}