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


informloci <- function(pop, min_ind = 2, quiet = FALSE){
  if(!is.genind(pop)){
    stop("This function only works on genind objects.")
  }
  if(pop@type == "PA"){
    
    locivals <- apply(pop@tab, 2, sum) %in% min_ind:(nInd(pop) - min_ind)
    if(all(locivals == TRUE)){
      cat("No sites found with fewer than", min_ind, 
          "different individuals.\n", fill = 80)
    }
    else{
      cat(sum(!locivals), "uninformative", 
          ifelse(sum(!locivals) > 1, "loci", "locus"), "found:", 
          pop@loc.names[!locivals],"\n", fill = 80)
    }
    return(pop[, locivals])
  }
  else{
    locivals <- apply(as.loci(pop)[-1], 2, testable, min_ind)
    if(all(locivals == TRUE)){
      cat("No sites found with fewer than", min_ind, 
          "different individuals.\n", fill = 80)
    }
    else{
      cat(sum(!locivals), "uninformative", 
          ifelse(sum(!locivals) > 1, "loci", "locus"), "found:", 
          pop@loc.names[!locivals],"\n", fill = 80)
    }
    return(pop[, loc = names(pop@loc.names[locivals])])
  }
}

testable <- function(loc, min_ind){
  tab <- table(loc)
  if(length(tab) > 2){
    return(TRUE)
  }
  else{
    if(length(tab) < 2){
      return(FALSE)
    }
    else{
      return(ifelse(any(tab < min_ind), FALSE, TRUE))        
    }
  }
}





testing_funk <- function(){
  cat("This test worked...maybe.\n")
}



.new.sampling <- function(pop,iterations,quiet="noisy",missing="ignore",type=type, method=1){ 
  METHODS = c("multilocus", "permute alleles", "parametric bootstrap",
      "non-parametric bootstrap")
  if(!is.list(pop)){
    if(type=="PA"){
      .Ia.Rd <- .PA.Ia.Rd
    }
  }
	sample.data <- NULL
	for (c in 1:iterations){
    IarD <- .Ia.Rd(total.shuffler(pop, method=method), missing=missing)
		sample.data <- rbind(sample.data, as.data.frame(list( 
                                      Ia=IarD[1],
                                      rbarD=IarD[2]
                                      )))
    if (quiet != TRUE){
      if(quiet == "noisy"){
        cat("Sample: ",c,"\n")
        cat("Index of Association: ", IarD[1],"\n")
        cat("Standardized Index of Association (rbarD): ", IarD[2],"\n")
      }
      else{
        if(c%%50 != 0){
          cat(".")
          if (c == iterations){
            cat("\n")
          }
#          else if( c == 666 ){
#            cat("\\m/")
#          }
        }
        else{
          cat(".")
          cat("\n")
        }      
      }
    }
	}
	return(sample.data)
}

new.poppr <- function(pop,total=TRUE,sublist=c("ALL"),blacklist=c(NULL), sample=0,
  method=1,missing="ignore", quiet="minimal",clonecorrect=FALSE,hier=c(1),dfname="hier",
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
    # Shannon-Weiner vegan:::diversity index.
    H <- vegan:::diversity(pop.mat)
    # E_1, Pielou's evenness.
    # J <- H / log(rowSums(pop.mat > 0))
    # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
    G <- vegan:::diversity(pop.mat, "inv")
    Hexp <- (N.vec/(N.vec-1))*vegan:::diversity(pop.mat, "simp")
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
    # Shannon-Weiner vegan:::diversity index.
    H <- vegan:::diversity(pop.mat)
    # E_1, Pielou's evenness.
    # J <- H / log(rowSums(pop.mat > 0))
    # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
    G <- vegan:::diversity(pop.mat, "inv")
    Hexp <- (N.vec/(N.vec-1))*vegan:::diversity(pop.mat, "simp")
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

.new.ia <- function(pop,sample=0,method=1,quiet="minimal",namelist=NULL,missing="ignore",
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
  .quiet(quiet=quiet, IarD=IarD, pop=namelist$population)
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
      if(require(ggplot2)){
        poppr.plot(samp, observed=IarD, pop=namelist$population,
                          file=namelist$File, pval=p.val, N=nrow(pop@tab))
      }
      else{      
        permut.histogram(samp, IarD, p.val[1], pop=namelist$population, 
                        file=namelist$File)
      }
    }
    result <- 1:4
    result[c(1,3)] <- IarD
    result[c(2,4)] <- p.val
    names(result) <- c("Ia","p.Ia","rbarD","p.rD")
    result <- c(result[1:2], sampstatIa, result[3:4], sampstatRd)
  } 
  return(final(Iout, result))
}



getinds <- function(x){
  # x is a vector indicating the size of loci in a data set (must be continuous).
  # The sum of x is equal to the number of observations.
  to <- cumsum(x)
  from <- c(1, to[-length(to)]+1)
  indices <- rbind(from, to)
  rownames(indices) <- NULL
  return(indices)
}

total.shuffler <- function(pop, method){
  return(lapply(pop, pop.sampler, method))
}
  


pop.sampler <- function(pop, method=1){
  pop_names <- names(table(pop@pop))
  total_inds <- 1:length(pop@pop)
  repop <- function(x, pop, pop_names, total_inds, method){
    pop_vec <- total_inds[pop@pop %in% pop_names[x]]
    if(all(is.na(pop@tab[pop_vec, ]))){
      pop@tab[pop_vec, ] <<- pop@tab[pop_vec, ]
      return(1)
    }
    cols <- which(colSums(pop[pop_vec, ]@tab, na.rm=TRUE) > 0)
    newpop <- pop
    if(method != 2){
      newpop@tab[pop_vec, -cols] <- 0
    }
    newpop@tab[pop_vec, cols] <- .locus.shuffler(popsub(pop, x), method)@tab
    if(method == 1){
      replacements <- is.na(rowSums(newpop@tab[pop_vec, cols]))
      if(any(replacements == TRUE)){
        #cat(pop_vec,"\t",replacements,"\n")
        newpop@tab[pop_vec[replacements], -cols] <- NA
      }
    }
    pop@tab[pop_vec, ] <<- newpop@tab[pop_vec, ]
    return(0)
  }
  invisible(lapply(1:length(pop_names), repop, pop, pop_names, total_inds, method))
#   nas <- function(x, pop, pop_names, total_inds, method){
#     pop_vec <- total_inds[pop@pop %in% pop_names[x]]
#     if(any(is.na(pop@tab[pop_vec, ]))){
#       newpop <- pop
#       replacena <- function(locus, newpop, pop_vec){
#         locus <- newpop@loc.fac %in% locus
#         truth <- is.na(rowSums(newpop@tab[pop_vec, locus]))
#         if(any(truth == TRUE)){
#           #cat(pop_vec,"\n",as.numeric(locus),"\n", truth, "\n")
#           newpop@tab[pop_vec[truth], locus] <<- NA
#         }
#       }
#       invisible(lapply(names(pop@loc.names), replacena, newpop, pop_vec))
#       pop@tab[pop_vec, ] <<- newpop@tab[pop_vec, ]
#     }
#   }
#   if(method == 1){
#     invisible(lapply(1:length(pop_names), nas, pop, pop_names, total_inds, method))
#   }
  return(pop)
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

genoid.bruvo.boot <- function(pop, replen=c(2), sample = 100, tree = "upgma", 
                       showtree=TRUE, cutoff=NULL, quiet=FALSE, ...) {
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
  if(any(!round(pop@tab,10) %in% c(0,((1:ploid)/ploid),1, NA))){
    pop@tab[!round(pop@tab,10) %in% c(0,((1:ploid)/ploid),1, NA)] <- NA
  }
  # Converting the genind object into a matrix with each allele separated by "/"
  bar <- as.matrix(genind2df(pop, sep="/", usepop=FALSE))
  # The bruvo algorithm will ignore missing data, coded as 0.
  bar[bar %in% c("", NA)] <- paste(rep(0, ploid), collapse="/")
  stopifnot(require(phangorn))
  # Steps: Create initial tree and then use boot.phylo to perform bootstrap
  # analysis, and then place the support labels on the tree.
  if(tree == "upgma"){
    newfunk <- match.fun(upgma)
  }
  else if(tree == "nj"){
    newfunk <- match.fun(nj)
  }
  tre <- newfunk(phylo.bruvo.dist(bar, replen=replen, ploid=ploid))
  if (any (tre$edge.length < 0)){
      tre$edge.length[tre$edge.length < 0] <- 0
    }
  if(quiet == FALSE){
    cat("\nBootstrapping... (note: calculation of node labels can take a while even after the progress bar is full)\n\n")
  }
  bp <- boot.phylo(tre, bar, FUN = function (x) newfunk(phylo.bruvo.dist(x, replen = replen, ploid = ploid)), B = sample, quiet=quiet, ...)
  tre$node.labels <- round(((bp / sample)*100))
  if (!is.null(cutoff)){
    if (cutoff < 1 | cutoff > 100){
      cat("Cutoff value must be between 0 and 100.\n")
      cutoff<- as.numeric(readline(prompt = "Choose a new cutoff value between 0 and 100:\n"))
    }
    tre$node.labels[tre$node.labels < cutoff]<-NA
  }
  tre$tip.label <- pop@ind.names
  if(showtree == TRUE){
    plot(tre, show.node.label=TRUE)
  }
  if(tree=="upgma"){
    axisPhylo(3)
  }
  return(tre)
}









javier<-function(x){
  cat ("http://www.youtube.com/watch?v=1-ctsxVXvO0")
}
