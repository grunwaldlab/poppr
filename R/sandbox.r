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



myTheme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

read_allele_columns <- function(x, ploidy = 2){
  clm <- ncol(x)
  if (length(ploidy) == 1){
    if (clm %% ploidy == 0){
      ploidy <- rep(ploidy, clm/ploidy)
    } else {
      stop("ploidy specified does not match number of columns.")
    }
  }
  loci  <- first_allele_col(ploidy)
  x2    <- x[, loci]
  lapply(1:length(loci), function(a) x2[, a] <<-
           pop_combiner(x, hier = loci[a]:(loci[a]+ploidy[a]-1), sep = "/"))
  
  x2 <- data.frame(lapply(x2, add_zeroes,  ploidy = max(ploidy)))
  x2[apply(x2, 2, function(na) grepl("NA", na))] <- NA
  return(x2)
}

first_allele_col <- function(ploidies) cumsum(ploidies) - ploidies + 1

# Function that will add extra zeroes onto genotypes that are deficient in the 
# right number of alleles.
add_zeroes <- function(x, ploidy = 3){
  extras <- ploidy - vapply(strsplit(x, "/"), length, 69)
  vapply(1:length(extras), zero_adder, character(1), extras, x)
}

zero_adder <- function(index, extras, df){
  if(extras[index] > 0){
    return(paste(c(rep(0, extras[index]), df[index]), collapse = "/"))
  } else {
    return(df[index])
  }
}




# Function that was to be used to calculate standard errors for bootstrapping
# diversity statistics. Probably will not use this one.
boot_se_table <- function(tab, n = 1000, ci = 95, total = TRUE, rarefy = FALSE, 
                          rare.val = NULL, ...){
  if (!is.matrix(tab) & is.genind(tab)){
    tab <- mlg.table(tab, total = total, plot = FALSE)
  }
  rareval <- NULL
  
  if (rarefy){
    rareval <- ifelse(is.null(rare.val), min(rowSums(tab)), rare.val)
  }
  res <- diversity_boot(tab, n, n.rare = rareval, ...)
  if (rarefy){
    orig <- get_boot_se(res, "mean")
  } else {
    orig <- get_boot_stats(res)
  }
  out  <- matrix(nrow = nrow(orig), ncol = ncol(orig)*2, 
                 dimnames = list(rownames(orig), NULL))
  out[, !evens(1:ncol(out))]  <- orig
  se_out <- get_boot_se(res, "sd")
  out[, evens(1:ncol(out))] <- se_out
  colnames(out) <- intersp(colnames(orig), colnames(se_out))
  return(out)
}

# Attempt at rarefcation for the index of association.
rare_ia <- function(x, n = 1000, rare = 10, obs = FALSE){
  if (is.genind(x) || is.clone(x) || is(x, "genlight")){
    xloc <- seploc(x)
    est <- bootjack(xloc, n, rare, progbar = NULL)
    se  <- vapply(est, sd, numeric(1), na.rm = TRUE)
    if (obs){
      res <- matrix(nrow = 4, ncol = 1,
                    dimnames = list(c("Ia", "Ia.se", "rbarD", "rbarD.se"), NULL))
      est <- vapply(est, mean, numeric(1), na.rm = TRUE)
      res[c("Ia", "rbarD"), ] <- est
    } else {
      res <- matrix(nrow = 2, ncol = 1,
                    dimnames = list(c("Ia.se", "rbarD.se"), NULL))
    } 
    res[c("Ia.se", "rbarD.se"), ] <- se
  } else {
    if (obs){
      res <- matrix(nrow = 4, ncol = length(x),
                    dimnames = list(c("Ia", "Ia.se", "rbarD", "rbarD.se"), NULL))
    } else {
      res <- matrix(nrow = 2, ncol = length(x),
                    dimnames = list(c("Ia.se", "rbarD.se"), NULL))
    }
    for (i in seq(length(x))){
      if (nInd(x[[i]]) > rare){
        iloc <- seploc(x[[i]])
        est  <- bootjack(iloc, n , rare, progbar = NULL)
        se   <- vapply(est, sd, numeric(1), na.rm = TRUE)
        if (obs) est <- vapply(est, mean, numeric(1), na.rm = TRUE)
      } else {
        if (obs) est <- ia(x[[i]])
        se  <- c(0, 0)
      }
      if (obs) res[c("Ia", "rbarD"), i] <- est
      res[c("Ia.se", "rbarD.se"), i] <- se
    }
  }
  return(res)
}


old_pair_ia <- function(pop){
  
  if(pop@type == "codom"){
    pop_loci <- seploc(pop)
    loci_pairs <- combn(locNames(pop), 2)
    pair_ia_vector <- apply(loci_pairs, 2, function(x) .Ia.Rd(pop_loci[x]))
    colnames(pair_ia_vector) <- apply(loci_pairs, 2, paste, collapse = ":")
  } else {
    loci_pairs <- combn(1:nLoc(pop), 2)
    pair_ia_vector <- apply(loci_pairs, 2, function(x) .PA.Ia.Rd(pop[, x], missing = "ignore"))
    colnames(pair_ia_vector) <- apply(combn(locNames(pop), 2), 2, paste, collapse = ":")
  }
  rownames(pair_ia_vector) <- c("Ia", "rbarD")    
  return(pair_ia_vector)
}




ia_pair_loc <- function(pair, V, np, progbar, iterations){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, as.numeric(pair[3])/iterations)
  }
  newV <- V[, pair[-3]]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = rowSums(newV)
  )
  return(jack.calc(V, np))
}


poppr.pair.ia <- function(pop){
  if(is.null(pop(pop))){
    return(pair.ia(pop))
  }
  pops       <- seppop(pop, drop = FALSE)
  loci_pairs <- choose(nLoc(pop), 2)
  res_mat    <- matrix(0.5, 2, loci_pairs)
  pops_array <- vapply(pops, pair.ia, res_mat)
  return(pops_array)
}

testing_funk <- function(){
  cat("This test worked...maybe.\n")
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

#==============================================================================#
# bootjack is a function that will calculate values of the index of association
# after sampling with or without replacement. The purpose of this function is
# to help determine distributions of I_A under non-random mating by creating 
# exact copies of samples by sampling with replacement and also to derive a
# distribution of the data by sampling a subset of the data (63% by default)
# without replacement. 
#
# Since the data itself is not being changed, we can use the distances observed.
# These distances are calculated per-locus and presented in a matrix (V) with the
# number of columns equal to the number of loci and the number of rows equal to
# choose(N, 2) where N is the number of samples in the data set. Sampling has
# been optimized with this strategy by first creating a lower triangle distance
# matrix containing indices from 1 to choose(N, 2). This is then converted to a
# square matrix (mat) and subset after sampling with or without replacement. The
# remaining indices in the distance matrix is used to subset the distance-per-
# locus matrix (V). Calculations are then performed on this matrix. It must be
# noted that for sampling with replacement, duplicated indices must be
# supplemented with rows of zeroes to indicate no distance.  
#==============================================================================#
bootjack <- function(gid, iterations = 999, half = NULL, progbar = NULL){
  if (!is.list(gid)){
    N <- nInd(gid)
    numLoci <- nLoc(gid)
  } else {
    N <- nInd(gid[[1]])
    numLoci <- length(gid)
  }
  # Step 1: make a distance matrix defining the indices for the pairwise 
  # distance matrix of loci.
  np  <- choose(N, 2)
  dis <- 1:np
  dis <- make_attributes(dis, N, 1:N, "dist", call("dist"))
  mat <- as.matrix(dis)
  # Step 2: calculate the pairwise distances for each locus. 
  V   <- jack.pair_diffs(gid, numLoci, np)
  if (!is.null(half)){
    np <- choose(half, 2)
    sample.data <- vapply(1:iterations, jackrun, numeric(2), V, mat, N, half, 
                          np, progbar, iterations)
  } else {
    sample.data <- vapply(1:iterations, bootrun, numeric(2), V, mat, N, np, 
                          progbar, iterations)
  }
  rownames(sample.data) <- c("Ia", "rbarD")
  return(data.frame(t(sample.data)))
}

#------------------------------------------------------------------------------#
# Variables for bootrun and jackrun
#
# count      = the current iteration
# V          = the matrix of pairwise distances per locus
# mat        = square matrix containing indices for V
# N          = number of observations in data set
# half       = the number of observations to be sampled without replacement
# np         = choose(N, 2)
# progbar    = a progress bar object
# iterations = the number of total samples
#------------------------------------------------------------------------------#
bootrun <- function(count, V, mat, N, np, progbar, iterations){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, count/iterations)
  }
  # Step 3: sample individuals with replacement and subset distance matrix
  # from step 1.
  inds    <- sample(N, replace = TRUE)
  newInds <- as.vector(as.dist(mat[inds, inds])) 
  
  # Step 4: Find the number of zeroes in from the distance matrix. This will be 
  # the number of duplicated genotypes.
  zeroes <- length(inds == 0)

  # Step 5: subset the pairwise locus matrix.
  newV <- V[newInds, ]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = c(rowSums(newV), rep(0, zeroes))
              )
  # Step 6: Calculate the index of association.
  return(jack.calc(V, np))
}

jackrun <- function(count, V, mat, N, half, np, progbar, iterations){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, count/iterations)
  }

  inds    <- sample(N, half)
  newInds <- as.vector(as.dist(mat[inds, inds])) 

  newV <- V[newInds, ]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = rowSums(newV)
              )
  return(jack.calc(V, np))
}

jack.calc <- function(V, np){
  varD <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector, PACKAGE = "poppr")
  sigVarj <- sum(vard.vector)
  rm(vard.vector)
  Ia <- (varD/sigVarj) - 1
  rbarD <- (varD - sigVarj)/(2 * sum(vardpair.vector))
  return(c(Ia, rbarD))
}

jack.pair_diffs <- function(pop, numLoci, np){
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  if (!is.list(pop)){
    ploid <- 2
    temp.d.vector <- vapply(seq(numLoci), 
                      function(x) as.vector(dist(pop@tab[,x])), 
                      temp.d.vector[,1])
    # checking for missing data and imputing the comparison to zero.
    if(any(is.na(temp.d.vector))){
      temp.d.vector[which(is.na(temp.d.vector))] <- 0
    }
  } else {
    temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", tab(x), PACKAGE = "poppr")/2, 
                            temp.d.vector[, 1])
    temp.d.vector <- ceiling(temp.d.vector)
  }
  return(temp.d.vector)
}


jackbootplot <- function(df, obs.df, index = c("rd", "ia", "both"), obsline = TRUE){
  ARGS                <- c("rd", "ia", "both")
  index               <- match.arg(index, ARGS)
  Indexfac            <- c("I[A]","bar(r)[d]")
  levels(df$variable) <- Indexfac
  obs.df$variable     <- Indexfac

  if (index == "rd"){
    df     <- df[df$variable == "bar(r)[d]", ]
    obs.df <- obs.df[obs.df$variable == "bar(r)[d]", ]
  } else if (index == "ia"){
    df     <- df[df$variable == "I[A]", ]
    obs.df <- obs.df[obs.df$variable == "I[A]", ]
  }

  distplot <- ggplot(df, aes_string(y = "value", x = "Distribution")) +
                geom_violin(aes_string(fill = "Distribution")) + 
                geom_boxplot(alpha = 0.25, width = 0.125)
  if (index == "both"){
    distplot <- distplot + facet_grid("variable ~ .", scales = "free", 
                                      labeller = label_parsed) 
  } 
  if (obsline){
    aesth    <- aes_string(yintercept = "value", group = "variable") 
    distplot <- distplot + geom_hline(aesth, data = obs.df, linetype = 2)
  }
  return(distplot)      
}

#' @importFrom reshape2 melt

jackbootcomp <- function(pop, sample = 999, quiet = FALSE, method = 1, 
                         jack = 0.63, plotindex = "rd"){
  ARGS <- c("rd", "ia", "both")
  plotindex <- match.arg(plotindex, ARGS)
  
  inds <- nInd(pop)
  message("Null Distribution...")

  null <- ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
             hist = FALSE, method = method)
  if(pop@type == "codom"){
    popx <- seploc(pop)
  } else {
    popx <- pop
  }

  message("Alternative Distribution (Bootstrap)...")
  if(!quiet){
    bootbar <- txtProgressBar(style = 3)
  } else {
    bootbar <- NULL
  }
  jackbar <- bootbar

  altboot <- bootjack(popx, sample, half = NULL, bootbar)

  message("Alternative Distribution (Jack Knife)...")
  half    <- round(inds*jack)
  altjack <- bootjack(popx, sample, half, jackbar)  
  datlist <- list(null = null$sample, bootstrap = altboot, `jack knife` = altjack)
  dat           <- melt(datlist, measure.vars = c("Ia", "rbarD"))
  names(dat)[3] <- "Distribution"
  message("Creating Plots")
  obs.df   <- data.frame(list(variable = names(null$samples), 
                              value = null$index[c(1,3)]))
  distplot <- jackbootplot(dat, obs.df, plotindex) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals,", half, "Sampled for Jack knife")
                     )
  if (plotindex == "both"){
    distplot <- distplot + theme(strip.text = element_text(size = rel(2)))    
  } else if (plotindex == "rd"){
    distplot <- distplot + ylab(expression(bar(r)[d]))
  } else {
    distplot <- distplot + ylab(expression(I[A]))
  }
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



javier<-function(x){
  cat ("http://www.youtube.com/watch?v=1-ctsxVXvO0")
}
