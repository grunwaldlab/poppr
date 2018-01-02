
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
#' Internal function to plot the results from ia() and poppr()
#'
#' @param sample either an object of class "ialist" or a list of ialists
#' @param pval a named vector specifying the p values to display
#' @param pop The name of the population
#' @param file The name of the source file
#' @param N The number of samples in the population
#' @param observed observed values of Ia and rbarD
#' @param index The index to plot (defaults to "rbarD")
#' @param labsize size of the in-plot label
#' @param linesize size of the in-plot line
#'
#' @return a ggplot2 object
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' data(Pinf)
#' x <- Pinf %>% seppop() %>% lapply(ia, sample = 99, valuereturn = TRUE, quiet = TRUE, plot = FALSE)
#' x
#' poppr:::poppr.plot(sample = x, file = "hey") # plots multiple populations
#' # plot.ialist takes care of the single populations.
#' for (i in x){
#'   print(plot(i))
#' }
#' }
#' 
poppr.plot <- function(sample, pval = c(Ia = 0.05, rbarD = 0.05), 
                       pop = NULL, file = NULL, N = NULL,
                       observed = c(Ia = 0, rbarD = 0), 
                       index = c("rbarD", "Ia"),
                       labsize = rel(3),
                       linesize = rel(1)){
  INDEX_ARGS <- c("rbarD", "Ia")
  index      <- match.arg(index, INDEX_ARGS)
  tidy_ialist <- function(ialist){
    res <- stats::reshape(ialist$samples,
                          direction = "long",
                          varying = 1:2,
                          times = c("Ia", "rbarD"),
                          timevar = "variable",
                          v.names = "value")[-3] %>%
      dplyr::as_tibble()
    res$variable <- as.factor(res$variable)
    res
  }
  if (!class(sample) %in% "ialist" & class(sample) %in% "list"){
    ggsamps <- dplyr::bind_rows(lapply(sample, tidy_ialist), .id = "population")
    ggvals <- vapply(sample, "[[", numeric(4), "index")
    ggvals <-  setNames(as.data.frame.table(ggvals, stringsAsFactors = FALSE), 
                        c("index", "population", "value"))
    ggsamps <- ggsamps[ggsamps$variable == index, ]
    ggvals  <- ggvals[grep(ifelse(index == "Ia", "I", "D"), ggvals$index), ]
    ggindex <- ggvals[ggvals$index == index, ]
    ggpval  <- ggvals[ggvals$index == ifelse(index == "Ia", "p.Ia", "p.rD"), ]
    
    ggsamps$population <- factor(ggsamps$population, names(sample))
    srange  <- range(ggsamps$value, na.rm = TRUE)
    binw    <- diff(srange/30)
    plot_title <- make_poppr_plot_title(sample[[1]][["samples"]][[index]], 
                                        file = file)
    thePlot <- ggplot(ggsamps, aes_string(x = "value", group = "population")) +
      geom_histogram(binwidth = binw, position = "identity") +
      geom_rug(alpha = 0.5) + 
      geom_vline(aes_string(xintercept = "value"), color = "blue", linetype = 2,
                 data = ggindex) +
      facet_wrap(~population, scales = "free_x") +
      ggtitle(plot_title)
    if (index == "rbarD"){
      thePlot <- thePlot + xlab(expression(paste(bar(r)[d])))
    } else if (index == "Ia"){
      thePlot <- thePlot + xlab(expression(paste(I[A])))
    }
    return(thePlot)
  }
  plot_title <- make_poppr_plot_title(sample[[index]], file, N, pop)
  if (all(is.nan(sample$rbarD))){
    oops <- ggplot(as.data.frame(list(x=-10:9)), aes_string(x = "x")) + 
      geom_histogram(binwidth=1, fill="orange") + 
      geom_text(aes(label="Warning:", x=0, y=0.8), color="black", size=rel(15)) + 
      geom_text(aes(label="Data contains only NaNs and\ncannot be displayed graphically", 
                    x=0, y=0.5, hjust=0.5), color="black", size=rel(10)) +
      theme(line = element_blank(), 
            axis.text = element_blank(), 
            axis.title = element_blank()) +
      ggtitle(plot_title)
    return(oops)
  }
  if (any(is.na(sample[[index]]))){
    sample[[index]][is.na(sample[[index]])] <- mean(sample[[index]], na.rm=TRUE)
  }
  srange <- range(sample[[index]])
  labs   <- c(Ia = "I[A]", rbarD = "bar(r)[d]")
  names(pval) <- names(labs)
  binw <- diff(srange/30)
  if (is.na(observed[index])){
    xval <- srange[2]
    just <- 1
    vline <- geom_blank()
    annot_color <- "red"
    obs <- "NaN"
  } else if (all(observed[index] > srange)){
    xval <- observed[index]
    just <- 1.1
  } else {
    maxobsmin <- c(srange[1], observed[index], srange[2])
    xjust <- which.max(abs(diff(maxobsmin)))
    xval  <- srange[xjust]
    just  <- ifelse(xjust < 2, 0, 1)
  }
  if (!exists("annot_color")){
    annot_color <- "blue"
    vline       <- geom_vline(xintercept = observed[index], color = "blue", 
                              linetype = 2, size = linesize)
    obs         <- signif(observed[index], 3)
  }
  obslab  <- paste(labs[index], ":", obs, sep = "")
  plab    <- paste("p =", signif(pval[index], 3))
  ggdata  <- tidy_ialist(list(samples = sample))
  ggdata  <- ggdata[ggdata$variable == index, ] 
  thePlot <- ggplot(ggdata, aes_string(x = "value"))
  thePlot <- thePlot + 
    geom_histogram(binwidth = binw, position = "identity") +
    geom_rug() + 
    vline +
    annotate(geom = "text", x = xval, y = Inf, label = obslab, color = annot_color, 
             vjust = 2, hjust = just, parse = TRUE, size = labsize) + 
    annotate(geom = "text", x = xval, y = Inf, label = plab, color = annot_color, 
             vjust = 4, hjust = just, size = labsize)
  if (index == "rbarD"){
    thePlot <- thePlot + xlab(expression(paste(bar(r)[d])))
  } else if (index == "Ia"){
    thePlot <- thePlot + xlab(expression(paste(I[A])))
  }
  thePlot <- thePlot + ggtitle(plot_title)
  return(thePlot)
}

#==============================================================================#
#
#' Create a minimum spanning network of selected populations using a distance 
#' matrix.
#' 
#' @param gid a \code{\link{genind}}, \code{\link{genclone}},
#'   \code{\link{genlight}}, or \code{\link{snpclone}} object
#'   
#' @param distmat a distance matrix that has been derived from your data set.
#'   
#' @param mlg.compute if the multilocus genotypes are set to "custom" (see 
#'   \code{\link{mll.custom}} for details) in your genclone object, this will 
#'   specify which mlg level to calculate the nodes from. See details.
#'   
#' @param palette a \code{vector} or \code{function} defining the color palette 
#'   to be used to color the populations on the graph. It defaults to 
#'   \code{\link{topo.colors}}. See examples for details.
#'   
#' @param sublist a \code{vector} of population names or indexes that the user 
#'   wishes to keep. Default to "ALL".
#'   
#' @param blacklist a \code{vector} of population names or indexes that the user
#'   wishes to discard. Default to \code{NULL}
#'   
#' @param vertex.label a \code{vector} of characters to label each vertex. There
#'   are two defaults: \code{"MLG"} will label the nodes with the multilocus 
#'   genotype from the original data set and \code{"inds"} will label the nodes 
#'   with the representative individual names.
#'   
#' @param gscale "grey scale". If this is \code{TRUE}, this will scale the color
#'   of the edges proportional to the observed distance, with the lines becoming
#'   darker for more related nodes. See \code{\link{greycurve}} for details.
#'   
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#'   the upper and lower limits for the \code{\link{gray}} function. Default is 
#'   0 (black) and 0.8 (20\% black). See \code{\link{greycurve}} for details.
#'   
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that 
#'   will serve as the exponent to the edge weight to scale the grey value to 
#'   represent that weight. See \code{\link{greycurve}} for details.
#'   
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale 
#'   will be weighted to emphasize the differences between closely related 
#'   nodes. If it is 2, the grey scale will be weighted to emphasize the 
#'   differences between more distantly related nodes. See 
#'   \code{\link{greycurve}} for details.
#'   
#' @param wscale "width scale". If this is \code{TRUE}, the edge widths will be 
#'   scaled proportional to the inverse of the observed distance , with the 
#'   lines becoming thicker for more related nodes.
#'   
#' @param showplot logical. If \code{TRUE}, the graph will be plotted. If 
#'   \code{FALSE}, it will simply be returned.
#'   
#' @param include.ties logical. If \code{TRUE}, the graph will include all edges
#'   that were arbitrarily passed over in favor of another edge of equal weight.
#'   If \code{FALSE}, which is the default, one edge will be arbitrarily 
#'   selected when two or more edges are tied, resulting in a pure minimum
#'   spanning network.
#'   
#' @param threshold numeric. By default, this is \code{NULL}, which will have no
#'   effect. Any threshold value passed to this argument will be used in
#'   \code{\link{mlg.filter}} prior to creating the MSN. If you have a data set
#'   that contains contracted MLGs, this argument will override the threshold in
#'   the data set. See Details.
#'
#' @param clustering.algorithm string. By default, this is \code{NULL}. If
#'   \code{threshold = NULL}, this argument will have no effect. When supplied
#'   with either "farthest_neighbor", "average_neighbor", or "nearest_neighbor",
#'   it will be passed to \code{\link{mlg.filter}} prior to creating the MSN. If
#'   you have a data set that contains contracted MLGs, this argument will
#'   override the algorithm in the data set. See Details.
#'
#' @param ... any other arguments that could go into plot.igraph
#'   
#' @return \item{graph}{a minimum spanning network with nodes corresponding to 
#'   MLGs within the data set. Colors of the nodes represent population 
#'   membership. Width and color of the edges represent distance.} 
#'   \item{populations}{a vector of the population names corresponding to the 
#'   vertex colors} \item{colors}{a vector of the hexadecimal representations of
#'   the colors used in the vertex colors}
#'   
#'   
#' @details The minimum spanning network generated by this function is generated
#'   via igraph's \code{\link[igraph]{minimum.spanning.tree}}. The resultant
#'   graph produced can be plotted using igraph functions, or the entire object
#'   can be plotted using the function \code{\link{plot_poppr_msn}}, which will
#'   give the user a scale bar and the option to layout your data.
#'   \subsection{node sizes}{
#'   The area of the nodes are representative of the number of samples. Because
#'   \pkg{igraph} scales nodes by radius, the node sizes in the graph are 
#'   represented as the square root of the number of samples.}
#'   \subsection{mlg.compute}{
#'   Each node on the graph represents a different multilocus genotype. 
#'   The edges on the graph represent genetic distances that connect the
#'   multilocus genotypes. In genclone objects, it is possible to set the
#'   multilocus genotypes to a custom definition. This creates a problem for
#'   clone correction, however, as it is very possible to define custom lineages
#'   that are not monophyletic. When clone correction is performed on these
#'   definitions, information is lost from the graph. To circumvent this, The
#'   clone correction will be done via the computed multilocus genotypes, either
#'   "original" or "contracted". This is specified in the \code{mlg.compute}
#'   argument, above.}
#'   \subsection{contracted multilocus genotypes}{
#'   If your incoming data set is of the class \code{\linkS4class{genclone}},
#'   and it contains contracted multilocus genotypes, this function will retain
#'   that information for creating the minimum spanning network. You can use the
#'   arguments \code{threshold} and \code{clustering.algorithm} to change the
#'   threshold or clustering algorithm used in the network. For example, if you
#'   have a data set that has a threshold of 0.1 and you wish to have a minimum
#'   spanning network without a threshold, you can simply add 
#'   \code{threshold = 0.0}, and no clustering will happen. 
#'   
#'   The \code{threshold} and \code{clustering.algorithm} arguments can also be
#'   used to filter un-contracted data sets.
#'   
#'   All filtering will use the distance matrix supplied in the argument 
#'   \code{distmat}.
#'   }
#'   
#' @note The edges of these graphs may cross each other if the graph becomes too
#'   large.
#'   
#' @seealso \code{\link{plot_poppr_msn}} \code{\link{nancycats}},
#'   \code{\link{upgma}}, \code{\link{nj}}, \code{\link{nodelabels}},
#'   \code{\link{tab}}, \code{\link{missingno}}, \code{\link{bruvo.msn}},
#'   \code{\link{greycurve}}
#' 
#' @export
#' @aliases msn.poppr
#' @author Javier F. Tabima, Zhian N. Kamvar, Jonah C. Brooks
#' @examples
#' 
#' # Load the data set and calculate the distance matrix for all individuals.
#' data(Aeut)
#' A.dist <- diss.dist(Aeut)
#' 
#' # Graph it.
#' A.msn <- poppr.msn(Aeut, A.dist, gadj = 15, vertex.label = NA)
#' 
#' # Find the sizes of the nodes (number of individuals per MLL):
#' igraph::vertex_attr(A.msn$graph, "size")^2
#' 
#' \dontrun{
#' # Set subpopulation structure.
#' Aeut.sub <- as.genclone(Aeut)
#' setPop(Aeut.sub) <- ~Pop/Subpop
#' 
#' # Plot respective to the subpopulation structure
#' As.msn <- poppr.msn(Aeut.sub, A.dist, gadj=15, vertex.label=NA)
#' 
#' # Show only the structure of the Athena population.
#' As.msn <- poppr.msn(Aeut.sub, A.dist, gadj=15, vertex.label=NA, sublist=1:10)
#' 
#' # Let's look at the structure of the microbov data set
#'
#' library("igraph")
#' data(microbov)
#' micro.dist <- diss.dist(microbov, percent = TRUE)
#' micro.msn <- poppr.msn(microbov, micro.dist, vertex.label=NA)
#' 
#' # Let's plot it and show where individuals have < 15% of their genotypes 
#' # different.
#' 
#' edge_weight <- E(micro.msn$graph)$weight
#' edge_labels <- ifelse(edge_weight < 0.15, round(edge_weight, 3), NA)
#' plot.igraph(micro.msn$graph, edge.label = edge_labels, vertex.size = 2, 
#' edge.label.color = "red")
#'
#' }
#' 
#==============================================================================#
poppr.msn <- function (gid, distmat, palette = topo.colors, mlg.compute = "original",
                       sublist = "All", blacklist = NULL, vertex.label = "MLG", 
                       gscale=TRUE, glim = c(0,0.8), gadj = 3, gweight = 1, 
                       wscale=TRUE, showplot = TRUE, 
                       include.ties = FALSE, threshold = NULL, 
                       clustering.algorithm = NULL, ...){
  
  # testing gid -------------------------------------------------------------
  if (!inherits(gid, c("genclone", "snpclone"))){
    if (inherits(gid, "genind")){
      gid <- as.genclone(gid)
    } else if (inherits(gid, "genlight")){
      gid <- as.snpclone(gid)
    } else {
      stop("gid must be a genind/genclone or genlight/snpclone object.")
    }
  }
  if (!inherits(gid@mlg, "MLG")){
    gid@mlg <- new("MLG", gid@mlg)
  }
  
  
  # Handling MLG type -------------------------------------------------------
  visible_mlg <- visible(gid@mlg)
  if (visible_mlg == "custom"){
    mll(gid) <- mlg.compute
  } else if (visible_mlg == "contracted"){
    mll(gid) <- "original"
    if (is.null(threshold)){
      threshold <- cutoff(gid@mlg)["contracted"]
    }
    
    if (is.null(clustering.algorithm)){
      clustering.algorithm <- distalgo(gid@mlg)
    } 
    
  }
  
  # testing dist ------------------------------------------------------------
  is_dist <- inherits(distmat, "dist")
  is_mat  <- inherits(distmat, "matrix")
  if (is_dist | is_mat){
    n       <- nInd(gid)
    eq_size <- if (is_dist) n == attr(distmat, "Size") else n == nrow(distmat)
    if (!eq_size){
      stop("The size of the distance matrix does not match the size of the data.\n")
    }
  } else {
    stop("The distance matrix is neither a dist object nor a matrix.\n")
  }
  distmat <- as.matrix(distmat)
  gadj   <- ifelse(gweight == 1, gadj, -gadj)
  
  # Subsetting the population -----------------------------------------------
  # This will subset both the population and the matrix. 
  if (toupper(sublist[1]) != "ALL" | !is.null(blacklist)){
    sublist_blacklist <- sub_index(gid, sublist, blacklist)
    distmat <- distmat[sublist_blacklist, sublist_blacklist, drop = FALSE]
    gid    <- popsub(gid, sublist, blacklist)
  }

  # Clone correcting the matrix ---------------------------------------------
  if (!is.null(threshold)){
    filtered <- filter_at_threshold(gid, 
                                    threshold, 
                                    indist = distmat,
                                    clustering.algorithm,
                                    bruvo_args = NULL)
    distmat <- filtered$indist
    cgid    <- filtered$cgid
    gid     <- filtered$gid
  } else {  
    cgid    <- gid[.clonecorrector(gid), ]
    singles <- !duplicated(mll(gid))
    distmat <- distmat[singles, singles, drop = FALSE]
  }
  rownames(distmat) <- indNames(cgid) -> colnames(distmat)
  poppr_msn_list <- msn_constructor(
    gid = gid,
    cgid = cgid,
    palette = palette,
    indist = distmat,
    include.ties = include.ties,
    mlg.compute = mlg.compute,
    vlab = vertex.label,
    visible_mlg = visible_mlg,
    wscale = wscale,
    gscale = gscale,
    glim = glim,
    gadj = gadj,
    showplot = showplot,
    ...)
  return(poppr_msn_list)
}




#==============================================================================#
#' Create a table summarizing missing data or ploidy information of a genind or
#' genclone object
#' 
#' @param gen a \linkS4class{genind} or \linkS4class{genclone} object.
#'   
#' @param type \code{character}. What information should be returned. Choices
#'   are "missing" (Default) and "ploidy". See Description.
#'   
#' @param percent \code{logical}. (ONLY FOR \code{type = 'missing'}) If
#'   \code{TRUE} (default), table and plot will represent missing data as a
#'   percentage of each cell. If \code{FALSE}, the table and plot will represent
#'   missing data as raw counts. (See details)
#'   
#' @param plot \code{logical}. If \code{TRUE}, a simple heatmap will be 
#'   produced. If \code{FALSE} (default), no heatmap will be produced.
#'   
#' @param df \code{logical}. If \code{TRUE}, the data will be returned as a long
#'   form data frame. If \code{FALSE} (default), a matrix with samples in rows
#'   and loci in columns will be returned.
#'   
#' @param returnplot \code{logical}. If \code{TRUE}, a list is returned with two
#'   elements: \code{table} - the normal output and \code{plot} - the ggplot 
#'   object. If \code{FALSE}, the table is returned.
#'   
#' @param low \code{character}. What color should represent no missing data or 
#'   lowest observed ploidy? (default: "blue")
#'   
#' @param high \code{character}. What color should represent the highest amount 
#'   of missing data or observed ploidy? (default: "red")
#'   
#' @param plotlab \code{logical}. (ONLY FOR \code{type = 'missing'}) If
#'   \code{TRUE} (default), values of missing data greater than 0\% will be
#'   plotted. If \code{FALSE}, the plot will appear un-appended.
#'   
#' @param scaled \code{logical}. (ONLY FOR \code{type = 'missing'}) This is for
#'   when \code{percent = TRUE}. If \code{TRUE} (default), the color specified
#'   in \code{high} will represent the highest observed value of missing data.
#'   If \code{FALSE}, the color specified in \code{high} will represent 100\%.
#'   
#' @return a matrix, data frame (\code{df = TRUE}), or a list (\code{returnplot 
#'   = TRUE}) representing missing data per population (\code{type = 'missing'})
#'   or ploidy per individual (\code{type = 'ploidy'}) in a \linkS4class{genind}
#'   or \linkS4class{genclone} object.
#' 
#' @details 
#'   Missing data is accounted for on a per-population level.\cr
#'   Ploidy is accounted for on a per-individual level.
#'   
#'   \subsection{For type = 'missing'}{
#'   This data is potentially useful for identifying areas of systematic missing
#'   data. There are a few caveats to be aware of. \itemize{ \item
#'   \strong{Regarding counts of missing data}: Each count represents the number
#'   of individuals with missing data at each locus. The last column, "mean" can
#'   be thought of as the average number of individuals with missing data per
#'   locus. \item \strong{Regarding percentage missing data}: This percentage is
#'   \strong{relative to the population and locus}, not to the entire data set.
#'   The last column, "mean" represents the average percent of the population
#'   with missing data per locus. }} 
#'   \subsection{For type = 'ploidy'}{
#'   This option is useful for data that has been imported with mixed ploidies.
#'   It will summarize the relative levels of ploidy per individual per locus.
#'   This is simply based off of observed alleles and does not provide any
#'   further estimates.}
#'   
#' @export
#' @keywords missing ploidy
#' @author Zhian N. Kamvar
#' @examples
#' data(nancycats)
#' nancy.miss <- info_table(nancycats, plot = TRUE, type = "missing")
#' data(Pinf)
#' Pinf.ploid <- info_table(Pinf, plot = TRUE, type = "ploidy")
#' 
#==============================================================================#
info_table <- function(gen, type = c("missing", "ploidy"), percent = TRUE, plot = FALSE, 
                       df = FALSE, returnplot = FALSE, low = "blue", 
                       high = "red", plotlab = TRUE, scaled = TRUE){
  datalabel <- as.character(match.call()[2])
  ARGS      <- c("missing", "ploidy")
  type      <- match.arg(type, ARGS)

  if (type == "missing"){

    valname    <- "Missing"
    pops       <- seppop(gen, drop = FALSE)
    pops$Total <- gen
    inds       <- 1
    if (percent){
      inds <- c(table(pop(gen)), nInd(gen))
    }
    data_table <- matrix(0, nrow = nLoc(gen) + 1, ncol = length(pops))
    data_table[1:nLoc(gen), ] <- vapply(pops, number_missing_locus, numeric(nLoc(gen)), 1)
    data_table[-nrow(data_table), ] <- t(apply(data_table[-nrow(data_table), ], 1, "/", inds))
    data_table[nrow(data_table), ]  <- colMeans(data_table[-nrow(data_table), ])
    rownames(data_table)         <- c(locNames(gen), "Mean")
    colnames(data_table)         <- names(pops)
    dimnames(data_table) <- list(Locus = c(locNames(gen), "Mean"), Population = names(pops))
    if (all(data_table == 0)){
      message("No Missing Data Found!")
      return(NULL)
    }
    if (plot){
      data_df      <- as.data.frame.table(data_table, 
                                          responseName = valname,
                                          stringsAsFactors = FALSE)
      leg_title    <- valname
      data_df[1:2] <- data.frame(lapply(data_df[1:2], 
                                       function(x) factor(x, levels = unique(x))))
      plotdf <- data_df -> textdf 
      if (percent) {
        plotdf$Missing <- round(plotdf$Missing*100, 2)
        textdf$Missing <- paste(plotdf$Missing, "%")
        miss           <- "0 %"
        title          <- paste("Percent missing data per locus and population of", 
                                datalabel)
        leg_title      <- paste("Percent", leg_title)
        if(!scaled){
          lims <- c(0, 100)
        }
      } else {
        textdf$Missing <- round(textdf$Missing, 2)
        miss           <- 0
        title          <- paste("Missing data per locus and population of", 
                                datalabel)
      }
      if (scaled | !percent){
        lims <- c(0, max(plotdf$Missing))
      }
      linedata <- data.frame(list(yint = 1.5, xint = nrow(data_table) - 0.5))
      textdf$Missing <- ifelse(textdf$Missing == miss, "", textdf$Missing)
      plotdf$Missing[plotdf$Locus == "Mean" & plotdf$Population == "Total"] <- NA

      outplot <- ggplot(plotdf, aes_string(x = "Locus", y = "Population")) + 
        geom_tile(aes_string(fill = valname)) +
        labs(list(title = title, x = "Locus", y = "Population")) +
        labs(fill = leg_title) + 
        scale_fill_gradient(low = low, high = high, na.value = "white", 
                            limits = lims) +
        geom_hline(aes_string(yintercept = "yint"), data = linedata) + 
        geom_vline(aes_string(xintercept = "xint"), data = linedata) 
      if (plotlab){
        outplot <- outplot + geom_text(aes_string(label = valname), 
                                       data = textdf)
      }
      outplot <- outplot +
        theme_classic() + 
        scale_x_discrete(expand = c(0, -1)) + 
        scale_y_discrete(expand = c(0, -1), 
                         limits = rev(unique(plotdf$Population))) + 
        theme(axis.text.x = element_text(size = 10, angle = -45, 
                                         hjust = 0, vjust = 1)) 
      print(outplot)
    }

  } else if (type == "ploidy"){

    valname    <- "Observed_Ploidy"
    data_table <- get_local_ploidy(gen)
    
    dimnames(data_table) <- list(Samples = indNames(gen), Loci = locNames(gen))
    if (plot){
      data_df      <- as.data.frame.table(data_table, 
                                          responseName = valname,
                                          stringsAsFactors = FALSE)
      data_df[1:2] <- data.frame(lapply(data_df[1:2], 
                                       function(x) factor(x, levels = unique(x))))
      vars <- aes_string(x = "Loci", y = "Samples", fill = valname)

      mytheme <- theme_classic() +  
                 theme(axis.text.x = element_text(size = 10, angle = -45, 
                                                  hjust = 0, vjust = 1)) 
      legvars <- unique(data_df[[valname]])
      title <- paste("Observed ploidy of", datalabel)
      outplot <- ggplot(data_df) + geom_tile(vars) + 
                   scale_fill_gradient(low = low, high = high, guide = "legend", 
                                       breaks = legvars) + 
                   scale_x_discrete(expand = c(0, -1)) + 
                   scale_y_discrete(expand = c(0, -1), 
                                    limits = rev(unique(data_df$Samples))) + 
                   labs(list(title = title, x = "Locus", y = "Sample", 
                             fill = "Observed\nPloidy")) +
                   mytheme 

      print(outplot)
    }
  } 
  if (df){
    if (!exists("data_df")){
      data_df <- as.data.frame.table(data_table, 
                                     responseName = valname,
                                     stringsAsFactors = FALSE)
    }
    data_table <- data_df
  } else {
    if (type == "missing"){
      data_table <- t(data_table)
    }
    class(data_table) <- c("locustable", "matrix")
  }
  if (returnplot & exists("outplot")){
    data_table <- list(table = data_table, plot = outplot)
  }
  return(data_table)
}

#==============================================================================#
#' Display a greyscale gradient adjusted to specific parameters
#' 
#' This function has one purpose. It is for deciding the appropriate scaling for
#' a grey palette to be used for edge weights of a minimum spanning network.
#'
#'
#' @param data a sequence of numbers to be converted to greyscale.
#' 
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#' the upper and lower limits for the \code{\link{gray}} function. Default is 0
#' (black) and 0.8 (20\% black). 
#' 
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that
#' will serve as the exponent to the edge weight to scale the grey value to
#' represent that weight.
#' 
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale
#' will be weighted to emphasize the differences between closely related nodes.
#' If it is 2, the grey scale will be weighted to emphasize the differences
#' between more distantly related nodes. 
#'
#' @param scalebar When this is set to \code{TRUE}, two scalebars will be
#' plotted. The purpose of this is for adding a scale bar to minimum spanning
#' networks produced in earlier versions of poppr. 
#' 
#' @return A plot displaying a grey gradient from 0.001 to 1 with minimum and 
#' maximum values displayed as yellow lines, and an equation for the correction 
#' displayed in red. 
#' 
#' @author Zhian N. Kamvar
#'
#' @examples
#' # Normal grey curve with an adjustment of 3, an upper limit of 0.8, and
#' # weighted towards smaller values.
#' greycurve()
#' \dontrun{
#' # 1:1 relationship grey curve.
#' greycurve(gadj=1, glim=1:0)
#' 
#' # Grey curve weighted towards larger values.
#' greycurve(gweight=2)
#' 
#' # Same as the first, but the limit is 1.
#' greycurve(glim=1:0)
#' 
#' # Setting the lower limit to 0.1 and weighting towards larger values.
#' greycurve(glim=c(0.1,0.8), gweight=2)
#' }
#' @export
#==============================================================================#
greycurve <- function(data = seq(0, 1, length = 1000), glim = c(0,0.8), 
                      gadj = 3, gweight = 1, scalebar = FALSE){
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  adjustcurve(data, glim, correction = gadj, show = TRUE, scalebar = scalebar)
}


#==============================================================================#
#' Plot minimum spanning networks produced in poppr.
#' 
#' This function allows you to take the output of poppr.msn and bruvo.msn and 
#' customize the plot by labeling groups of individuals, size of nodes, and 
#' adjusting the palette and scale bar.
#' 
#' @param x a \code{\linkS4class{genind}}, \code{\linkS4class{genclone}},
#'   \code{\linkS4class{genlight}}, or \code{\linkS4class{snpclone}} object from
#'   which \code{poppr_msn} was derived.
#'   
#' @param poppr_msn a \code{list} produced from either \code{\link{poppr.msn}} 
#'   or \code{\link{bruvo.msn}}. This list should contain a graph, a vector of 
#'   population names and a vector of hexadecimal color definitions for each 
#'   population.
#'   
#' @inheritParams greycurve
#' 
#' @inheritParams poppr.msn
#'   
#'   
#' @param nodescale a \code{numeric} indicating how to scale the node sizes (scales by area).
#' @param nodebase \strong{deprecated} a \code{numeric} indicating what base logarithm should be
#'   used to scale the node sizes. Defaults to 1.15. See details.
#'   
#' @param nodelab an \code{integer} specifying the smallest size of node to 
#'   label. See details.
#'   
#' @param inds a \code{character} or \code{numeric} vector indicating which
#'   samples or multilocus genotypes to label on the graph. See details.
#'   
#' @param mlg \code{logical} When \code{TRUE}, the nodes will be labeled by
#'   multilocus genotype. When \code{FALSE} (default), nodes will be labeled by
#'   sample names.
#'   
#' @param quantiles \code{logical}. When set to \code{TRUE} (default), the scale
#'   bar will be composed of the quantiles from the observed edge weights. When 
#'   set to \code{FALSE}, the scale bar will be composed of a smooth gradient 
#'   from the minimum edge weight to the maximum edge weight.
#'   
#' @param cutoff a number indicating the longest distance to display in your 
#'   graph. This is performed by removing edges with weights greater than this 
#'   number.
#'   
#' @param palette a function or character corresponding to a specific palette 
#'   you want to use to delimit your populations. The default is whatever 
#'   palette was used to produce the original graph.
#'   
#' @param layfun a function specifying the layout of nodes in your graph. It 
#'   defaults to \code{\link[igraph]{layout.auto}}.
#'   
#' @param beforecut if \code{TRUE}, the layout of the graph will be computed 
#'   before any edges are removed with \code{cutoff}. If \code{FALSE} (Default),
#'   the layout will be computed after any edges are removed.
#'   
#' @param pop.leg if \code{TRUE}, a legend indicating the populations will
#'   appear in the top right corner of the graph, but will not overlap. Setting
#'   \code{pop.leg = FALSE} disables this legend. See details.
#' 
#' @param size.leg if \code{TRUE}, a legend displyaing the number of samples per
#'   node will appear either below the population legend or in the top right
#'   corner of the graph. Setting \code{size.leg = FALSE} disables this legend.
#'   
#' @param scale.leg if \code{TRUE}, a scale bar indicating the distance will
#'   appear under the graph. Setting \code{scale.leg = FALSE} suppresses this
#'   bar. See details.
#'   
#' @param ... any other parameters to be passed on to 
#'   \code{\link[igraph]{plot.igraph}}.
#'   
#' @details The previous incarnation of msn plotting in poppr simply plotted the
#'   minimum spanning network with the legend of populations, but did not 
#'   provide a scale bar and it did not provide the user a simple way of 
#'   manipulating the layout or labels. This function allows the user to 
#'   manipulate many facets of graph creation, making the creation of minimum 
#'   spanning networks ever so slightly more user friendly. 
#'   
#'   This function must have both the source data and the output msn to work. 
#'   The source data must contain the same population structure as the graph. 
#'   Every other parameter has a default setting.
#'   
#'   \subsection{Parameter details}{ \itemize{ 
#'   \item \code{inds} By default, the graph will label each node (circle) with
#'   all of the samples (individuals) that are contained within that node. As
#'   each node represents a single multilocus genotype (MLG) or individuals (n
#'   >= 1), this argument is designed to allow you to selectively label the
#'   nodes based on query of sample name or MLG number. If the option \code{mlg
#'   = TRUE}, the multilocus genotype assignment will be used to label the node.
#'   If you do not want to label the nodes by individual or multilocus genotype,
#'   simply set this to a name that doesn't exist in your data.
#'   \item \code{nodescale} The nodes (circles) on the graph represent different
#'   multilocus genotypes. The area of the nodes represent the number of
#'   individuals. Setting nodescale will scale the area of the nodes.
#'   \item \code{nodelab} If a node is not labeled by individual, this will
#'   label the size of the nodes greater than or equal to this value. If you
#'   don't want to label the size of the nodes, simply set this to a very high
#'   number.
#'   \item \code{cutoff} This is useful for when you want to investigate groups
#'   of multilocus genotypes separated by a specific distance or if you have two
#'   distinct populations and you want to physically separate them in your
#'   network.
#'   \item \code{beforecut} This is an indicator useful if you want to maintain
#'   the same position of the nodes before and after removing edges with the
#'   \code{cutoff} argument. This works best if you set a seed before you run
#'   the function.
#'   }}
#'   
#'   \subsection{mlg.compute}{
#'   Each node on the graph represents a different multilocus genotype. 
#'   The edges on the graph represent genetic distances that connect the
#'   multilocus genotypes. In genclone objects, it is possible to set the
#'   multilocus genotypes to a custom definition. This creates a problem for
#'   clone correction, however, as it is very possible to define custom lineages
#'   that are not monophyletic. When clone correction is performed on these
#'   definitions, information is lost from the graph. To circumvent this, The
#'   clone correction will be done via the computed multilocus genotypes, either
#'   "original" or "contracted". This is specified in the \code{mlg.compute}
#'   argument, above.}
#'   
#'   \subsection{legends}{ To avoid drawing the legend over the graph, legends 
#'   are separated by different plotting areas. This means that if legends are 
#'   included, you cannot plot multiple MSNs in a single plot. The scale bar (to
#'   be added in manually) can be obtained from \code{\link{greycurve}} and the
#'   legend can be plotted with \code{\link{legend}}.}
#' 
#' @return the modified msn list, invisibly. 
#' 
#' @seealso \code{\link[igraph]{layout.auto}} \code{\link[igraph]{plot.igraph}}
#' \code{\link{poppr.msn}} \code{\link{bruvo.msn}} \code{\link{greycurve}}
#' \code{\link[igraph]{delete.edges}} \code{\link{palette}}
#' 
#' @author Zhian N. Kamvar
#' @export
#' @examples
#' # Using a data set of the Aphanomyces eutieches root rot pathogen.
#' data(Aeut)
#' adist <- diss.dist(Aeut, percent = TRUE)
#' amsn <- poppr.msn(Aeut, adist, showplot = FALSE)
#' 
#' # Default
#' library("igraph") # To get all the layouts.
#' set.seed(500)
#' plot_poppr_msn(Aeut, amsn, gadj = 15)
#' 
#' \dontrun{
#' 
#' # Different layouts (from igraph) can be used by supplying the function name.
#' set.seed(500)
#' plot_poppr_msn(Aeut, amsn, gadj = 15, layfun = layout_with_kk)
#' 
#' # Removing link between populations (cutoff = 0.2) and labelling no individuals
#' set.seed(500)
#' plot_poppr_msn(Aeut, amsn, inds = "none", gadj = 15, beforecut = TRUE, cutoff = 0.2)
#' 
#' # Labelling individual #57 because it is an MLG that crosses populations
#' # Showing clusters of MLGS with at most 5% variation
#' # Notice that the Mt. Vernon population appears to be more clonal
#' set.seed(50) 
#' plot_poppr_msn(Aeut, amsn, gadj = 15, cutoff = 0.05, inds = "057")
#' 
#' 
#' data(partial_clone)
#' pcmsn <- bruvo.msn(partial_clone, replen = rep(1, 10))
#' 
#' # You can plot using a color palette or a vector of named colors
#' # Here's a way to define the colors beforehand
#' pc_colors <- nPop(partial_clone) %>% 
#'   RColorBrewer::brewer.pal("Set2") %>% 
#'   setNames(popNames(partial_clone))
#'
#' pc_colors
#' 
#' # Labelling the samples contained in multilocus genotype 9
#' set.seed(999)
#' plot_poppr_msn(partial_clone, pcmsn, palette = pc_colors, inds = 9)
#' 
#' # Doing the same thing, but using one of the sample names as input.
#' set.seed(999)
#' plot_poppr_msn(partial_clone, pcmsn, palette = pc_colors, inds = "sim 20")
#' 
#' # Note that this is case sensitive. Nothing is labeled. 
#' set.seed(999)
#' plot_poppr_msn(partial_clone, pcmsn, palette = pc_colors, inds = "Sim 20")
#' 
#' # Something pretty
#' data(microbov)
#' mdist <- diss.dist(microbov, percent = TRUE)
#' micmsn <- poppr.msn(microbov, mdist, showplot = FALSE)
#' 
#' plot_poppr_msn(microbov, micmsn, palette = "terrain.colors", inds = "n", 
#'   quantiles = FALSE)
#' plot_poppr_msn(microbov, micmsn, palette = "terrain.colors", inds = "n", 
#'   cutoff = 0.3, quantiles = FALSE)
#'   
#' ### Utilizing vectors for palettes
#' 
#' data(Pram)
#' Pram_sub <- popsub(Pram, blacklist = c("Nursery_CA", "Nursery_OR"))
#' 
#' # Creating the network for the forest
#' min_span_net_sub <- bruvo.msn(Pram_sub, replen = other(Pram)$REPLEN, 
#'                               add = TRUE, loss = TRUE, showplot = FALSE, 
#'                               include.ties = TRUE)
#'                               
#' # Creating the network with nurseries
#' min_span_net     <- bruvo.msn(Pram, replen = other(Pram)$REPLEN, 
#'                               add = TRUE, loss = TRUE, showplot = FALSE, 
#'                               include.ties = TRUE)
#'
#' # Only forest genotypes
#' set.seed(70)
#' plot_poppr_msn(Pram,
#'                min_span_net_sub,
#'                inds = "ALL",
#'                mlg = TRUE,
#'                gadj = 9,
#'                nodescale = 5,
#'                palette = other(Pram)$comparePal,
#'                cutoff = NULL,
#'                quantiles = FALSE,
#'                beforecut = TRUE)
#' 
#' # With Nurseries
#' set.seed(70)
#' plot_poppr_msn(Pram,
#'                min_span_net,
#'                inds = "ALL",
#'                mlg = TRUE,
#'                gadj = 9,
#'                nodescale = 5,
#'                palette = other(Pram)$comparePal,
#'                cutoff = NULL,
#'                quantiles = FALSE,
#'                beforecut = TRUE)
#' }
#==============================================================================#
#' @importFrom igraph layout.auto delete.edges
plot_poppr_msn <- function(x,
                           poppr_msn,
                           gscale = TRUE,
                           gadj = 3,
                           mlg.compute = "original",
                           glim = c(0, 0.8),
                           gweight = 1,
                           wscale = TRUE,
                           nodescale = 10,
                           nodebase = NULL,
                           nodelab = 2,
                           inds = "ALL",
                           mlg = FALSE,
                           quantiles = TRUE,
                           cutoff = NULL,
                           palette = NULL,
                           layfun = layout.auto,
                           beforecut = FALSE,
                           pop.leg = TRUE,
                           size.leg = TRUE,
                           scale.leg = TRUE,
                           ...) {
  
  if (!is(x, "genlight") && !is.genind(x)){
    stop(paste(substitute(x), "is not a genind or genclone object."))
  }
  if (!identical(names(poppr_msn), c("graph", "populations", "colors"))){
    stop("graph not compatible")
  }
  if (!is.null(palette)){
    pal       <- palette
    newpal    <- palette_parser(pal, length(poppr_msn$populations), poppr_msn$populations)
    poppr_msn <- update_poppr_graph(poppr_msn, function(x) newpal)
  }
  # Making sure incoming data matches so that the individual names match.
  if (!all(is.na(poppr_msn$populations))){
    if (nPop(x) > 0 && nPop(x) >= length(poppr_msn$populations)){
      x <- popsub(x, sublist = poppr_msn$populations)
    } else {
      warning("populations in graph don't match data. Setting to none.")
    }
  }
  
  if (beforecut){
    LAYFUN <- match.fun(layfun)
    lay <- LAYFUN(poppr_msn$graph)
  } else {
    lay <- match.fun(layfun)
  }
  # delete.edges      <- match.fun(igraph::delete.edges)
  if (!is.null(cutoff) && !is.na(cutoff)){
    if (all(cutoff < E(poppr_msn$graph)$weight)){
      msg <- paste0("Cutoff value (", cutoff, ") is below the minimum observed",
                    " distance. Edges will not be removed.")
      warning(msg)
    } else {
      E_above_cutoff  <- E(poppr_msn$graph)[E(poppr_msn$graph)$weight >= cutoff]
      poppr_msn$graph <- delete.edges(poppr_msn$graph, E_above_cutoff)
    }
  }
  # Adjusting color scales. This will replace any previous scaling contained in
  # poppr_msn.
  weights <- E(poppr_msn$graph)$weight
  wmin    <- min(weights)
  wmax    <- max(weights)
  gadj    <- ifelse(gweight == 1, gadj, -gadj)
  poppr_msn$graph <- update_edge_scales(poppr_msn$graph, wscale, gscale, glim, gadj)
  
  x.mlg <- mlg.vector(x)
  if (is.factor(x.mlg)){
    x.mlg <- mll(x, mlg.compute)
    custom <- TRUE
    # labs <- unique(x.mlg)
    labs  <- correlate_custom_mlgs(x, mlg.compute)
  } else {
    labs  <- unique(x.mlg)    
    custom <- FALSE
  }

  

  # Highlighting only the names of the submitted genotypes and the matching
  # isolates.
  
  # The labels in the graph are organized by MLG, so we will use that to
  # extract the names we need.
  if (length(inds) == 1 & toupper(inds[1]) == "ALL"){
    x.input <- unique(x.mlg)
  } else if (is.character(inds)){
    x.input <- unique(x.mlg[indNames(x) %in% inds])
  } else if (is.numeric(inds)){
    x.input <- unique(x.mlg[x.mlg %in% inds])
  }
  
  # Remove labels that are not specified.
  if (!custom){
    labs[!labs %in% x.input] <- NA    
  } else {
    labs[!labs %in% labs[as.character(x.input)]] <- NA
  }


  if (!isTRUE(mlg)){
    # Combine all the names that match with each particular MLG in x.input.
    combined_names <- vapply(x.input, function(mlgname)
                             paste(rev(indNames(x)[x.mlg == mlgname]),
                                   collapse = "\n"),
                             character(1))
    labs[!is.na(labs)] <- combined_names
  } else if (is.numeric(x.mlg) && !custom){
    labs[!is.na(labs)] <- paste("MLG", labs[!is.na(labs)], sep = ".")
  } else {
    labs <- as.character(labs)
  }
  if (any(is.na(labs))){
    sizelabs <- V(poppr_msn$graph)$size
    sizelabs <- ifelse(sizelabs >= nodelab, sizelabs * sizelabs, NA)
    labs     <- ifelse(is.na(labs), sizelabs, labs)
  }

  if (!is.null(nodebase)){
    warnw    <- min(.9 * getOption("width"), 80)
    nodewarn <- paste("\n",
                      "The parameter nodebase has been deprecated.",
                      "It is currently being kept for backwards compatibility,",
                      "but will be removed in a future version of poppr.",
                      "\n\nPlease use the new parameter nodescale instead.")
    nodewarn <- strwrap(nodewarn, width = warnw)
    warning(paste(nodewarn, collapse = "\n"))
    # Change the size of the vertices to a log scale.
    if (nodebase == 1) {
      nodebase <- 1.15
      nodewarn <- paste0("Cannot set nodebase = 1:\n",
        "Log base 1 is undefined, reverting to nodebase = 1.15")
      warning(nodewarn)
    }
    vsize <- log(V(poppr_msn$graph)$size * V(poppr_msn$graph)$size, base = nodebase) + 3
  } else {
    vsize <- sqrt(V(poppr_msn$graph)$size * V(poppr_msn$graph)$size * nodescale)
  }
  # Plotting parameters.
  def.par <- par(no.readonly = TRUE)
  # Return top level plot to defaults.
  on.exit({
    graphics::layout(matrix(1, ncol = 1, byrow = TRUE))
    graphics::par(def.par)
  })
  
  pop.leg <- pop.leg && !all(is.na(poppr_msn$populations))
  
  if (!pop.leg && !size.leg && !scale.leg) {
    plot.igraph(poppr_msn$graph, vertex.label = labs, vertex.size = vsize, 
                layout = lay, ...)
  } else {
    
    # Setting up the matrix for plotting. One Vertical panel of width 1 and height
    # 5 for the legend, one rectangular panel of width 4 and height 4.5 for the
    # graph, and one horizontal panel of width 4 and height 0.5 for the greyscale.
    if (pop.leg || size.leg) {
      if (scale.leg) {
        mat <- matrix(c(1,4,
                        3,2), 
                      ncol = 2, 
                      byrow = TRUE)
        layout(mat, widths = c(1, 4), heights = c(4.5, 0.5))        
      } else {
        layout(matrix(c(1, 2), ncol = 2), widths = c(1, 4))
      }

      # mar = bottom left top right
      
      ## LEGEND
      par(mar = c(0, 0, 1, 0) + 0.5)
      a <- NULL
      main <- if (pop.leg) 'POPULATION' else NULL
      too_many_pops   <- as.integer(ceiling(nPop(x)/30))
      pops_correction <- ifelse(too_many_pops > 1, -1, 1)
      yintersperse    <- ifelse(too_many_pops > 1, 0.51, 0.62)
      graphics::plot(c(-1, 1), 
                     c(-0.5, 0.5), 
                     type = 'n', 
                     axes = FALSE, 
                     xlab = '', 
                     ylab = '',
                     main = main)
      if (pop.leg) {
        
        a <- graphics::legend("topleft", 
                              bty = "n", 
                              cex = 1.2^pops_correction,
                              legend = poppr_msn$populations, 
                              fill = poppr_msn$color, 
                              border = NULL,
                              ncol = too_many_pops, 
                              x.intersp = 0.45, 
                              y.intersp = yintersperse)
      }
      if (size.leg) {
        N <- V(poppr_msn$graph)$size * V(poppr_msn$graph)$size
        make_circle_legend(a = a, 
                           mlg_number = N, 
                           scale = sqrt(nodescale), 
                           cex = 1.2 ^ pops_correction,
                           radmult = 4,
                           xspace = 0.45,
                           font = 2,
                           pos = 0,
                           x = 0,
                           y = 0.55,
                           txt = 0.05)
      }

    } else {
      graphics::layout(matrix(c(2, 1), nrow = 2), heights = c(4.5, 0.5))
    }
    if (scale.leg) {
      ## SCALE BAR
      if (quantiles) {
        scales <- sort(weights)
      } else {
        scales <- seq(wmin, wmax, l = 1000)
      }
      greyscales <- grDevices::gray(adjustcurve(scales, show=FALSE, glim=glim, correction=gadj))
      legend_image <- grDevices::as.raster(matrix(greyscales, nrow=1))
      graphics::par(mar = c(0, 1, 0, 1) + 0.5)
      graphics::plot.new()
      graphics::rasterImage(legend_image, 0, 0.5, 1, 1)
      graphics::polygon(c(0, 1, 1), c(0.5, 0.5, 0.8), col = "white", border = "white", lwd = 2)
      graphics::axis(3, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(scales), 3))
      graphics::text(0.5, 0, labels = "DISTANCE", font = 2, cex = 1.5, adj = c(0.5, 0))
    }
    ## PLOT
    if ((pop.leg || size.leg) && scale.leg) frame()
    par(mar = c(0,0,0,0))
    plot.igraph(poppr_msn$graph, 
                vertex.label = labs, 
                vertex.size = vsize, 
                layout = lay, 
                ...)
    
  }
  return(invisible(poppr_msn))
}


#==============================================================================#
#' Produce a genotype accumulation curve
#' 
#' Genotype accumulation curves are useful for determining the minimum number of
#' loci necessary to discriminate between individuals in a population. This 
#' function will randomly sample loci without replacement and count the number 
#' of multilocus genotypes observed.
#' 
#' @param gen a \code{\linkS4class{genclone}}, \code{\linkS4class{genind}}, or
#'   \code{\link[pegas]{loci}} object.
#'   
#' @param sample an \code{integer} defining the number of times loci will be 
#'   resampled without replacement.
#'   
#' @param maxloci the maximum number of loci to sample. By default, 
#'   \code{maxloci = 0}, which indicates that n - 1 loci are to be used. Note 
#'   that this will always take min(n - 1, maxloci)
#'   
#' @param quiet if \code{FALSE} (default), Progress of the iterations will be 
#'   displayed. If \code{TRUE}, nothing is printed to screen as the function 
#'   runs.
#'   
#' @param thresh a number from 0 to 1. This will draw a line at that fraction of
#'   multilocus genotypes, rounded. Defaults to 1, which will draw a line at the
#'   maximum number of observable genotypes.
#'   
#' @param plot if \code{TRUE} (default), the genotype curve will be plotted via 
#'   ggplot2. If \code{FALSE}, the resulting matrix will be visibly returned.
#'   
#' @param drop if \code{TRUE} (default), monomorphic loci will be removed before
#'   analysis as these loci affect the shape of the curve.
#'
#' @param dropna if \code{TRUE} (default) and \code{drop = TRUE}, NAs will be
#'   ignored when determining if a locus is monomorphic. When \code{FALSE},
#'   presence of NAs will result in the locus being retained. This argument has
#'   no effect when \code{drop = FALSE}
#'   
#' @return (invisibly by deafuls) a matrix of integers showing the results of
#'   each randomization. Columns represent the number of loci sampled and rows 
#'   represent an independent sample.
#'   
#' @details Internally, this function works by converting the data into a 
#'   \code{\link[pegas]{loci}} object, which represents genotypes as a data 
#'   frame of factors. Random samples are taken of 1 to n-1 columns of the 
#'   matrix and the number of unique rows are counted to determine the number of
#'   multilocus genotypes in that random sample. This function does not take 
#'   into account any definitions of MLGs via \code{\link{mlg.filter}} or 
#'   \code{\link{mll.custom}}.
#'   
#' @author Zhian N. Kamvar
#' @export
#' @examples
#' data(nancycats)
#' nan_geno <- genotype_curve(nancycats)
#' \dontrun{
#' 
#' # Marker Type Comparison --------------------------------------------------
#' # With AFLP data, it is often necessary to include more markers for resolution
#' data(Aeut)
#' Ageno <- genotype_curve(Aeut)
#' 
#' # Many microsatellite data sets have hypervariable markers
#' data(microbov)
#' mgeno <- geotype_curve(microbov)
#' 
#' # Adding a trendline ------------------------------------------------------
#' 
#' # Trendlines: you can add a smoothed trendline with geom_smooth()
#' library("ggplot2")
#' p <- last_plot()
#' p + geom_smooth()
#' 
#' # Producing Figures for Publication ---------------------------------------
#'
#' # This data set has been pre filtered
#' data(monpop)
#' mongeno <- genotype_curve(monpop)
#' 
#' # Here, we add a curve and a title for publication
#' p <- last_plot()
#' mytitle <- expression(paste("Genotype Accumulation Curve for ", 
#'                             italic("M. fructicola")))
#' p + geom_smooth() + 
#'   theme_bw() + 
#'   theme(text = element_text(size = 12, family = "serif")) +
#'   theme(title = element_text(size = 14)) +
#'   ggtitle(mytitle)
#' }
#' 
#==============================================================================#
#' @importFrom pegas loci2genind
genotype_curve <- function(gen, sample = 100, maxloci = 0L, quiet = FALSE, 
                           thresh = 1, plot = TRUE, drop = TRUE, dropna = TRUE){
  datacall <- match.call()
  if (!inherits(gen, c("genind", "genclone", "loci"))){
    stop(paste(datacall[2], "must be a genind or loci object"))
  }
  quiet <- should_poppr_be_quiet(quiet)
  if (inherits(gen, "loci")){
    genloc <- gen
    gen    <- pegas::loci2genind(gen)
  } else {
    genloc <- pegas::as.loci(gen)
  }
  if (!is.genclone(gen)) gen <- as.genclone(gen)
  if (nLoc(gen) == 1){
    stop("genotype_curve needs data with at least two loci.")
  }
  the_loci <- attr(genloc, "locicol")
  res      <- integer(nrow(genloc))
  suppressWarnings(genloc <- vapply(genloc[the_loci], as.integer, res))
  if (drop){
    nas <- if (dropna) !is.na(genloc) else TRUE
    polymorphic <- colSums(!apply(genloc, 2, duplicated) & nas) > 1
    lnames      <- locNames(gen)
    gen         <- gen[loc = polymorphic]
    genloc      <- genloc[, polymorphic, drop = FALSE]
    if (sum(!polymorphic) > 0){
      msg <- paste("Dropping monomorphic loci:", paste(lnames[!polymorphic], collapse = ", "))
      message(msg)
    }
  }
  nloci  <- min(maxloci, nLoc(gen) - 1)
  nloci  <- as.integer(ifelse(nloci > 0, nloci, nLoc(gen) - 1))
  sample <- as.integer(sample)
  report <- ifelse(quiet, 0L, as.integer(sample/100))
  out    <- .Call("genotype_curve_internal", 
                  mat     = genloc, 
                  iter    = sample, 
                  maxloci = nloci, 
                  report  = report, 
                  PACKAGE = "poppr")
  if (!quiet) cat("\n")
  colnames(out)        <- seq(nloci)
  rownames(out)        <- seq(sample)
  names(dimnames(out)) <- c("sample", "NumLoci")
  # Visibly return the data if plot = FALSE
  if (!plot) {
    return(out)
  }
  suppressWarnings(max_obs <- nmll(gen, "original"))
  threshdf                 <- data.frame(x = round(max_obs*thresh))

  outmelt <- as.data.frame.table(out, 
                                 responseName = "MLG",
                                 stringsAsFactors = FALSE)
  outmelt$sample  <- as.integer(outmelt$sample)
  outmelt$NumLoci <- as.integer(outmelt$NumLoci)
  aesthetics      <- aes_string(x = "NumLoci", y = "MLG")
  outplot <- ggplot(outmelt, aesthetics) + 
             geom_boxplot(aes_string(group = "factor(NumLoci)")) + 
             labs(list(title = paste("Genotype accumulation curve for", datacall[2]), 
                       y           = "Number of multilocus genotypes",
                       x           = "Number of loci sampled")) +
             scale_x_continuous(breaks = seq(nloci), expand = c(0, 0.125))
  if (!is.null(thresh)) {
    outbreaks <- sort(c(pretty(0:max_obs), threshdf$x))
    tjust     <- ifelse(thresh > 0.9, 1.5, -1)
    outplot   <- outplot + 
                 geom_hline(aes_string(yintercept = "x"), data = threshdf, 
                            color = "red", linetype = 2) + 
                 annotate("text", x = 1, y = threshdf$x, vjust = tjust, 
                          label = paste0(thresh*100, "%"), color = "red",
                          hjust = 0) +
                 scale_y_continuous(breaks = outbreaks, limits = c(0, max_obs))
  }
  print(outplot)
  return(invisible(out))
}

