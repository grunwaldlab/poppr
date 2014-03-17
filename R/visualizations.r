
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
#
# This will make a histogram of permutations of a dataset and show a red line
# indicating p = 0.05 and a green line indicating the observed value. 
# 
#==============================================================================#
permut.histogram.maker <- function(index="index",sampled,observed,pop){
  if(all(is.nan(sampled[[index]]))){
    warning("Cannot create histogram due to NaN's")
    next
  }
	xmin <- min(observed, sampled[[index]])
	xmax <- max(observed, sampled[[index]])
  # R does not like it when the max value for xlim is infinity, so this will fix
  # the culprit where it occurs, usually in the rbarD caluclation. I will be
  # surprised if I see it in the Ia calculation.
  if(!is.nan(xmax) & xmax==Inf & index=="rbarD"){
    xmax <- 1
  }
  if(is.nan(xmin)){
    xmin <- -0.5
  }
  if(is.nan(xmax)){
    xmax <- ifelse(index=="Ia", 18, 1)
  }
	hist(sampled[[index]], xlim=c(xmin, xmax),
        main=c("Population:",as.character(pop)),
        xlab=sprintf("%s (%d permutations)",names(sampled[index]),
                      length(sampled[[index]])), col="grey")
	abline(v=observed, col="green")
  perc95 <- quantile(sampled[[index]], 0.95, na.rm=TRUE)[[1]]
  abline(v=perc95, col="red")
}
#==============================================================================#
# histogram will create two histograms showing the sampled distributions of the
# index of association and the standardized index of association along with an
# informative legend between them.
#==============================================================================#
permut.histogram <- function(sampled, observed, pval, pop="pop", file="file"){
  par(xpd=TRUE,mfrow=c(3,1))
  permut.histogram.maker(index="Ia",sampled,observed[1],pop)
  if (!is.na(pval)){
    pval <- ifelse(pval==0, sprintf("< %g", 1/length(sampled[["Ia"]])), pval)
    frame()
    legend('center', col=c("red", "green"), lty=c(1,1), c("p = 0.05",
              paste("Observed\n","(p-value  ",pval,")", sep="")), 
                        horiz=TRUE, title=paste("File : ",file,sep=""))
  }
  permut.histogram.maker(index="rbarD",sampled,observed[2],pop)
}

poppr.plot <- function(sample, pval = c("0.05", "0.05"), pop="pop", 
                        observed = observed, file="file", N=NA){
  #   if (!is.na(pval[1])){
  #     pval[1] <- ifelse(pval[1]==0, sprintf("< %g", 1/length(sample[["Ia"]])), pval[1])
  #   }
  #   if (!is.na(pval[2])){
  #     pval[2] <- ifelse(pval[2]==0, sprintf("< %g", 1/length(sample[["rbarD"]])), pval[2])
  #   }
  #````````````````````````````````````````````````````````````````````````````#
  # In the case that the sample contains all NaNs, a graph cannot be displayed.
  # In place of the graph will be an orange box with a warning message.
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
  
  if (all(is.nan(sample$rbarD))){
    warning(paste("The data from ",file,", population: ",pop," contains only missing values and cannot be displayed graphically", sep=""))
    background_stuff <- theme(panel.grid.major.y = element_line(size=0)) +
      theme(panel.grid.minor.y = element_line(size=0)) +
      theme(panel.grid.major.x = element_line(size=0)) +
      theme(panel.background = element_rect(fill="grey95")) +
      #theme(plot.background = element_rect(size=1, color="black")) +
      theme(axis.ticks.y = element_line(size=0)) +
      theme(axis.text.y = element_text(size=0)) +
      theme(axis.ticks.x = element_line(size=0)) +
      theme(axis.text.x = element_text(size=0)) +
      theme(axis.title.y = element_text(size=rel(0)))+
      theme(axis.title.x = element_text(size=rel(0)))
    oops <- ggplot(as.data.frame(list(x=-10:9)), aes_string(x = "x")) + 
      geom_histogram(binwidth=1, fill="orange") + 
      geom_text(aes(label="Warning:", x=0, y=0.8), color="black", size=rel(15)) + 
      geom_text(aes(label="Data contains only NaNs and\ncannot be displayed graphically", 
                    x=0, y=0.5, hjust=0.5), color="black", size=rel(10)) +
      labs(title=paste("Population: ", pop, "; N: ", N, "\nPermutations: ", 
                       length(sample$Ia), "\nFile: ", file, sep="")) + 
      theme(plot.title = element_text(vjust=1, size=rel(2), face="bold")) +
      background_stuff
    print(oops)
  }
  else {
    #``````````````````````````````````````````````````````````````````````````#
    # Normal Cases
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
    if(any(is.nan(sample[["rbarD"]])))
      sample[["rbarD"]][which(is.nan(sample[["rbarD"]]))] <- mean(sample[["rbarD"]], na.rm=TRUE)
    if(any(is.nan(sample[["Ia"]])))
      sample[["Ia"]][which(is.nan(sample[["Ia"]]))] <- mean(sample[["Ia"]], na.rm=TRUE)
    
    # Transforming the data into a data frame that ggplot can understand and put
    # into facets. It will have the following rows: Value, Index, and Quant.
    Indexfac <- factor(1:2, levels=1:2, labels=c("I[A]","bar(r)[d]"))
    infodata <- as.data.frame(list(Value=c(sample[["Ia"]], sample[["rbarD"]]), 
                                   Index=rep(Indexfac, each=length(sample[["Ia"]]))))
    
    # Setting up the observed values for formatting the overlaying lines.
    # It will have the following rows: Observed, Index, Quant, Sam, P, min, max
    
    obsdata <- data.frame(list(Observed=observed[1:2], Index=Indexfac))
    obsdata[["P"]] <- paste("Observed \n(p-value: ", pval,")",sep="")
    obsdata[["median"]] <- c(median(infodata[infodata[["Index"]] == Indexfac[1], "Value"]),
                        median(infodata[infodata[["Index"]] == Indexfac[2], "Value"]))


    obsdata[["label"]] <- paste("Observed: ",obsdata[["Observed"]], sep="")
    if(any(is.na(observed))){
      warning(paste("The Index of Association values from ",file,", population: ",pop," contain missing values and cannot be displayed graphically", sep=""))
      iard_plot <- ggplot(infodata, aes_string(x = "Value")) + 
        # Giving the data over to the histogram creating function and removing
        # all of the lines from each bar, so it's displayed as a solid area.
        geom_histogram(linetype="blank", #alpha=0.8, 
                       data = infodata[infodata[["Index"]] == Indexfac[1], ], 
                       position="identity",
                       binwidth=diff(range(infodata[infodata[["Index"]] == Indexfac[1], "Value"]))/30) + 
        geom_histogram(linetype="blank", #alpha=0.8, 
                       data = infodata[infodata[["Index"]] == Indexfac[2], ], 
                       position="identity",
                       binwidth=diff(range(infodata[infodata[["Index"]] == Indexfac[2], "Value"]))/30) + 
        geom_rug(alpha = 0.5) + #theme(legend.position = "none") + 
        # The label for the observed line is a bit more difficult to code as
        # it has the ability to appear anywhere on the chart. Here, I'm
        # forcing it to flip to one side or the other based on which side of
        # the mean that the observed value falls on.
        geom_text(aes_string(label="label", x=0, y=Inf, vjust=1.5),
                  data=obsdata, angle=0, color="red") + 
        
        # Splitting the data into separate plots with free x-axes.
        facet_grid(". ~ Index", scales="free_x", labeller=label_parsed) +
        
        # Title of the plot.
        labs(title=paste("Population: ", pop, "; N: ", N, "\nPermutations: ", 
                         length(sample$Ia), "\nFile: ", file, sep="")) + 
        theme_classic() %+replace%
        theme(plot.title = element_text(vjust=1, size=rel(2), face="bold")) +
        theme(panel.background = element_rect(fill="grey98")) +
        
        # Making the Index titles bigger. 
        theme(strip.text.x = element_text(size=rel(3), face="bold"))
    }
    else{
      iard_plot <- ggplot(infodata, aes_string(x = "Value")) + 
        # Giving the data over to the histogram creating function and removing
        # all of the lines from each bar, so it's displayed as a solid area.
        geom_histogram(linetype="blank", #alpha=0.8, 
                       data = infodata[infodata[["Index"]] == Indexfac[1], ], 
                       position="identity",
                       binwidth=diff(range(infodata[infodata[["Index"]] == Indexfac[1], "Value"]))/30) + 
        geom_histogram(linetype="blank", #alpha=0.8, 
                       data = infodata[infodata[["Index"]] == Indexfac[2], ], 
                       position="identity",
                       binwidth=diff(range(infodata[infodata[["Index"]] == Indexfac[2], "Value"]))/30) + 
        geom_rug(alpha = 0.5) + #theme(legend.position = "none") + 
        # Positioning the observed line and labeling it.
        geom_vline(aes_string(xintercept = "Observed"), data = obsdata, color="blue", 
                   show_guide=TRUE, linetype="dashed") +
        
        # The label for the observed line is a bit more difficult to code as
        # it has the ability to appear anywhere on the chart. Here, I'm
        # forcing it to flip to one side or the other based on which side of
        # the mean that the observed value falls on.
        #geom_text(aes_string(label=paste("Observed \n(p-value: ", P,")",sep=""),
        geom_text(aes_string(label="P", 
                      x="Observed",y=Inf,vjust=2,hjust="ifelse(Observed > median, 1.01, -0.01)"),
                  data=obsdata, angle=0, color="blue") + 
        
        # Splitting the data into separate plots with free x-axes.
        facet_grid(". ~ Index", scales="free_x", labeller=label_parsed) +
        
        # Title of the plot.
        labs(title=paste("Population: ", pop, "; N: ", N, "\nPermutations: ", 
                         length(sample$Ia), "\nFile: ", file, sep="")) + 
        theme_classic() %+replace%
        theme(plot.title = element_text(vjust=1, size=rel(2), face="bold")) +
        theme(panel.background = element_rect(fill="grey98")) +
        
        # Making the Index titles bigger. 
        theme(strip.text.x = element_text(size=rel(3), face="bold"))
    }
    print(iard_plot)
  }
} 

#==============================================================================#
#
#' Create a minimum spanning network of selected populations using a distance
#' matrix.
#'
#' @param pop a \code{\link{genind}} object
#'
#' @param distmat a distance matrix that has been derived from your data set.
#' 
#' @param palette a \code{function} defining the color palette to be used to
#' color the populations on the graph. It defaults to \code{\link{topo.colors}},
#' but you can easily create new schemes by using \code{\link{colorRampPalette}}
#' (see examples for details)
#' 
#' @param sublist a \code{vector} of population names or indexes that the user
#' wishes to keep. Default to "ALL".
#'
#' @param blacklist a \code{vector} of population names or indexes that the user
#' wishes to discard. Default to \code{NULL}
#' 
#' @param vertex.label a \code{vector} of characters to label each vertex. There
#' are two defaults: \code{"MLG"} will label the nodes with the multilocus genotype
#' from the original data set and \code{"inds"} will label the nodes with the 
#' representative individual names.
#' 
#' @param gscale "grey scale". If this is \code{TRUE}, this will scale the color
#' of the edges proportional to the observed distance, with the lines becoming 
#' darker for more related nodes. See \code{\link{greycurve}} for details.
#' 
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#' the upper and lower limits for the \code{\link{gray}} function. Default is 0
#' (black) and 0.8 (20\% black). See \code{\link{greycurve}} for details.
#' 
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that
#' will serve as the exponent to the edge weight to scale the grey value to
#' represent that weight. See \code{\link{greycurve}} for details.
#' 
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale
#' will be weighted to emphasize the differences between closely related nodes.
#' If it is 2, the grey scale will be weighted to emphasize the differences
#' between more distantly related nodes. See \code{\link{greycurve}} for details.
#' 
#' @param wscale "width scale". If this is \code{TRUE}, the edge widths will be
#' scaled proportional to the inverse of the observed distance , with the lines 
#' becoming thicker for more related nodes.
#' 
#' @param ... any other arguments that could go into plot.igraph
#'
#' @return 
#' \item{graph}{a minimum spanning network with nodes corresponding to MLGs within
#' the data set. Colors of the nodes represent population membership. Width and
#' color of the edges represent distance.}
#' \item{populations}{a vector of the population names corresponding to the 
#' vertex colors}
#' \item{colors}{a vector of the hexadecimal representations of the colors used
#' in the vertex colors}
#' 
#' @note The edges of these graphs may cross each other if the graph becomes too
#' large. 
#'
#' @seealso \code{\link{nancycats}}, \code{\link{upgma}}, \code{\link{nj}},
#' \code{\link{nodelabels}}, \code{\link{na.replace}}, \code{\link{missingno}}, 
#' \code{\link{bruvo.msn}}, \code{\link{greycurve}}.
#' 
#' @export
#' @aliases msn.poppr
#' @author Javier F. Tabima, Zhian N. Kamvar
#' @examples
#' 
#' # Load the data set and calculate the distance matrix for all individuals.
#' data(Aeut)
#' A.dist <- diss.dist(Aeut)
#' 
#' # Graph it.
#' A.msn <- poppr.msn(Aeut, A.dist, gadj=15, vertex.label=NA)
#' 
#' \dontrun{
#' # Set subpopulation structure.
#' Aeut.sub <- splitcombine(Aeut, method=2, hier=c("Pop", "Subpop"))
#' 
#' # Plot respective to the subpopulation structure
#' As.msn <- poppr.msn(Aeut.sub, A.dist, gadj=15, vertex.label=NA)
#' 
#' # Show only the structure of the Athena population.
#' As.msn <- poppr.msn(Aeut.sub, A.dist, gadj=15, vertex.label=NA, sublist=1:10)
#' 
#' # Let's look at the structure of the microbov data set
#' data(microbov)
#' micro.dist <- diss.dist(microbov)
#' micro.msn <- poppr.msn(microbov, diss.dist(microbov), vertex.label=NA)
#' 
#' Let's plot it and show where individuals have < 15% of their genotypes 
#' different.
#' 
#' plot.igraph(micro.msn$graph, edge.label = ifelse(E(micro.msn$graph)$weight < 0.15, 
#' round(E(micro.msn$graph)$weight, 3), NA), vertex.size=2, edge.label.color="red")
#' 
#' }
#' 
#==============================================================================#
poppr.msn <- function (pop, distmat, palette = topo.colors, 
                       sublist = "All", blacklist = NULL, vertex.label = "MLG", 
                       gscale=TRUE, glim = c(0,0.8), gadj = 3, gweight = 1, 
                       wscale=TRUE, ...){
  if (class(distmat) != "dist"){
    if (is.matrix(distmat)){
      if (any(nInd(pop) != dim(distmat))){
        stop("The size of the distance matrix does not match the size of the data.\n")
      }
      distmat <- as.dist(distmat)
    } else {
      stop("The distance matrix is neither a dist object nor a matrix.\n")
    }
  }
  if (nInd(pop) != attr(distmat, "Size")){
    stop("The size of the distance matrix does not match the size of the data.\n")
  }
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  # Storing the MLG vector into the genind object
  if (!is.genclone(pop)){
    # Storing the MLG vector into the genind object
    pop$other$mlg.vec <- mlg.vector(pop)  
  }
  bclone <- as.matrix(distmat)

  # The clone correction of the matrix needs to be done at this step if there
  # is only one or no populations. 
  if (is.null(pop(pop)) | length(pop@pop.names) == 1){
    if (is.genclone(pop)){
      mlgs <- pop@mlg
    } else {
      mlgs <- pop$other$mlg.vec
    }
    bclone <- bclone[!duplicated(mlgs), !duplicated(mlgs)]
    return(singlepop_msn(pop, vertex.label, distmat = bclone, gscale = gscale, 
                         glim = glim, gadj = gadj, wscale = wscale, 
                         palette = palette))
  }
  # This will subset both the population and the matrix. 
  if(sublist[1] != "ALL" | !is.null(blacklist)){
    sublist_blacklist <- sub_index(pop, sublist, blacklist)
    bclone <- bclone[sublist_blacklist, sublist_blacklist]
    pop <- popsub(pop, sublist, blacklist)
  }
  cpop <- pop[.clonecorrector(pop), ]
  if (is.genclone(pop)){
    mlgs <- pop@mlg
    cmlg <- cpop@mlg
  } else {
    mlgs <- pop$other$mlg.vec
    cmlg <- cpop$other$mlg.vec
  }
  # This will clone correct the incoming matrix. 
  bclone <- bclone[!duplicated(mlgs), !duplicated(mlgs)]
  
  if (is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop_msn(pop, vertex.label, distmat = bclone, gscale = gscale, 
                         glim = glim, gadj = gadj, wscale = wscale, 
                         palette = palette))
  }
  # Obtaining population information for all MLGs
  mlg.cp <- mlg.crosspop(pop, mlgsub=1:mlg(pop, quiet=TRUE), quiet=TRUE)

  names(mlg.cp) <- paste0("MLG.", sort(unique(mlgs)))
  # This will determine the size of the nodes based on the number of individuals
  # in the MLG. Subsetting by the MLG vector of the clone corrected set will
  # give us the numbers and the population information in the correct order.
  # Note: rank is used to correctly subset the data
  mlg.number <- table(mlgs)[rank(cmlg)]
  mlg.cp     <- mlg.cp[rank(cmlg)]
  rownames(bclone) <- cpop$pop
  colnames(bclone) <- cpop$pop
  
  g   <- graph.adjacency(bclone, weighted=TRUE, mode="undirected")
  mst <- minimum.spanning.tree(g, algorithm="prim", weights=E(g)$weight)
  
  if (!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if(toupper(vertex.label) == "MLG"){
      vertex.label <- paste("MLG.", cmlg, sep="")
    }
    else if(toupper(vertex.label) == "INDS"){
      vertex.label <- cpop$ind.names
    }
  }
  ###### Color schemes #######  
  # The pallete is determined by what the user types in the argument. It can be 
  # rainbow, topo.colors, heat.colors ...etc.
  palette <- match.fun(palette)
  color   <- palette(length(pop@pop.names))
  
  ###### Edge adjustments ######
  # Grey Scale Adjustment weighting towards more diverse or similar populations.
  if (gscale == TRUE){
    E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correction=gadj, 
                                     show=FALSE))
  } else {
    E(mst)$color <- rep("black", length(E(mst)$weight))
  }
  
  # Width scale adjustment to avoid extremely large widths.
  # by adding 0.08 to entries, the max width is 12.5 and the min is 0.9259259
  edgewidth <- 2
  if (wscale==TRUE){
    edgewidth <- 1/(E(mst)$weight)
    if (any(E(mst)$weight < 0.08)){
      edgewidth <- 1/(E(mst)$weight + 0.08)
    }
  }
  
  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pop@pop.names %in% names(x)])
  
  plot.igraph(mst, edge.width = edgewidth, edge.color = E(mst)$color, 
       vertex.size = mlg.number*3, vertex.shape = "pie", vertex.pie = mlg.cp, 
       vertex.pie.color = mlg.color, vertex.label = vertex.label, ...)
  legend(-1.55 ,1 ,bty = "n", cex = 0.75, legend = pop$pop.names, 
         title = "Populations", fill=color, border=NULL)
  E(mst)$width     <- edgewidth
  V(mst)$size      <- mlg.number
  V(mst)$shape     <- "pie"
  V(mst)$pie       <- mlg.cp
  V(mst)$pie.color <- mlg.color
  V(mst)$label     <- vertex.label
  return(list(graph = mst, populations = pop$pop.names, colors = color))
}


#==============================================================================#
#' Create a table summarizing missing data a genind or genclone object
#' 
#' @param x a \linkS4class{genind} or \linkS4class{genclone} object.
#' 
#' @param percent \code{logical}. If \code{TRUE} (default), table and plot will
#'   represent missing data as a percentage of each cell. If \code{FALSE}, the 
#'   table and plot will represent missing data as raw counts. (See details)
#' 
#' @param plot \code{logical}. If \code{TRUE}, a simple heatmap will be
#'   produced. If \code{FALSE} (default), no heatmap will be produced.
#' 
#' @param df \code{logical}. If \code{TRUE}, the data will be returned as a long
#'   form data frame. If \code{FALSE} (default), a matrix with populations in
#'   rows and loci in columns will be returned.
#' 
#' @param returnplot \code{logical}. If \code{TRUE}, a list is returned with two
#'   elements: \code{table} - the normal output and \code{plot} - the ggplot
#'   object. If \code{FALSE}, the missing table is returned.
#' 
#' @param low \code{character}. What color should represent no missing data?
#'   (default: "blue")
#' 
#' @param high \code{character}. What color should represent the highest amount
#'   of missing data? (default: "red")
#'   
#' @param plotlab \code{logical}. If \code{TRUE} (default), values of missing
#'   data greater than 0\% will be plotted. If \code{FALSE}, the plot will
#'   appear unappended.
#' 
#' @param scaled \code{logical}. This is for when \code{percent = TRUE}. If
#'   \code{TRUE} (default), the color specified in \code{high} will represent
#'   the highest observed value of missing data. If \code{FALSE}, the color
#'   specified in \code{high} will represent 100\%.
#'   
#' @return a matrix, data frame (\code{df = TRUE}), or a list (\code{returnplot = TRUE}) representing missing data in a \linkS4class{genind} or \linkS4class{genclone} object.
#' 
#' @details This data is potentially useful for identifying areas of systematic
#' missing data. It is important to note that when visualizing missing data as a
#' percentage, that the percentage of missing data is \strong{relative to the
#' population and locus}, not to the entire data set. Likewise, the last column
#' represents an average amount of missing data over all loci for each
#' population, and the last row, "Total", represents the entire data set.
#' 
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' data(nancycats)
#' nancy.miss <- missing_table(nancycats, plot = TRUE)
#' 
#==============================================================================#
missing_table <- function(x, percent = TRUE, plot = FALSE, df = FALSE, 
                          returnplot = FALSE, low = "blue", high = "red", 
                          plotlab = TRUE, scaled = TRUE){
  pops       <- seppop(x, drop = FALSE)
  pops$Total <- x
  inds       <- 1
  if (percent){
    inds <- c(table(pop(x)), nInd(x))
  }
  misstab              <- matrix(0, nrow = nLoc(x) + 1, ncol = length(pops))
  misstab[1:nLoc(x), ] <- vapply(pops, number_missing_locus, numeric(nLoc(x)), 1)
  misstab[-nrow(misstab), ] <- t(apply(misstab[-nrow(misstab), ], 1, "/", inds))
  misstab[nrow(misstab), ]  <- colMeans(misstab[-nrow(misstab), ])
  rownames(misstab)         <- c(x@loc.names, "Mean")
  colnames(misstab)         <- names(pops)
  if (all(misstab == 0)){
    cat("No Missing Data Found!")
    return(NULL)
  }
  if (plot){
    missdf      <- melt(misstab, varnames = c("Locus", "Population"), 
                        value.name = "Missing")
    leg_title   <- "Missing"
    missdf[1:2] <- data.frame(lapply(missdf[1:2], 
                                     function(x) factor(x, levels = unique(x))))
    # levels(missdf[[1]]) <- rev(levels(missdf[[1]]))
    plotdf <- textdf <- missdf
    if (percent) {
      plotdf$Missing <- round(plotdf$Missing*100, 2)
      textdf$Missing <- paste(plotdf$Missing, "%")
      miss           <- "0 %"
      title          <- paste("Percent missing data per locus and population of", 
                              as.character(substitute(x)))
      leg_title      <- paste("Percent", leg_title)
      if(!scaled){
        lims <- c(0, 100)
      }
    } else {
      textdf$Missing <- round(textdf$Missing, 2)
      miss           <- 0
      title          <- paste("Missing data per locus and population of", 
                              as.character(substitute(x)))
    }
    if (scaled | !percent){
      lims <- c(0, max(plotdf$Missing))
    }
    linedata       <- data.frame(list(yint = 1.5, #ncol(misstab) - 0.5, 
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
      scale_y_discrete(expand = c(0, -1), 
                       limits = rev(unique(plotdf$Population))) + 
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
    class(misstab) <- c("locustable", "matrix")
  }
  if (returnplot){
    misstab <- list(table = misstab, plot = outplot)
  }
  return(misstab)
}

#==============================================================================#
#' Display a greyscale gradient adjusted to specific parameters
#' 
#' This function has one purpose. It is for deciding the appropriate scaling for
#' a grey palette to be used for edge weights of a minimum spanning network.
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
greycurve <- function(glim = c(0,0.8), gadj = 3, gweight = 1){
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  adjustcurve(seq(0.001, 1, 0.001), glim, correction=gadj, show=TRUE)
}
