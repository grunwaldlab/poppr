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
#' Filter a distance matrix and genind/genlight object at a given threshold
#'
#' @param gid a genind/genclone or genlight/snpclone object
#' @param threshold a numeric value specifying the threshold at which to filter
#' @param indist a distance matrix OR NULL. If NULL, the distance will be
#'   calculated with Bruvo's distance
#' @param clustering.algorithm the clustering algorithm to be passed to
#'   mlg.filter
#' @param bruvo_args a list of arguments for Bruvo's distance.
#'
#' @return a list containing the clone-corrected data (gid) and distance matrix
#'   (gid)
#' @noRd
filter_at_threshold <- function(gid, threshold, indist, clustering.algorithm,
                                bruvo_args){
  if (is.null(clustering.algorithm)){
    clustering.algorithm <- "farthest_neighbor"
  }
  if (is.null(indist)){
    filter.stats <- mlg.filter(gid, 
                               threshold, 
                               distance = bruvo.dist, 
                               algorithm = clustering.algorithm,
                               stats="ALL", 
                               replen = bruvo_args$replen,
                               add =  bruvo_args$add,
                               loss = bruvo_args$loss)    
  } else {
    filter.stats <- mlg.filter(gid, 
                               threshold, 
                               distance = indist, 
                               algorithm = clustering.algorithm,
                               stats="ALL") 
  }

  # TODO: The following two lines should be a product of mlg.filter
  visible(gid$mlg) <- "contracted"
  gid$mlg[]        <- filter.stats[["MLGS"]] 
  # Obtaining population information for all MLGs
  cgid   <- gid[.clonecorrector(gid), ]
  indist <- filter.stats[["DISTANCES"]]
  indist <- if (!is.matrix(indist)) as.matrix(indist) else indist
  # Fix issue #66
  rownames(indist) <- indNames(cgid) -> colnames(indist)
  return(list(gid = cgid, indist = indist))
}

#' Construct Minimum Spanning Networks
#'
#' @param gid a genind/genclone or genlight/snpclone
#' @param cgid the clone corrected version of that above
#' @param palette a character vector or palette function
#' @param indist a square distance matrix
#' @param include.ties logical. When TRUE, this will include exact ties in the
#'   minimum spanning network
#' @param mlg.compute character. Either "original" or "contracted". This is how
#'   to compute the MLGs when the input data is set to "custom"
#' @param vlab vertex labels. Either "INDS", "MLGS" or a vector of labels
#' @param visible_mlg This is the mlg status of the original data
#' @param wscale width scale
#' @param gscale grey scale
#' @param glim grey limit
#' @param gadj grey adjust
#' @param showplot logical whether to show the plot
#' @param ... params to be passed to igraph
#'
#' @return a list containing a graph, population names, and colors
#' @noRd
msn_constructor <-
  function(gid, cgid, palette, indist, include.ties, mlg.compute, vlab, 
           visible_mlg, wscale, gscale, glim, gadj, showplot, ...) {
    
  mlgs <- mll(gid)
  cmlg <- mll(cgid)
  if (!is.numeric(mlgs)){
    mlgs <- as.character(mlgs)
    cmlg <- as.character(cmlg)
  }
  
  
  # MSN colors --------------------------------------------------------------
  ## Color schemes 
  # The pallete is determined by what the user types in the argument. It can be 
  # rainbow, topo.colors, heat.colors ...etc.
  npop   <- nPop(gid)
  pnames <- popNames(gid)
  color  <- palette_parser(palette, npop, pnames)
  # This will determine the size of the nodes based on the number of
  # individuals in the MLG. Subsetting by the MLG vector of the clone
  # corrected set will give us the numbers and the population information in
  # the correct order. Note: rank is used to correctly subset the data
  mlg.number <- table(mlgs)[rank(cmlg)]
  # The MSN should not be drawn as a pie if there is a single population or
  # there is no population structure.
  piece_of_pie <- !is.null(pnames) && npop > 1
  if (piece_of_pie){
    # Obtaining population information for all MLGs
    subs   <- sort(unique(mlgs))
    mlg.cp <- mlg.crosspop(gid, mlgsub = subs, quiet=TRUE)
    if (is.numeric(mlgs)){
      names(mlg.cp) <- paste0("MLG.", sort(unique(mlgs)))    
    }
    mlg.cp <- mlg.cp[rank(cmlg)]  
    # This creates a list of colors corresponding to populations.
    mlg.color <- lapply(mlg.cp, function(x) color[pnames %in% names(x)])
  }


  # Creating the Minimum Spanning Network -----------------------------------
  g <- graph.adjacency(indist, weighted = TRUE, mode = "undirected")

  if (length(cmlg) > 1){
    mst <- minimum.spanning.tree(g, algorithm = "prim", weights = E(g)$weight)
    # Add any relevant edges that were cut from the mst while still being tied
    # for the title of optimal edge.
    if (include.ties){
      mst <- add_tied_edges(mst, indist, tolerance = .Machine$double.eps ^ 0.5)
    }
  } else { # if there's only one clone
    mst <- minimum.spanning.tree(g)
  }

  # Handling vertex labels --------------------------------------------------
  if (!is.na(vlab[1]) & length(vlab) == 1){
    if (toupper(vlab) == "MLG"){
       if (visible_mlg == "custom"){
        mll(gid) <- visible_mlg
        vlab     <- correlate_custom_mlgs(gid, mlg.compute)
      } else if (is.numeric(cmlg)){
        vlab <- paste0("MLG.", cmlg)        
      }
    } else if (toupper(vlab) == "INDS"){
      vlab <- indNames(cgid)
    }
  }

  # Adjusting edge color/thickness ------------------------------------------
  if (length(cmlg) > 1){
    mst <- update_edge_scales(mst, wscale, gscale, glim, gadj)
  }

  if (showplot){
    if (piece_of_pie){
      plot.igraph(
        mst,
        edge.width = E(mst)$width,
        edge.color = E(mst)$color,
        vertex.size = mlg.number * 3,
        vertex.shape = "pie",
        vertex.pie = mlg.cp,
        vertex.pie.color = mlg.color,
        vertex.label = vlab,
        ...
      )      
    } else {
      plot.igraph(
        mst,
        edge.width = E(mst)$width,
        edge.color = E(mst)$color,
        vertex.label = vlab,
        vertex.size = mlg.number * 3,
        vertex.color = color,
        ...
      )
    }
    graphics::legend(
      -1.55,
      1,
      bty = "n",
      cex = 0.75,
      legend = pnames,
      title = "Populations",
      fill = color,
      border = NULL
    )
  }
  
  V(mst)$size  <- mlg.number
  if (piece_of_pie){
    V(mst)$shape     <- "pie"
    V(mst)$pie       <- mlg.cp
    V(mst)$pie.color <- mlg.color    
  } else {
    V(mst)$color     <- color
  }
  V(mst)$label <- vlab
  return(list(graph = mst, populations = pnames, colors = color))
}
