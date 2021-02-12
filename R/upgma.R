#' UPGMA
#'
#' UPGMA clustering. Just a wrapper function around \code{\link[stats]{hclust}}.
#'
#' @param d A distance matrix.
#' @return A phylogenetic tree of class \code{phylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{hclust}}, \code{\link{as.phylo}}
#' @importFrom ape as.phylo.hclust as.phylo
#' @importFrom stats hclust as.dist
#' @keywords cluster
#' @examples
#'
#' library(ape)
#' data(woodmouse)
#' dm <- dist.dna(woodmouse)
#' tree <- upgma(dm)
#' plot(tree)
#'
#' @rdname upgma
#' @export
"upgma" <- function(d) as.phylo(hclust(as.dist(d), method = "average"))
