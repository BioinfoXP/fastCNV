#' Construct a Phylogenetic Tree from a Copy Number Variation (CNV) Matrix
#'
#' This function constructs a phylogenetic tree based on a given copy number variation (CNV) matrix.
#' It adds a baseline "Normal" profile only to root the tree, which is not shown in final output.
#' First, it computes pairwise distances between profiles using Euclidean distance, and then applies a specified tree-building function
#' (e.g., Neighbor-Joining) to construct the tree.
#'
#' @param cnv_matrix A matrix representing copy number variation, where rows correspond to samples
#' and columns correspond to genomic regions. Each value represents the CNV at a given region in a sample.
#' @param tree_function A function to construct the phylogenetic tree from a distance matrix. The default is
#' `nj` (Neighbor-Joining). Other functions (e.g., `upgma`, `wpgma`) can also be used.
#' @param dist_method The distance method to be used.
#'
#' @return A rooted phylogenetic tree (of class \code{phylo})
#'
#' @examples
#' # Example usage with Neighbor-Joining (default)
# cnv_matrix <- matrix(c(2, 3, 2, 4, 3, 5), nrow = 2, byrow = TRUE,
#                      dimnames = list(c("Clone1", "Clone2"), c("Region1", "Region2", "Region3")))
# tree <- buildCNVTree(cnv_matrix)
# plot(tree)
#'
#' @importFrom ape nj root drop.tip
#' @importFrom phangorn wpgma upgma Ancestors
#'
#' @export


buildCNVTree <- function(cnv_matrix, tree_function = nj, dist_method = "euclidean") {
  normal_root <- matrix(rep(0, ncol(cnv_matrix)), nrow = 1,
                        dimnames = list("Normal", colnames(cnv_matrix)))
  cnv_mat <- rbind(cnv_matrix, normal_root)

  cnv_dist <- dist(cnv_mat, method = dist_method)
  tree <- tree_function(cnv_dist)
  tree <- root(tree, "Normal", resolve.root = TRUE)
  tree <- drop.tip(tree, "Normal")
  return(tree)
}
