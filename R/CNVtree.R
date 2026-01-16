#' Build, annotate and plot a Phylogenetic Tree from a seurat Object containing the CNV results from `fastCNV()`
#'
#' @param seuratObj A Seurat object containing CNV data and metadata.
#' @param healthyClusters A numeric vector or `NULL`. If provided, clusters specified in this vector
#' will be labeled as "Benign" instead of "Clone". Default is `NULL`.
#' @param values one of 'scores' or 'calls'. 'scores' returns the mean CNV score per cluster,
#' while 'calls' uses `cnv_thresh` to establish a cut-off for gains and losses, returning a matrix
#' of CNV calls (0=none, 1=gain, -1=loss).
#' @param cnv_thresh A numeric threshold to filter significant CNV events. Default is 0.15.
#' @param tree_function A function to construct the phylogenetic tree from a distance matrix. The default is
#' `nj` (Neighbor-Joining). Other functions (e.g., `upgma`, `wpgma`) can also be used.
#' @param dist_method The distance method to be used.
#' @param clone_cols a color palette to color the clones. If NULL, points are
#' not colored. If TRUE, clones are colored using default color palette. If a
#' palette is given, clones are colored following the palette, with
#' values passed to `scale_color_manual`.
#'
#' @export


CNVTree <- function(seuratObj,
                    healthyClusters = NULL,
                    values = "scores",
                    cnv_thresh = 0.15,
                    tree_function = nj,
                    dist_method = "euclidean",
                    clone_cols = TRUE){

  #Get a unique CNV score per cluster and per chromosomic arm
  cnv_matrix_clusters <- generateCNVClonesMatrix(seuratObj = seuratObj,
                                                 healthyClusters = healthyClusters,
                                                 values = values,
                                                 cnv_thresh = cnv_thresh)

  #Build the CNV tree
  tree <- buildCNVTree(cnv_matrix_clusters,
                       tree_function = tree_function,
                       dist_method = dist_method)

  #Annotate the CNV tree
  tree_data <- annotateCNVTree(tree,
                               cnv_matrix_clusters,
                               cnv_thresh = cnv_thresh)

  #Plot the CNV tree
  plot(plotCNVTree(tree_data,
              clone_cols = clone_cols))

  #Return the annotated CNV tree
  return(tree_data)

}
