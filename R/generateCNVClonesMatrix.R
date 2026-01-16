#' Generate CNV Matrix for CNV Clusters by Chromosome Arm
#'
#' This function generates a matrix of metacells where each metacell corresponds to a CNV cluster.
#' The CNV matrix is calculated by chromosome arm. If specified, certain clusters will be labeled as "Benign"
#' rather than "Clone".
#'
#' @param seuratObj A Seurat object containing CNV data and metadata.
#' @param healthyClusters A numeric vector or `NULL`. If provided, clusters specified in this vector
#' will be labeled as "Benign" instead of "Clone". Default is `NULL`.
#' @param values one of 'scores' or 'calls'. 'scores' returns the mean CNV score per cluster,
#' while 'calls' uses `cnv_thresh` to establish a cut-off for gains and losses, returning a matrix
#' of CNV calls (0=none, 1=gain, -1=loss).
#' @param cnv_thresh A numeric threshold to filter significant CNV events. Default is 0.15.
#'
#' @return A matrix of CNVs with row names corresponding to the clone or benign labels and columns representing
#' the chromosome arms, with values corresponding to CNV scores or CNV calls.
#'
#' @export

generateCNVClonesMatrix <- function(seuratObj,  healthyClusters = NULL, values = "scores", cnv_thresh = 0.15) {
  # Extract CNV matrix based on chromosome arms
  cnv_matrix <- as.matrix(seuratObj[[which(names(seuratObj@meta.data) == "1.p_CNV") : which(names(seuratObj@meta.data) == "X.q_CNV")]])
  cnv_matrix <- cnv_matrix[Seurat::Cells(seuratObj),]

  # Initialize the cnv_matrix_clusters
  cnv_matrix_clusters <- matrix(nrow = 0, ncol = ncol(cnv_matrix))

  # Get unique clusters, excluding NA
  unique_clusters <- unique(seuratObj[["cnv_clusters"]][[1]])
  unique_clusters <- unique_clusters[!is.na(unique_clusters)]

  # Loop through valid clusters only
  for (cluster in unique_clusters) {
    cells <- rownames(seuratObj@meta.data)[which(seuratObj[["cnv_clusters"]] == cluster)]
    cnv_matrix_clusters <- rbind(cnv_matrix_clusters, colMeans(cnv_matrix[cells, , drop = FALSE]))
  }

  rownames(cnv_matrix_clusters) <- unique_clusters

  # Label clusters as "Clone X" or "Benign X"
  rownames(cnv_matrix_clusters) <- paste0("Clone ", rownames(cnv_matrix_clusters))
  if (!is.null(healthyClusters)) {
    for (hc in healthyClusters) {
      rownames(cnv_matrix_clusters)[rownames(cnv_matrix_clusters) == paste0("Clone ", as.character(hc))] <- paste0("Benign ", hc)
    }
  }

  if(values == "scores") {
    # Return the CNV matrix
    return(cnv_matrix_clusters)
  } else if (values == "calls") {
    alt_matrix <- matrix(0, nrow = nrow(cnv_matrix_clusters), ncol = ncol(cnv_matrix_clusters),
                         dimnames = dimnames(cnv_matrix_clusters))
    alt_matrix[cnv_matrix_clusters >= cnv_thresh] <- 1
    alt_matrix[cnv_matrix_clusters <= -cnv_thresh] <- -1
    return(alt_matrix)
  } else {
    stop("Non supported `values` value. Must be one of: 'scores', 'calls'")
  }

}
