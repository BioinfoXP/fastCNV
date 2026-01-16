#' Merge CNV Clusters in a Seurat Object
#'
#' This function merges CNV clusters in a Seurat object based on the correlation
#' of their average CNV profiles across chromosome arms. Clusters with correlation
#' greater than a user-specified threshold are merged into a single cluster.
#'
#' @param seuratObj A Seurat object containing fastCNV's results by chromosome arm, and CNV clustering.
#' @param mergeThreshold A numeric value between 0 and 1. Clusters with correlation
#'   greater than this threshold will be merged. Default is 0.98.
#'
#' @return A Seurat Object with updated CNV clusters, where highly correlated clusters have been merged.
#'
#' @export

mergeCNVClusters <- function(seuratObj, mergeThreshold = 0.98){

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

  cnv_matrix_clusters_clean <- cnv_matrix_clusters[,colSums(cnv_matrix_clusters != 0) > 0 ]
  corrmat <- cor(t(cnv_matrix_clusters_clean))

  cluster_names <- rownames(corrmat)
  n <- nrow(corrmat)
  adj <- (corrmat > mergeThreshold) * 1
  groups <- 1:n

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (adj[i, j] == 1) {
        gmin <- min(groups[i], groups[j])
        gmax <- max(groups[i], groups[j])
        groups[groups == gmax] <- gmin
      }
    }
  }

  merged <- split(cluster_names, groups)

  mapping <- unlist(lapply(merged, function(grp) {
    setNames(rep(paste(grp, collapse = "_"), length(grp)), as.character(grp))
  }))

  names(mapping) <- sub(".*\\.", "", names(mapping))

  orig <- as.character(seuratObj@meta.data[["cnv_clusters"]])

  mapped <- mapping[orig]

  mapped[is.na(mapped)] <- orig[is.na(mapped)]
  names(mapped) <- rownames(seuratObj@meta.data)

  seuratObj@meta.data[["cnv_clusters"]] <- sub("_.*", "", mapped)
  seuratObj@meta.data[["cnv_clusters"]] <- as.numeric(factor(seuratObj@meta.data[["cnv_clusters"]]))

  return(seuratObj)
}
