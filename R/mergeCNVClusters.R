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

  # --- 1. 稳健的数据提取 (保持之前的修复) ---
  chrom_arms_standard <- c(
    paste0(rep(1:22, each=2), c(".p", ".q")),
    "X.p", "X.q"
  )
  cnv_cols <- paste0(chrom_arms_standard, "_CNV")

  # 检查列是否存在
  available_cols <- intersect(cnv_cols, colnames(seuratObj@meta.data))

  if(length(available_cols) == 0) {
    warning("mergeCNVClusters: No per-chromosome arm CNV columns found in metadata. Skipping merge.")
    return(seuratObj)
  }

  # 按列名提取矩阵
  cnv_matrix <- as.matrix(seuratObj@meta.data[, available_cols, drop=FALSE])
  rownames(cnv_matrix) <- rownames(seuratObj@meta.data)
  cnv_matrix <- cnv_matrix[Seurat::Cells(seuratObj), , drop=FALSE]

  # --- 2. 计算平均 CNV (修复类型报错的核心) ---
  if (!"cnv_clusters" %in% colnames(seuratObj@meta.data)) {
    warning("mergeCNVClusters: 'cnv_clusters' column not found. Skipping.")
    return(seuratObj)
  }

  # 【关键修改】使用 $ 获取向量，而不是 [[ ]] 获取数据框
  clusters_vec <- seuratObj$cnv_clusters

  unique_clusters <- unique(clusters_vec)
  unique_clusters <- unique_clusters[!is.na(unique_clusters)]
  unique_clusters <- sort(unique_clusters)

  cnv_matrix_clusters <- matrix(0, nrow = length(unique_clusters), ncol = ncol(cnv_matrix))
  rownames(cnv_matrix_clusters) <- unique_clusters
  colnames(cnv_matrix_clusters) <- colnames(cnv_matrix)

  for (i in seq_along(unique_clusters)) {
    cluster <- unique_clusters[i]

    # 【关键修改】使用向量进行比较 (clusters_vec == cluster)
    cells <- rownames(seuratObj@meta.data)[which(clusters_vec == cluster)]

    if(length(cells) > 1) {
      cnv_matrix_clusters[i, ] <- colMeans(cnv_matrix[cells, , drop = FALSE], na.rm = TRUE)
    } else if(length(cells) == 1) {
      cnv_matrix_clusters[i, ] <- cnv_matrix[cells, ]
    }
  }

  # --- 3. 计算相关性并合并 (保持之前的修复) ---
  # 移除全为0的列
  keep_cols <- colSums(abs(cnv_matrix_clusters)) > 0

  if(sum(keep_cols) < 2) {
    return(seuratObj)
  }

  cnv_matrix_clusters_clean <- cnv_matrix_clusters[, keep_cols, drop=FALSE]

  # 计算相关性 (处理 NA)
  corrmat <- cor(t(cnv_matrix_clusters_clean), use = "pairwise.complete.obs")
  corrmat[is.na(corrmat)] <- 0

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

  # --- 4. 更新 Cluster ID ---
  merged <- split(cluster_names, groups)

  mapping <- unlist(lapply(merged, function(grp) {
    new_name <- grp[1]
    setNames(rep(new_name, length(grp)), grp)
  }))

  # 这里也可以用 seuratObj$cnv_clusters，更安全
  orig <- as.character(seuratObj$cnv_clusters)
  mapped <- mapping[orig]
  mapped[is.na(mapped)] <- orig[is.na(mapped)]

  seuratObj@meta.data[["cnv_clusters"]] <- as.numeric(as.factor(mapped))

  return(seuratObj)
}
