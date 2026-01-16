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
generateCNVClonesMatrix <- function(seuratObj, healthyClusters = NULL, values = "scores", cnv_thresh = 0.15) {

  # --- 1. 稳健的数据提取 (修复核心报错) ---
  chrom_arms_standard <- c(
    paste0(rep(1:22, each=2), c(".p", ".q")),
    "X.p", "X.q"
  )
  cnv_cols <- paste0(chrom_arms_standard, "_CNV")

  # 检查列是否存在
  available_cols <- intersect(cnv_cols, colnames(seuratObj@meta.data))

  if(length(available_cols) == 0) {
    stop("generateCNVClonesMatrix: No per-chromosome arm CNV columns found in metadata. Please run CNVPerChromosomeArm() first.")
  }

  # 按列名提取矩阵 (避免位置索引导致的错位)
  cnv_matrix <- as.matrix(seuratObj@meta.data[, available_cols, drop=FALSE])
  rownames(cnv_matrix) <- rownames(seuratObj@meta.data)
  cnv_matrix <- cnv_matrix[Seurat::Cells(seuratObj), , drop=FALSE]

  # --- 2. 计算 Cluster 均值 ---
  if(!"cnv_clusters" %in% colnames(seuratObj@meta.data)) {
    stop("generateCNVClonesMatrix: 'cnv_clusters' column not found in metadata.")
  }

  # 使用 $ 获取向量，避免数据框类型报错
  clusters_vec <- seuratObj$cnv_clusters
  unique_clusters <- sort(unique(clusters_vec[!is.na(clusters_vec)]))

  # 初始化结果矩阵
  cnv_matrix_clusters <- matrix(0, nrow = length(unique_clusters), ncol = ncol(cnv_matrix))
  rownames(cnv_matrix_clusters) <- unique_clusters
  colnames(cnv_matrix_clusters) <- colnames(cnv_matrix)

  for (i in seq_along(unique_clusters)) {
    cluster <- unique_clusters[i]
    # 使用向量比较
    cells <- rownames(seuratObj@meta.data)[which(clusters_vec == cluster)]

    if(length(cells) > 1) {
      cnv_matrix_clusters[i, ] <- colMeans(cnv_matrix[cells, , drop = FALSE], na.rm = TRUE)
    } else if(length(cells) == 1) {
      cnv_matrix_clusters[i, ] <- cnv_matrix[cells, ]
    }
  }

  # --- 3. 标签重命名 ---
  new_rownames <- paste0("Clone ", rownames(cnv_matrix_clusters))

  if (!is.null(healthyClusters)) {
    for (hc in healthyClusters) {
      # 确保匹配正确 (转换为字符比较)
      idx <- which(rownames(cnv_matrix_clusters) == as.character(hc))
      if(length(idx) > 0) {
        new_rownames[idx] <- paste0("Benign ", hc)
      }
    }
  }
  rownames(cnv_matrix_clusters) <- new_rownames

  # --- 4. 返回结果 ---
  if(values == "scores") {
    return(cnv_matrix_clusters)
  } else if (values == "calls") {
    alt_matrix <- matrix(0, nrow = nrow(cnv_matrix_clusters), ncol = ncol(cnv_matrix_clusters),
                         dimnames = dimnames(cnv_matrix_clusters))
    # 向量化赋值
    alt_matrix[cnv_matrix_clusters >= cnv_thresh] <- 1
    alt_matrix[cnv_matrix_clusters <= -cnv_thresh] <- -1
    return(alt_matrix)
  } else {
    stop("Non supported `values` value. Must be one of: 'scores', 'calls'")
  }
}
