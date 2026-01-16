#' Perform CNV Clustering with Seurat Object
#'
#' The `CNVcluster` function performs hierarchical clustering on a genomic score matrix extracted from a Seurat object.
#' It provides options for plotting a dendrogram, an elbow plot for optimal cluster determination,
#' and cluster visualization on the dendrogram. The resulting cluster assignments are stored in the Seurat object.
#'
#' @param seuratObj A Seurat object containing a "genomicScores" assay with a matrix of genomic scores for clustering.
#' @param k Optional. The number of clusters to cut the dendrogram into. If `NULL`, the optimal number of clusters is determined automatically using the elbow method.
#' @param h Optional. The height at which to cut the dendrogram for clustering. If both `k` and `h` are provided, `k` takes precedence.
#'
#' @details
#' The function computes a Manhattan distance matrix and performs hierarchical clustering using the Ward.D2 method.
#' If `k` is not provided, the elbow method is applied to determine the optimal number of clusters based on the within-cluster sum of squares (WSS).
#'
#' The clusters are assigned to the Seurat object under the metadata column `cnv_clusters`.
#'
#' @return A Seurat object with an additional metadata column, `cnv_clusters`, containing the cluster assignments.
#'
#' @importFrom proxy dist as.dist
#' @importFrom utils tail
#' @importFrom graphics abline text
#' @importFrom stats hclust cutree rect.hclust
#'
#' Perform CNV Clustering with Seurat Object
#' @export
CNVCluster <- function(seuratObj,
                       k = NULL,
                       h = NULL) {

  if (is.null(k)){kDetection = "automatic"}
  if (!is.null(k)){kDetection = "manual"}

  # --- 核心修复：兼容 v4/v5 的数据获取 ---
  # 尝试获取 genomicScores 的数据矩阵
  # 优先尝试 v4 标准 (slot="data")，失败则尝试 v5 (layer="data")
  # 同时也作为保险，尝试 "counts" 如果 data 为空
  mat <- tryCatch({
    as.matrix(Seurat::GetAssayData(seuratObj, assay = "genomicScores", slot = "data"))
  }, error = function(e) {
    tryCatch({
      as.matrix(Seurat::GetAssayData(seuratObj, assay = "genomicScores", layer = "data"))
    }, error = function(e2) {
      # 如果 data slot 也没拿到，尝试 counts
      as.matrix(Seurat::GetAssayData(seuratObj, assay = "genomicScores", slot = "counts"))
    })
  })

  # 转置矩阵 (Rows=Cells, Cols=Features) 用于距离计算
  genomicMatrix <- t(mat)

  # 根据细胞数量选择聚类策略
  if (nrow(genomicMatrix) < 30000) {
    # 小数据量：使用层次聚类 (Hierarchical Clustering)

    # 计算曼哈顿距离
    dist_cos <- proxy::dist(genomicMatrix, method = "Manhattan")
    dist_matrix <- proxy::as.dist(dist_cos)

    # Ward.D2 聚类
    hc <- stats::hclust(dist_matrix, method = "ward.D2")

    # 自动确定 k (Elbow Method)
    if(kDetection == "automatic") {
      dist_matrix_full <- as.matrix(dist_matrix)

      find_elbow <- function(x, y, sensitivity = 0.05) {
        x_norm <- (x - min(x)) / (max(x) - min(x))
        y_norm <- (y - min(y)) / (max(y) - min(y))
        slopes <- diff(y_norm) / diff(x_norm)
        slope_changes <- diff(slopes)
        significant_changes <- which(abs(slope_changes) > sensitivity * max(abs(slope_changes)))
        if (length(significant_changes) == 0) return(x[ceiling(length(x)/2)])

        first_quarter <- ceiling(length(x) / 4)
        elbow_point <- significant_changes[significant_changes > first_quarter][1]
        if (is.na(elbow_point)) elbow_point <- utils::tail(significant_changes, 1)
        return(x[elbow_point + 1])
      }

      # 测试 k=1 到 15
      k_values <- 1:min(15, nrow(genomicMatrix)-1)
      if(length(k_values) < 2) {
        k <- 1
      } else {
        wss <- sapply(k_values, function(k_val) {
          cl <- stats::cutree(hc, k = k_val)
          # 计算类内平方和 (WSS) 的近似替代
          sum(sapply(unique(cl), function(c) {
            idx <- which(cl == c)
            if(length(idx) > 1) sum(dist_matrix_full[idx, idx]^2) / (2 * length(idx)) else 0
          }))
        })
        k <- find_elbow(k_values, wss)
      }
    }

    # 剪枝获取聚类结果
    clusters <- stats::cutree(hc, k = k, h = h)
    seuratObj$cnv_clusters <- as.factor(clusters)

  } else {
    # 大数据量：使用 PCA + K-means

    # 运行 PCA
    # 注意：prcomp 期望 Rows=Samples(Cells), Cols=Features
    pca_res <- prcomp(genomicMatrix, center = TRUE, scale. = TRUE)
    X <- pca_res$x[, 1:min(30, ncol(pca_res$x))]

    set.seed(42)

    # 自动确定 k
    if (kDetection == "automatic") {
      wss <- numeric(15)
      for (ktest in 1:15) {
        km_temp <- kmeans(X, centers = ktest, iter.max = 100)
        wss[ktest] <- km_temp$tot.withinss
      }
      # 简单的 Elbow 检测：找二阶导数最大点
      diff1 <- diff(wss)
      diff2 <- diff(diff1)
      # 如果 wss 变化平缓，默认给一个值，否则找拐点
      if(length(diff2) > 0) k <- which.min(diff2) + 1 else k <- 3
    }

    # K-means 聚类
    km <- kmeans(X, centers = k, iter.max = 100)
    seuratObj$cnv_clusters <- NA
    seuratObj$cnv_clusters[rownames(genomicMatrix)] <- km$cluster
  }

  return(seuratObj)
}
