#' Aggregate Observations by Cell Type for CNV Analysis
#'
#' Aggregates observations with the same cell types to increase counts per observation,
#' improving Copy Number Variation (CNV) computation. This function is compatible with
#' both Seurat v4 and v5.
#'
#' @param seuratObj A Seurat object containing the data.
#' @param sampleName A character string specifying the sample name.
#' @param referenceVar The name of the metadata column in the Seurat object that contains reference annotations.
#' @param aggregateByVar Logical. If `TRUE` (default), aggregates observations based on `referenceVar` annotations.
#' @param aggregFactor Integer. The target number of counts per observation (default = `15 000`).
#' @param seuratClusterResolution Numeric. The resolution used for Seurat clustering (default = `0.8`).
#' @param reClusterSeurat Logical. If `TRUE`, re-runs clustering on the Seurat object.
#'
#' @return A Seurat object with:
#' - A new assay called **"AggregatedCounts"** containing the modified count matrix.
#' - **Seurat clusters** stored in the metadata.
#' - The "AggregatedCounts" assay is forced to Seurat v3 structure to ensure compatibility with downstream CNV tools.
#'
#' @import Seurat
#' @importFrom crayon black cyan
#' @importFrom Matrix Matrix
#'
#' @examples
#' \dontrun{
#' # Load your Seurat object
#' my_seurat <- readRDS("sample_data.rds")
#'
#' # Prepare counts for CNV (aggregating by cell type annotation)
#' my_seurat <- prepareCountsForCNVAnalysis(
#'   seuratObj = my_seurat,
#'   sampleName = "Patient1",
#'   referenceVar = "cell_type_annotation",
#'   aggregFactor = 15000
#' )
#'
#' # The object now has an "AggregatedCounts" assay ready for CNVAnalysis
#' }
#'
#' @export
prepareCountsForCNVAnalysis <- function(seuratObj,
                                        sampleName = NULL,
                                        referenceVar = NULL,
                                        aggregateByVar = T,
                                        aggregFactor=15000,
                                        seuratClusterResolution = 0.8,
                                        reClusterSeurat = F ){

  # 1. 获取当前默认 Assay (通常是 RNA 或 Spatial)
  assay <- Seurat::Assays(seuratObj)[1]

  # --- 内部辅助函数：兼容 v4/v5 的矩阵获取 ---
  # 目的：无论是在 v4 还是 v5 环境，都能安全地拿到原始 Count 矩阵
  get_counts_matrix <- function(obj, assay_name) {
    mat <- tryCatch({
      # 优先尝试 v4 方式 (slot)
      as.matrix(Seurat::GetAssayData(obj, assay = assay_name, slot = "counts"))
    }, error = function(e) {
      # 失败则尝试 v5 方式 (layer)
      as.matrix(Seurat::GetAssayData(obj, assay = assay_name, layer = "counts"))
    })
    return(mat)
  }
  # -------------------------------------

  # 2. 聚类逻辑 (如果需要)
  need_clustering <- (!"seurat_clusters" %in% colnames(seuratObj[[]])) | reClusterSeurat
  if (need_clustering) {
      msg <- if (!is.null(sampleName)) paste0("Processing sample ", sampleName, "...") else "Processing..."
      message(crayon::black(msg))
      # 确保 SCTransform 运行顺利
      seuratObj <- Seurat::SCTransform(seuratObj, assay = assay, verbose = F)
      seuratObj <- Seurat::RunPCA(seuratObj, assay = "SCT", verbose = F)
      seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = "pca", dims = 1:10, verbose = F)
      seuratObj <- Seurat::FindClusters(seuratObj, resolution = seuratClusterResolution, verbose = F)
      message(crayon::black("Clustering done."))
  }

  # 3. 获取用于聚合的基础矩阵
  # 如果指定了 referenceVar，我们通常希望基于原始数据聚合；否则基于当前 assay
  target_assay_name <- if(aggregateByVar && !is.null(referenceVar)) Seurat::Assays(seuratObj)[1] else assay
  countsMat <- get_counts_matrix(seuratObj, target_assay_name)

  # 4. 构建聚合分组逻辑
  # Case A: 仅按 Seurat Cluster 聚合
  if (is.null(referenceVar) || aggregateByVar == F) {
      group_split <- list(split(Seurat::Cells(seuratObj), Seurat::FetchData(seuratObj, vars = "seurat_clusters")))
  }
  # Case B: 按 Reference + Cluster 双重聚合
  else {
      refs <- split(Seurat::Cells(seuratObj), Seurat::FetchData(seuratObj, vars = referenceVar))
      group_split <- lapply(refs, function(x) split(x, Seurat::FetchData(seuratObj, vars = "seurat_clusters", cells = x)))
      group_split <- lapply(group_split, function(x) Filter(function(y) length(y) > 0, x)) # 移除空组
  }

  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = 0, col.name = "metaSpots")
  nbMetaSpots <- 1

  # 5. 执行聚合循环
  for (top_group in group_split) {
    for (cells in top_group) {
      if (length(cells) == 0) next

      # 获取这组细胞的 nCount 并排序 (用于贪婪聚合算法)
      nc_data <- Seurat::FetchData(seuratObj, vars = paste0("nCount_", assay), cells = cells)

      # 边界检查：防止空数据
      if(ncol(nc_data) == 0) next

      nc_vec <- nc_data[,1]
      names(nc_vec) <- rownames(nc_data)
      nc_vec <- sort(nc_vec)

      barcodes <- c()
      taille <- 0

      for(i in seq_along(nc_vec)) {
        taille <- taille + nc_vec[i]
        barcodes <- c(barcodes, names(nc_vec)[i])

        # 判断是否达到聚合阈值 (aggregFactor) 或 是一组中最后一个细胞
        if (taille >= aggregFactor || i == length(nc_vec)) {
          if (length(barcodes) > 1) {
             # 确保 barcodes 在矩阵列名中存在
             valid_bcs <- intersect(barcodes, colnames(countsMat))
             if(length(valid_bcs) > 1) {
               # 核心聚合：计算行和
               row_sum <- rowSums(countsMat[, valid_bcs, drop=FALSE])
               # 赋值回这些列 (所有这些细胞在矩阵中变得相同，模拟“MetaSpot”)
               countsMat[, valid_bcs] <- row_sum
             }
          }
          # 更新元数据标记
          seuratObj@meta.data[barcodes, "metaSpots"] <- nbMetaSpots
          nbMetaSpots <- nbMetaSpots + 1

          # 重置计数器
          taille <- 0
          barcodes <- c()
        }
      }
    }
  }

  # --- 6. 【核心兼容性修复】创建 Assay 对象 ---
  # 步骤 A: 必须转为稀疏矩阵 (dgCMatrix)
  sparse_counts <- Matrix::Matrix(countsMat, sparse = TRUE)

  # 步骤 B: 创建标准 Assay 对象
  aggregAssay <- Seurat::CreateAssayObject(counts = sparse_counts)

  # 步骤 C: 如果环境是 Seurat v5 (生成了 Assay5)，强制降级为 v3 格式
  # 这解决了 "subscript out of bounds" 错误，因为下游工具可能还在用旧的索引方式
  if (inherits(aggregAssay, "Assay5")) {
      message(crayon::cyan("Note: Converting AggregatedCounts to Seurat v3 assay structure for compatibility."))
      aggregAssay <- Seurat::ConvertAssay(aggregAssay, convert.to = "v3")
  }

  # 步骤 D: 同时填充 data slot (有些旧代码会读取 @data 而不是 @counts)
  aggregAssay <- Seurat::SetAssayData(aggregAssay, slot = "data", new.data = sparse_counts)

  # 步骤 E: 添加到对象
  seuratObj[["AggregatedCounts"]] <- aggregAssay

  invisible(gc())

  if (!is.null(sampleName)) {
    if (Seurat::Project(seuratObj) == "SeuratProject") {
       Seurat::Project(seuratObj) <- sampleName
    }}
  return(seuratObj)
}