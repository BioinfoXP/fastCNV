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

  assay <- Seurat::Assays(seuratObj)[1]

  # --- 1. 获取矩阵 (兼容 v5 layer / v4 slot) ---
  get_matrix <- function(obj, target_assay) {
    tryCatch({
      as.matrix(Seurat::GetAssayData(obj, assay = target_assay, slot = "counts"))
    }, error = function(e) {
      # v5 fallback
      as.matrix(Seurat::GetAssayData(obj, assay = target_assay, layer = "counts"))
    })
  }

  # --- 2. 聚类流程 ---
  # 检查是否需要聚类
  need_clustering <- (!"seurat_clusters" %in% colnames(seuratObj[[]])) | reClusterSeurat
  if (need_clustering) {
    msg <- if (!is.null(sampleName)) paste0("Processing sample ", sampleName, "...") else "Processing..."
    message(crayon::black(msg))

    seuratObj <- Seurat::SCTransform(seuratObj, assay = assay, verbose = F)
    seuratObj <- Seurat::RunPCA(seuratObj, assay = "SCT", verbose = F)
    seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = "pca", dims = 1:10, verbose = F)
    seuratObj <- Seurat::FindClusters(seuratObj, resolution = seuratClusterResolution, verbose = F)
    message(crayon::black("Clustering done."))
  }

  # --- 3. 准备聚合数据 ---
  source_assay <- if(aggregateByVar && !is.null(referenceVar)) Seurat::Assays(seuratObj)[1] else assay

  countsMat <- get_matrix(seuratObj, source_assay)
  # 确保矩阵和对象细胞一致 (SCT可能过滤细胞)
  common_cells <- intersect(colnames(countsMat), Seurat::Cells(seuratObj))
  countsMat <- countsMat[, common_cells, drop=FALSE]
  seuratObj <- subset(seuratObj, cells = common_cells)

  # --- 4. 聚合逻辑 ---
  if (is.null(referenceVar) || aggregateByVar == F) {
    # 策略A: 仅按Cluster
    group_split <- list(split(Seurat::Cells(seuratObj), Seurat::FetchData(seuratObj, vars = "seurat_clusters")))
  } else {
    # 策略B: Ref + Cluster
    refs <- split(Seurat::Cells(seuratObj), Seurat::FetchData(seuratObj, vars = referenceVar))
    group_split <- lapply(refs, function(x) split(x, Seurat::FetchData(seuratObj, vars = "seurat_clusters", cells = x)))
    group_split <- lapply(group_split, function(x) Filter(function(y) length(y) > 0, x))
  }

  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = 0, col.name = "metaSpots")
  nbMetaSpots <- 1

  # 执行聚合
  for (top_group in group_split) {
    for (cells in top_group) {
      if (length(cells) == 0) next

      # 获取排序后的nCount
      nc_data <- Seurat::FetchData(seuratObj, vars = paste0("nCount_", assay), cells = cells)
      if(ncol(nc_data) == 0) next # 防御空数据
      nc_vec <- sort(nc_data[,1])
      names(nc_vec) <- rownames(nc_data)

      barcodes <- c()
      taille <- 0

      for(i in seq_along(nc_vec)) {
        taille <- taille + nc_vec[i]
        barcodes <- c(barcodes, names(nc_vec)[i])

        # 达到阈值或最后一组
        if (taille >= aggregFactor || i == length(nc_vec)) {
          if (length(barcodes) > 1) {
            # 聚合计算：行求和
            row_sum <- rowSums(countsMat[, barcodes, drop=FALSE])
            countsMat[, barcodes] <- row_sum
          }
          seuratObj@meta.data[barcodes, "metaSpots"] <- nbMetaSpots
          nbMetaSpots <- nbMetaSpots + 1
          taille <- 0
          barcodes <- c()
        }
      }
    }
  }

  # --- 5. 创建纯净的 v3 Assay (关键修复) ---
  sparse_counts <- Matrix::Matrix(countsMat, sparse = TRUE)
  aggregAssay <- Seurat::CreateAssayObject(counts = sparse_counts)

  # 强制降级检查：如果是 Assay5，转为 v3
  if (inherits(aggregAssay, "Assay5")) {
    aggregAssay <- Seurat::ConvertAssay(aggregAssay, convert.to = "v3")
  }

  # 填充 data slot (有些旧代码只读 data)
  aggregAssay <- Seurat::SetAssayData(aggregAssay, slot = "data", new.data = sparse_counts)
  seuratObj[["AggregatedCounts"]] <- aggregAssay

  invisible(gc())

  if (!is.null(sampleName)) {
    if (Seurat::Project(seuratObj) == "SeuratProject") {
      Seurat::Project(seuratObj) <- sampleName
    }}
  return(seuratObj)
}
