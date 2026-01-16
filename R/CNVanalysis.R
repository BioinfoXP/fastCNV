#' CNVAnalysis
#'
#' Runs Copy Number Variation (CNV) analysis on a Seurat object or a list of Seurat objects.
#'
#' This function performs CNV analysis by calculating genomic scores. It is designed to work
#' downstream of `prepareCountsForCNVAnalysis`. It automatically detects if an "AggregatedCounts"
#' assay exists and uses it.
#'
#' @param object A Seurat object or a list of Seurat objects containing the data for CNV analysis. Each object
#' can be either **single-cell** or **spatial transcriptomics** data.
#' @param referenceVar The name of the metadata column in the Seurat object that contains reference annotations.
#' @param referenceLabel The label within `referenceVar` that specifies the reference population (can be a vector of labels, e.g. `c("T_cells", "B_cells")`).
#' @param pooledReference Logical. If `TRUE` (default), builds a pooled reference across all samples.
#' @param scaleOnReferenceLabel Logical. If `TRUE` (default), scales the results based on the reference population.
#' @param assay Name of the assay to run the CNV analysis on. Defaults to "AggregatedCounts" if available, otherwise uses the DefaultAssay.
#' @param thresholdPercentile Numeric. Specifies the quantile range to consider (e.g., `0.01` keeps values between the 1st and 99th percentiles). Higher values filter out more background noise.
#' @param geneMetadata A dataframe containing gene metadata, typically from Ensembl.
#' @param windowSize Integer. Defines the size of genomic windows for CNV analysis (default = 150 genes).
#' @param windowStep Integer. Specifies the step size between genomic windows (default = 10 genes).
#' @param saveGenomicWindows Logical. If `TRUE`, saves genomic window information in the current directory (default = `FALSE`).
#' @param topNGenes Integer. The number of top-expressed genes to retain in the analysis.
#' @param chrArmsToForce A chromosome arm (e.g., `"8p"`, `"3q"`) or a list of chromosome arms to force into the analysis.
#' @param genesToForce A list of genes to force into the analysis.
#' @param regionToForce Chromosome region to force into the analysis.
#'
#' @return If given a **single** Seurat object, returns the same object with:
#' - An **additional assay** containing genomic scores per genomic window.
#' - A new **CNV fraction column** added to the object’s metadata.
#' If given a **list** of Seurat objects, returns the modified list.
#'
#' @import Seurat
#' @importFrom crayon yellow green cyan red
#'
#' @examples
#' \dontrun{
#' # Assuming you have run prepareCountsForCNVAnalysis
#' # Run CNV analysis using "T_cells" and "B_cells" as the reference normal
#' my_seurat <- CNVAnalysis(
#'   object = my_seurat,
#'   referenceVar = "cell_type",
#'   referenceLabel = c("T_cells", "B_cells"),
#'   windowSize = 150
#' )
#' }
#'
#' @export
CNVAnalysis <- function(object,
                        referenceVar = NULL,
                        referenceLabel = NULL,
                        pooledReference = TRUE,
                        scaleOnReferenceLabel = TRUE,
                        assay = NULL,
                        thresholdPercentile = 0.01,
                        geneMetadata=getGenes(),
                        windowSize=150,
                        windowStep=10,
                        saveGenomicWindows = FALSE,
                        topNGenes=7000,
                        chrArmsToForce = NULL,
                        genesToForce = NULL,
                        regionToForce = NULL) {

  message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Running CNV analysis...")))

  # --- 1. 智能选择 Assay 的辅助函数 ---
  pick_assay <- function(obj, user_assay) {
      # 如果用户手动指定了，直接使用
      if(!is.null(user_assay)) return(user_assay)

      # 自动检测：如果存在聚合后的 Assay，优先使用
      if("AggregatedCounts" %in% Seurat::Assays(obj)) {
          message(crayon::cyan("Note: Using 'AggregatedCounts' assay for analysis."))
          return("AggregatedCounts")
      }

      # 否则使用默认
      d_assay <- Seurat::DefaultAssay(obj)
      message(crayon::cyan(paste0("Note: 'AggregatedCounts' not found. Using default assay '", d_assay, "'.")))
      return(d_assay)
  }

  # --- 2. 单样本运行包装器 ---
  run_single <- function(obj) {
      use_assay <- pick_assay(obj, assay)

      # 安全检查：确保 Reference Label 在该样本中存在
      # 避免传入错误的 Label 导致后续矩阵切片越界
      if(!is.null(referenceVar) && !is.null(referenceLabel)) {
          if(!referenceVar %in% colnames(obj@meta.data)) {
               stop(paste0("Error: referenceVar '", referenceVar, "' not found in metadata."))
          }
          if(sum(obj[[referenceVar]] %in% referenceLabel) == 0) {
               stop(paste0("Error: No cells found for reference label(s): ", paste(referenceLabel, collapse=", "), "."))
          }
      }

      # 调用核心计算函数
      CNVCalling(obj,
                 assay = use_assay,
                 referenceVar = referenceVar,
                 referenceLabel = referenceLabel,
                 scaleOnReferenceLabel = scaleOnReferenceLabel,
                 thresholdPercentile = thresholdPercentile,
                 geneMetadata = geneMetadata,
                 windowSize = windowSize,
                 windowStep = windowStep,
                 saveGenomicWindows = saveGenomicWindows,
                 topNGenes = topNGenes)
  }

  # --- 3. 主流程控制 (支持 Single Object 或 List) ---
  if (!is.list(object)) {
      # 单个对象
      object <- run_single(object)
  } else {
      # 列表对象
      if (length(object) == 1) {
          object <- list(run_single(object[[1]]))
      } else {
          if (pooledReference) {
             # 如果是 Pooled Reference，通常 CNVCallingList 会处理跨样本合并
             # 我们假设所有样本结构一致，使用第一个样本来检测 Assay 名称
             use_assay <- pick_assay(object[[1]], assay)

             object <- CNVCallingList(object,
                                      assay = use_assay,
                                      referenceVar = referenceVar,
                                      referenceLabel = referenceLabel,
                                      scaleOnReferenceLabel = scaleOnReferenceLabel,
                                      thresholdPercentile = thresholdPercentile,
                                      geneMetadata = geneMetadata,
                                      windowSize = windowSize,
                                      windowStep = windowStep,
                                      saveGenomicWindows = saveGenomicWindows,
                                      topNGenes = topNGenes)
          } else {
             # 逐个独立运行
             object <- lapply(object, run_single)
          }
      }
  }

  invisible(gc())
  message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  return(object)
}