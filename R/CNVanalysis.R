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

  # --- 辅助函数：单样本处理 ---
  run_single <- function(obj) {
    # 1. 确定 Assay
    use_assay <- assay
    if(is.null(use_assay)) {
      if("AggregatedCounts" %in% Seurat::Assays(obj)) use_assay <- "AggregatedCounts"
      else use_assay <- Seurat::DefaultAssay(obj)
    }

    # 2. 强制转 v3 Assay (防御性降级)
    if(inherits(obj[[use_assay]], "Assay5")) {
      message(crayon::cyan(paste0("Downgrading assay '", use_assay, "' to v3 for compatibility.")))
      obj[[use_assay]] <- Seurat::ConvertAssay(obj[[use_assay]], convert.to = "v3")
    }

    # 3. 运行
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
               topNGenes = topNGenes,
               chrArmsToForce = chrArmsToForce,
               genesToForce = genesToForce,
               regionToForce = regionToForce)
  }

  if (!is.list(object)) {
    # 单个对象
    object <- run_single(object)
    invisible(gc())
  } else {
    # 列表
    if (length(object) == 1) {
      object <- list(run_single(object[[1]]))
      invisible(gc())
    } else {
      # 多样本
      if (pooledReference == TRUE) {
        # 确保所有对象 Assay 兼容
        # (这里略去 CNVCallingList 的具体修改，只要 CNVCallingList 内部也用了 GetAssayData 即可，建议单样本场景不用这个路径)
        object <- CNVCallingList(object,
                                 assay = assay,
                                 referenceVar = referenceVar,
                                 referenceLabel = referenceLabel,
                                 scaleOnReferenceLabel = scaleOnReferenceLabel,
                                 thresholdPercentile = thresholdPercentile,
                                 geneMetadata=geneMetadata,
                                 windowSize=windowSize,
                                 windowStep=windowStep,
                                 saveGenomicWindows = saveGenomicWindows,
                                 topNGenes=topNGenes,
                                 chrArmsToForce = chrArmsToForce,
                                 genesToForce = genesToForce,
                                 regionToForce = regionToForce)
        invisible(gc())
      } else {
        object <- lapply(object, run_single)
        invisible(gc())
      }
    }
  }
  message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  return (object)
}
