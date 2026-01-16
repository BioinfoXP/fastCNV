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
#' CNVAnalysis
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
  
  # --- 1. 智能选择 Assay ---
  pick_assay <- function(obj, user_assay) {
    if(!is.null(user_assay)) return(user_assay)
    if("AggregatedCounts" %in% Seurat::Assays(obj)) {
      message(crayon::cyan("Note: Using 'AggregatedCounts' assay for analysis."))
      return("AggregatedCounts")
    }
    d_assay <- Seurat::DefaultAssay(obj)
    message(crayon::cyan(paste0("Note: 'AggregatedCounts' not found. Using default assay '", d_assay, "'.")))
    return(d_assay)
  }
  
  # --- 2. 核心运行函数 ---
  run_single <- function(obj) {
    use_assay <- pick_assay(obj, assay)
    
    # --- 【关键修复】更稳健的 Reference 检查 ---
    if(!is.null(referenceVar) && !is.null(referenceLabel)) {
      # A. 检查列是否存在
      if(!referenceVar %in% colnames(obj@meta.data)) {
        stop(paste0("Error: Metadata column '", referenceVar, "' not found in the Seurat object."))
      }
      
      # B. 安全获取 Metadata 向量 (强制转为 character 避免因子问题)
      # 直接从 @meta.data 获取，避开 Seurat [[ ]] 访问器的潜在版本差异
      meta_vals <- as.character(obj@meta.data[[referenceVar]])
      ref_vals <- as.character(referenceLabel)
      
      # C. 检查是否有匹配
      # 使用 any(...) 检查是否至少有一个细胞匹配
      if(!any(meta_vals %in% ref_vals)) {
        # 如果报错，打印出当前列里前5个值，帮助 debug
        found_vals <- paste(head(unique(meta_vals), 5), collapse = ", ")
        stop(paste0("Error: No cells found for reference label(s): ", paste(ref_vals, collapse=", "), 
                    ".\n  > Inside '", referenceVar, "' column, found values include: ", found_vals))
      }
    }
    # -----------------------------------------------
    
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
  
  # --- 3. 列表/单对象处理流程 ---
  if (!is.list(object)) {
    object <- run_single(object)
  } else {
    if (length(object) == 1) {
      object <- list(run_single(object[[1]]))
    } else {
      if (pooledReference) {
        # 假设列表所有对象结构一致
        use_assay <- pick_assay(object[[1]], assay)
        
        # 对列表中的每个对象进行简单的 Reference 预检查 (可选)
        for(i in seq_along(object)) {
             if(!is.null(referenceVar) && !referenceVar %in% colnames(object[[i]]@meta.data)) {
                 warning(paste0("Warning: Object ", i, " is missing referenceVar '", referenceVar, "'"))
             }
        }

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
        object <- lapply(object, run_single)
      }
    }
  }
  
  invisible(gc())
  message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  return(object)
}