#' CNVAnalysis
#' Runs Copy Number Variation (CNV) analysis on a Seurat object or a list of Seurat objects.
#'
#' This function performs CNV analysis by calculating genomic scores, applying optional denoising, and optionally
#' scaling the results based on a reference population. It processes single-cell or spatial transcriptomics data,
#' generating an additional assay with genomic scores and adding a new metadata column for CNV fractions.
#'
#' @param object A Seurat object or a list of Seurat objects containing the data for CNV analysis. Each object
#' can be either **single-cell** or **spatial transcriptomics** data.
#' @param referenceVar The name of the metadata column in the Seurat object that contains reference annotations.
#' @param referenceLabel The label within `referenceVar` that specifies the reference population (can be any type of annotation).
#' @param pooledReference Logical. If `TRUE` (default), builds a pooled reference across all samples.
#' @param scaleOnReferenceLabel Logical. If `TRUE` (default), scales the results based on the reference population.
#' @param assay Name of the assay to run the CNV analysis on. Defaults to the results of `prepareCountsForCNVAnalysis` if available.
#' @param thresholdPercentile Numeric. Specifies the quantile range to consider (e.g., `0.01` keeps values between the 1st and 99th percentiles). Higher values filter out more background noise.
#' @param geneMetadata A dataframe containing gene metadata, typically from Ensembl.
#' @param windowSize Integer. Defines the size of genomic windows for CNV analysis.
#' @param windowStep Integer. Specifies the step size between genomic windows.
#' @param saveGenomicWindows Logical. If `TRUE`, saves genomic window information in the current directory (default = `FALSE`).
#' @param topNGenes Integer. The number of top-expressed genes to retain in the analysis.
#' @param chrArmsToForce A chromosome arm (e.g., `"8p"`, `"3q"`) or a list of chromosome arms (e.g., `c("3q", "8p", "17p")`) to force into the analysis.
#' If specified, all genes within the given chromosome arm(s) will be included.
#' @param genesToForce A list of genes to force into the analysis (e.g. `c("FOXP3","MUC16","SAMD15")`).
#' @param regionToForce Chromosome region to force into the analysis (vector containing chr, start, end).
#'
#' @return If given a **single** Seurat object, returns the same object with:
#' - An **additional assay** containing genomic scores per genomic window.
#' - A new **CNV fraction column** added to the object’s metadata.
#' If given a **list** of Seurat objects, returns the modified list.
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

  # --- 内部辅助函数：自动选择正确的 Assay ---
  check_assay <- function(obj, current_assay) {
    if (is.null(current_assay)) {
      if ("AggregatedCounts" %in% Seurat::Assays(obj)) {
        message(crayon::cyan("Note: Using 'AggregatedCounts' assay for CNV analysis."))
        return("AggregatedCounts")
      } else {
        d_assay <- Seurat::DefaultAssay(obj)
        message(crayon::cyan(paste0("Note: 'AggregatedCounts' not found. Using default assay: '", d_assay, "'")))
        return(d_assay)
      }
    }
    return(current_assay)
  }

  # --- 内部辅助函数：检查 Reference 是否存在 ---
  check_reference <- function(obj, ref_var, ref_labels) {
    if (!is.null(ref_var)) {
      if (!ref_var %in% colnames(obj@meta.data)) {
        stop(paste0("Error: referenceVar '", ref_var, "' not found in metadata."))
      }
      if (!is.null(ref_labels)) {
        found_cells <- sum(obj@meta.data[[ref_var]] %in% ref_labels)
        if (found_cells == 0) {
          stop(paste0("Error: No cells found matching referenceLabel(s): ", paste(ref_labels, collapse=", "), 
                      " in column '", ref_var, "'."))
        }
      }
    }
  }
  # ----------------------------------------------

  if (!is.list(object)) {
    # 单个对象处理
    assay <- check_assay(object, assay)
    check_reference(object, referenceVar, referenceLabel)
    
    object <- CNVCalling(object,
                         assay = assay,
                         referenceVar = referenceVar,
                         referenceLabel = referenceLabel,
                         scaleOnReferenceLabel = scaleOnReferenceLabel,
                         thresholdPercentile = thresholdPercentile,
                         geneMetadata=geneMetadata,
                         windowSize=windowSize,
                         windowStep=windowStep,
                         saveGenomicWindows = saveGenomicWindows,
                         topNGenes=topNGenes)
    invisible(gc())
    
  } else {
    # 列表对象处理
    # 假设所有对象结构相似，检查第一个对象的 Assay
    if (length(object) > 0) {
      assay <- check_assay(object[[1]], assay)
    }

    if (length(object) == 1) {
      check_reference(object[[1]], referenceVar, referenceLabel)
      object <- list(CNVCalling(object[[1]],
                                assay = assay,
                                referenceVar = referenceVar,
                                referenceLabel = referenceLabel,
                                scaleOnReferenceLabel = scaleOnReferenceLabel,
                                thresholdPercentile = thresholdPercentile,
                                geneMetadata=geneMetadata,
                                windowSize=windowSize,
                                windowStep=windowStep,
                                saveGenomicWindows = saveGenomicWindows,
                                topNGenes=topNGenes))
      invisible(gc())
      
    } else {
      # 批量处理
      if (pooledReference == TRUE) {
        # 检查所有对象的 reference (可选，防止中间报错)
        for (i in seq_along(object)) {
             check_reference(object[[i]], referenceVar, referenceLabel)
        }
        
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
                                 topNGenes=topNGenes)
        invisible(gc())
        
      } else {
        object <- lapply(object, function(x) {
          check_reference(x, referenceVar, referenceLabel)
          CNVCalling(x,
                     assay = assay,
                     referenceVar = referenceVar,
                     referenceLabel = referenceLabel,
                     scaleOnReferenceLabel = scaleOnReferenceLabel,
                     thresholdPercentile = thresholdPercentile,
                     geneMetadata=geneMetadata,
                     windowSize=windowSize,
                     windowStep=windowStep,
                     saveGenomicWindows = saveGenomicWindows,
                     topNGenes=topNGenes) } )
        invisible(gc())
      }
    }
  }
  
  message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  return (object)
}