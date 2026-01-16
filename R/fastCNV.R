#' fastCNV calls all of the internal functions needed to compute the putative CNV on a Seurat object or a list of Seurat objects
#'
#' This function orchestrates the CNV analysis on a Seurat object (or multiple objects). It calls internal functions such as
#' `prepareCountsForCNVAnalysis`, `CNVAnalysis`, `CNVPerChromosomeArm`, `CNVcluster`, and `PlotCNVResults` to compute the CNVs,
#' perform clustering, and generate heatmaps. The results are saved in the metadata of the Seurat object(s), with options for
#' generating and saving plots.
#'
#' @param seuratObj Seurat object or list of Seurat objects to perform the CNV analysis on.
#' @param sampleName Name of the sample or a list of names corresponding to the samples in the `seuratObj`.
#' @param referenceVar The variable name of the annotations in the Seurat metadata to be used as reference.
#' @param referenceLabel The label given to the observations you want as reference (can be any type of annotation).
#' @param assay Name of the assay to run the CNV on. Takes the results of `prepareCountsForCNVAnalysis` by default if available.
#' @param prepareCounts If `FALSE`, will not run the `prepareCountsForCNVAnalysis` function (default = `TRUE`).
#' @param aggregFactor The number of counts per spot desired (default = 15 000). If less than 1,000, will not run the `prepareCountsForCNVAnalysis` function.
#' @param seuratClusterResolution The resolution wanted for the Seurat clusters (default = 0.8).
#' @param aggregateByVar If `referenceVar` is given, determines whether to use it to pool the observations (default = `TRUE`).
#' @param reClusterSeurat Whether to re-cluster if the Seurat object given already has a `seurat_clusters` slot in its metadata (default = `FALSE`).
#' @param pooledReference Default is `TRUE`. Will build a pooled reference across all samples if `TRUE`.
#' @param denoise If `TRUE`, the denoised data will be used in the heatmap (default = `TRUE`).
#' @param scaleOnReferenceLabel If `TRUE`, scales the results depending on the normal observations (default = `TRUE`).
#' @param thresholdPercentile Which quantiles to take (default 0.01). For example, `0.01` will take quantiles between 0.01-0.99. Background noise appears with higher numbers.
#' @param geneMetadata List of genes and their metadata (default uses genes from Ensembl version 113).
#' @param chrArmsToForce A chromosome arm (e.g., `"8p"`, `"3q"`) or a list of chromosome arms (e.g., `c("3q", "8p", "17p")`) to force into the analysis.
#' @param genesToForce A list of genes to force into the analysis (e.g. `c("FOXP3","MUC16","SAMD15")`).
#' @param regionToForce Chromosome region to force into the analysis (vector containing chr, start, end).
#' @param windowSize Size of the genomic windows for CNV analysis (default = 150).
#' @param windowStep Step between the genomic windows (default = 10).
#' @param saveGenomicWindows If `TRUE`, saves the information of the genomic windows in the current directory (default = `FALSE`).
#' @param topNGenes Number of top expressed genes to keep (default = 7000).
#' @param getCNVPerChromosomeArm If `TRUE`, will save the CNV per chromosome arm into the metadata.
#' @param getCNVClusters If `TRUE`, will perform clustering on the CNV scores and save them in the metadata of the Seurat object as `cnv_clusters`.
#' @param k_clusters Optional. Number of clusters to cut the dendrogram into. If `NULL`, the optimal number of clusters is determined automatically using the elbow method.
#' @param h_clusters Optional. The height at which to cut the dendrogram for clustering. If both `k` and `h` are provided, `k` takes precedence.
#' @param mergeCNV Logical. Whether to merge the highly correlated CNV clusters.
#' @param mergeThreshold A numeric value between 0 and 1. Clusters with correlation greater than this threshold will be merged. Default is 0.98.
#' @param doPlot If `TRUE`, will build a heatmap for each of the samples (default = `TRUE`).
#' @param printPlot If `TRUE`, the heatmap will be printed in the console (default = `FALSE`, the plot will only be saved in a PDF).
#' @param savePath Path to save the heatmap plot. If `NULL`, the plot won't be saved (default = `.`).
#' @param outputType Specifies the file format for saving the plot, either `"png"` or `"pdf"` (default = `"png"`).
#' @param clustersVar The variable name of the clusters in the Seurat metadata (default = `"cnv_clusters"`).
#' @param splitPlotOnVar The name of the metadata column to split the observations during the `plotCNVResults` step, if different from `referenceVar`.
#' @param referencePalette The color palette that should be used for `referenceVar` (default = `"default"`).
#' @param clusters_palette The color palette that should be used for `clustersVar` (default = `"default"`).
#'
#' @return A list of Seurat objects after all the analysis is complete. Heatmaps of the CNVs for every object in `seuratObj` are generated and saved in the specified path (default = current working directory).
#'
#' @importFrom crayon red yellow green black
#'
#' @export
fastCNV <- function (seuratObj, sampleName, referenceVar = NULL, referenceLabel = NULL, assay = NULL,
                     prepareCounts = TRUE, aggregFactor = 15000, seuratClusterResolution = 0.8,
                     aggregateByVar = TRUE, reClusterSeurat = FALSE, pooledReference = TRUE,
                     scaleOnReferenceLabel = TRUE, thresholdPercentile = 0.01, geneMetadata = getGenes(),
                     windowSize = 150, windowStep = 10, saveGenomicWindows = FALSE, topNGenes = 7000,
                     chrArmsToForce = NULL, genesToForce = NULL, regionToForce = NULL,
                     getCNVPerChromosomeArm = TRUE, getCNVClusters = TRUE, k_clusters = NULL, h_clusters = NULL,
                     mergeCNV = TRUE, mergeThreshold = 0.98, doPlot = TRUE, denoise = TRUE, printPlot = FALSE,
                     savePath = ".", outputType = "png", clustersVar = "cnv_clusters", splitPlotOnVar = clustersVar,
                     referencePalette = "default", clusters_palette = "default"){

  if(!length(seuratObj)==length(sampleName)) stop("seuratObj & sampleName length mismatch")
  options(future.globals.maxSize = 8000*1024^2)

  # 列表标准化
  if(!is.list(seuratObj) || inherits(seuratObj, "Seurat")){
    seuratObj <- list(seuratObj)
    names(seuratObj) <- sampleName
  }
  for (i in 1:length(seuratObj)) seuratObj[[i]]@project.name = sampleName[[i]]

  # 追踪 Assay
  use_assay <- assay

  # 1. 运行 PrepareCounts
  if (prepareCounts == TRUE & aggregFactor >= 1000) {
    message(crayon::yellow(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Aggregating counts matrix...")))
    for (i in 1:length(seuratObj)) {
      seuratObj[[i]] <- prepareCountsForCNVAnalysis(seuratObj[[i]], sampleName = sampleName[[i]],
                                                    referenceVar = referenceVar, aggregateByVar = aggregateByVar,
                                                    aggregFactor = aggregFactor, seuratClusterResolution = seuratClusterResolution,
                                                    reClusterSeurat = reClusterSeurat)
      invisible(gc())
    }
    use_assay <- "AggregatedCounts" # 明确更新
    message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  }

  # 2. 运行 CNV 分析
  # 【单样本保护】如果是单样本，强制关闭 pooledReference
  if(length(seuratObj) == 1) {
    pooledReference <- FALSE
  }

  seuratObj <- CNVAnalysis(seuratObj,
                           referenceVar = referenceVar,
                           referenceLabel = referenceLabel,
                           pooledReference = pooledReference, # 动态调整
                           scaleOnReferenceLabel = scaleOnReferenceLabel,
                           assay = use_assay, # 明确传递
                           thresholdPercentile = thresholdPercentile,
                           geneMetadata = geneMetadata,
                           windowSize = windowSize,
                           windowStep = windowStep,
                           saveGenomicWindows = saveGenomicWindows,
                           topNGenes = topNGenes,
                           chrArmsToForce = chrArmsToForce,
                           genesToForce = genesToForce,
                           regionToForce = regionToForce)
  invisible(gc())

  # 3. 后续步骤 (Chromosome Arm, Clustering, Plotting)
  # (保持原样，省略以节省空间，直接调用即可)

  if (getCNVPerChromosomeArm) {
    message(crayon::yellow("Computing CNV per chromosome arm..."))
    for (i in seq_along(seuratObj)) seuratObj[[i]] <- CNVPerChromosomeArm(seuratObj[[i]])
  }

  if (getCNVClusters) {
    message(crayon::yellow("Clustering CNVs..."))
    for (i in seq_along(seuratObj)) seuratObj[[i]] <- CNVCluster(seuratObj[[i]], k = k_clusters, h = h_clusters)
    if (mergeCNV) {
      for (i in seq_along(seuratObj)) seuratObj[[i]] <- mergeCNVClusters(seuratObj = seuratObj[[i]], mergeThreshold = mergeThreshold)
    }
  }

  if (doPlot) {
    message(crayon::yellow("Plotting CNV heatmap..."))
    for (i in seq_along(seuratObj)) {
      # 确保 Project Name
      if(Seurat::Project(seuratObj[[i]]) == "SeuratProject") Seurat::Project(seuratObj[[i]]) <- paste0("Sample", i)

      # 确定 Split 变量
      p_var <- splitPlotOnVar
      if ("cnv_clusters" %in% names(seuratObj[[i]]@meta.data)) p_var <- "cnv_clusters"

      plotCNVResults(seuratObj[[i]], referenceVar = referenceVar, splitPlotOnVar = p_var,
                     clustersVar = clustersVar, savePath = savePath, printPlot = printPlot,
                     referencePalette = referencePalette, clusters_palette = clusters_palette,
                     outputType = outputType, denoise = denoise)
    }
    message(crayon::green("Done !"))
  }

  if (length(seuratObj) == 1) return(seuratObj[[1]])
  return(seuratObj)
}
