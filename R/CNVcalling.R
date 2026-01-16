#' CNVCalling
#' Performs Copy Number Variation (CNV) analysis on a Seurat object.
#'
#' @param seuratObj A Seurat object containing the data for CNV analysis.
#' Can be either **single-cell** or **spatial transcriptomics** data.
#' @param assay Name of the assay to run the CNV analysis on. Defaults to the results of `prepareCountsForCNVAnalysis` if available.
#' @param referenceVar The name of the metadata column in the Seurat object that contains reference annotations.
#' @param referenceLabel The label within `referenceVar` that specifies the reference population (can be any type of annotation).
#' @param scaleOnReferenceLabel Logical. If `TRUE` (default), scales the results based on the reference population.
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
#' @return The same Seurat object provided in `seuratObj`, with:
#' - An **additional assay** containing genomic scores per genomic window.
#' - A new **CNV fraction column** added to the object’s metadata.
#'
#' @importFrom stats na.omit median
#' @import Seurat
#' @import SeuratObject
#' @importFrom crayon black
#'
#' @export
CNVCalling <- function(seuratObj,
                       assay = NULL,
                       referenceVar = NULL,
                       referenceLabel = NULL,
                       scaleOnReferenceLabel = TRUE,
                       thresholdPercentile = 0.01,
                       geneMetadata=getGenes(),
                       windowSize=150,
                       windowStep=10,
                       saveGenomicWindows = FALSE,
                       topNGenes=7000,
                       chrArmsToForce = NULL,
                       genesToForce = NULL,
                       regionToForce = NULL) {

  # --- 0. 确定 Assay ---
  if (is.null(assay)) {
    if ("AggregatedCounts" %in% Seurat::Assays(seuratObj)) {
      assay = "AggregatedCounts"
    } else {
      assay = Seurat::Assays(seuratObj)[1]
    }
  }

  # --- 1. 获取矩阵 (兼容 v4/v5 核心修复) ---
  rawCounts <- tryCatch({
    as.matrix(Seurat::GetAssayData(seuratObj, assay = assay, slot = "counts"))
  }, error = function(e) {
    as.matrix(Seurat::GetAssayData(seuratObj, assay = assay, layer = "counts"))
  })

  # 获取 Reference Cells (简化逻辑)
  if (is.null(referenceVar) || is.null(referenceLabel)){
    message(crayon::black("referenceVar/referenceLabel not found. Computing without reference."))
    scaleOnReferenceLabel = FALSE
  } else {
    if (length(referenceLabel) == 1) {
      referenceCells <- Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == referenceLabel)]
    } else {
      referenceCells <- list()
      for (i in referenceLabel){
        cells <- Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == i)]
        if(length(cells) >= 5) referenceCells[[i]] <- cells
      }
      if(length(referenceCells)==0) scaleOnReferenceLabel = FALSE
    }
    if (length(unlist(referenceCells)) == 0) scaleOnReferenceLabel = FALSE
  }

  # 准备基因信息
  geneMetadata <- geneMetadata[which(geneMetadata$gene_biotype %in% c("protein_coding","lncRNA") & geneMetadata$chromosome_name %in% c(1:22,"X") & geneMetadata$hgnc_symbol !=""),]
  geneMetadata$chromosome_num <- geneMetadata$chromosome_name
  geneMetadata$chromosome_num[which(geneMetadata$chromosome_num=="X")]<- 23
  geneMetadata$chromosome_num <- as.numeric(geneMetadata$chromosome_num)
  geneMetadata2 <- unique(geneMetadata[,c("hgnc_symbol","chromosome_num","start_position", "end_position", "chr_arm")])

  funTrim <- function(normcounts,lo=-3,up=3){
    t(apply(normcounts,1, function(z) {
      z[which(z < lo)] <- lo ; z[which(z > up)]<-up;z } ))
  }

  # --- 2. 基因过滤与匹配 ---
  if (nrow(rawCounts) < topNGenes) { topNGenes = nrow(rawCounts) }

  commonGenes <- intersect(rownames(rawCounts), geneMetadata2$hgnc_symbol)

  # 【防御性检查】
  if(length(commonGenes) < 10) {
    stop(paste0("Error: Low gene match. Matrix rows: ", paste(head(rownames(rawCounts),3), collapse=","),
                ". Metadata: ", paste(head(geneMetadata2$hgnc_symbol,3), collapse=",")))
  }

  rawCounts <- rawCounts[commonGenes, ]
  invisible(gc())

  # 计算平均表达量
  if (scaleOnReferenceLabel){
    ref_c <- intersect(unlist(referenceCells), colnames(rawCounts))
    if(length(ref_c)>0) averageExpression <- rowMeans(rawCounts[, ref_c, drop=FALSE])
    else averageExpression <- rowMeans(rawCounts)
  } else {
    averageExpression <- rowMeans(rawCounts)
  }

  topExprGenes <- commonGenes[order(averageExpression, decreasing = T)[1:topNGenes]]

  # 强制包含基因
  if(!is.null(genesToForce)) topExprGenes <- union(topExprGenes, intersect(commonGenes, genesToForce))
  if(!is.null(regionToForce)) {
    region_genes <- geneMetadata2 %>%
      filter(.data$chromosome_name == regionToForce[1], .data$start_position >= regionToForce[2], .data$end_position <= regionToForce[3]) %>%
      pull(.data$hgnc_symbol) %>% unique() %>% setdiff("")
    topExprGenes <- union(topExprGenes, intersect(commonGenes, region_genes))
  }

  # 确保染色体臂基因充足
  topExprGenes_metadata <- geneMetadata2[geneMetadata2$hgnc_symbol %in% topExprGenes, ]
  topExprGenes_metadata$chr_arm_full <- paste0(topExprGenes_metadata$chromosome_num, topExprGenes_metadata$chr_arm)
  genes_by_arm <- split(topExprGenes_metadata$hgnc_symbol, topExprGenes_metadata$chr_arm_full)

  for (arm in unique(geneMetadata2$chr_arm)) {
    if (!(arm %in% names(genes_by_arm))) genes_by_arm[[arm]] <- character(0)
    if (length(genes_by_arm[[arm]]) < 200) {
      remaining_genes <- commonGenes[!commonGenes %in% genes_by_arm[[arm]]]
      top_arm_genes <- remaining_genes[order(averageExpression[commonGenes %in% remaining_genes], decreasing = TRUE)[1:200]]
      genes_by_arm[[arm]] <- unique(c(genes_by_arm[[arm]], top_arm_genes))
    }
  }

  # 强制染色体臂
  if(!is.null(chrArmsToForce)){
    if (length(chrArmsToForce > 1)){
      for (chr in chrArmsToForce){
        genesToAdd <- geneMetadata2$hgnc_symbol[which(paste0(geneMetadata2$chromosome_num,geneMetadata2$chr_arm) == chr)]
        genes_by_arm[[chr]] <- genesToAdd
      }
    } else {
      genesToAdd <- geneMetadata2$hgnc_symbol[which(paste0(geneMetadata2$chromosome_num,geneMetadata2$chr_arm) == chrArmsToForce)]
      genes_by_arm[[chrArmsToForce]] <- genesToAdd
    }
  }

  final_selected_genes <- intersect(unique(unlist(genes_by_arm)), rownames(rawCounts))
  rawCounts <- rawCounts[final_selected_genes, ]
  invisible(gc())

  # 归一化
  normCounts <- log2(1+rawCounts)
  rm(rawCounts); invisible(gc())
  normCounts <- scale(normCounts, scale = FALSE)

  # Scale on Reference
  if (scaleOnReferenceLabel) {
    if (length(referenceLabel) == 1) {
      ref_c <- intersect(unlist(referenceCells), colnames(normCounts))
      scaleFactor <- rowMeans(normCounts[, ref_c, drop=FALSE])
    } else {
      scaleFactor <- list()
      for (i in referenceLabel) {
        ref_c <- intersect(referenceCells[[i]], colnames(normCounts))
        if(length(ref_c)>0) scaleFactor[[i]] <- rowMeans(normCounts[, ref_c, drop=FALSE])
      }
      if(length(scaleFactor) > 0) {
        scaleFactor <- do.call(rbind, scaleFactor)
        scaleFactor <- na.omit(scaleFactor)
        scaleFactor <- apply(scaleFactor, 2, median)
      } else {
        scaleFactor <- rowMeans(normCounts)
      }
    }
  } else {
    scaleFactor <- rowMeans(normCounts)
  }

  normCounts <- normCounts - scaleFactor
  normCounts <- funTrim(normCounts, lo = -3, up = 3)
  invisible(gc())

  # Genomic Windows
  geneMetadata2 <- geneMetadata2[which(geneMetadata2$hgnc_symbol %in% rownames(normCounts)),]
  geneMetadata2 <- geneMetadata2[order(geneMetadata2$chromosome_num,geneMetadata2$start_position),]

  genomicWindows <- lapply(c(1:23), function(chrom) {
    genesC <- geneMetadata2[which(geneMetadata2$chromosome_num == chrom),]
    chr_arms <- unique(genesC$chr_arm)
    chrom_windows <- list()
    for (arm in chr_arms) {
      genesArm <- genesC[which(genesC$chr_arm == arm),]
      N <- nrow(genesArm)
      iter <- round(windowSize / 2)
      if (N > windowSize) {
        gw <- lapply(seq(iter + 1, N - iter, by = windowStep), function(i) {
          as.character(unlist(genesArm[(i - iter):(i + iter), "hgnc_symbol"]))
        })
        names(gw) <- paste0(chrom, ".", arm, 1:length(gw))
      } else {
        gw <- list(as.character(unlist(genesArm$"hgnc_symbol")))
        names(gw) <- paste0(chrom, ".", arm, 1)
      }
      chrom_windows <- c(chrom_windows, gw)
    }
    return(chrom_windows)
  })
  genomicWindows <- unlist(genomicWindows,recursive=F)

  # 计算 Genomic Scores
  genomicScores <- sapply(genomicWindows, function(g) {
    g_safe <- intersect(g, rownames(normCounts))
    if (length(g_safe) == 0) return(rep(0, ncol(normCounts)))
    if (length(g_safe) == 1) normCounts[g_safe, ]
    else colMeans(normCounts[g_safe, ])
  })
  if(nrow(genomicScores) == ncol(normCounts)) genomicScores <- t(genomicScores)

  rm(normCounts); invisible(gc())

  # 截断处理
  if (scaleOnReferenceLabel) {
    if (length(referenceLabel) == 1) {
      ref_c <- intersect(unlist(referenceCells), colnames(genomicScores))
      genomicScoresReferenceLabel <- t(genomicScores[, ref_c, drop=FALSE])
    } else {
      genomicScoresReferenceLabel <- list()
      for (i in referenceLabel){
        ref_c <- intersect(referenceCells[[i]], colnames(genomicScores))
        if(length(ref_c)>0) genomicScoresReferenceLabel[[i]] <- t(genomicScores[, ref_c, drop=FALSE])
      }
      genomicScoresReferenceLabel <- do.call(rbind, genomicScoresReferenceLabel)
    }
    Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, "probs"=c(0+thresholdPercentile,1-thresholdPercentile))
  } else {
    Q01Q99 <- apply(genomicScores, 1, stats::quantile, "probs"=c(0+thresholdPercentile, 1-thresholdPercentile))
  }

  genomicScoresTrimmed <- genomicScores
  for(i in 1:nrow(genomicScores)) {
    low <- Q01Q99[1, i]; high <- Q01Q99[2, i]
    row_vals <- genomicScores[i, ]
    row_vals[row_vals >= low & row_vals <= high] <- 0
    genomicScoresTrimmed[i, ] <- row_vals
  }

  if (saveGenomicWindows){
    save(genomicWindows, file = paste0("genomicWindows_size",windowSize,"_step",windowStep,".RData"))
  }

  # --- 【关键修复】输出 v3 Assay ---
  rawGenomicAssay <- Seurat::CreateAssayObject(counts = genomicScores)
  if(inherits(rawGenomicAssay, "Assay5")) rawGenomicAssay <- Seurat::ConvertAssay(rawGenomicAssay, convert.to="v3")
  suppressWarnings({seuratObj[["rawGenomicScores"]] <- rawGenomicAssay})

  genomicAssay <- Seurat::CreateAssayObject(counts = genomicScoresTrimmed)
  # 强制转为 v3 (同时填充 data slot)
  if(inherits(genomicAssay, "Assay5")) genomicAssay <- Seurat::ConvertAssay(genomicAssay, convert.to="v3")
  genomicAssay <- Seurat::SetAssayData(genomicAssay, slot="data", new.data=as.matrix(genomicScoresTrimmed)) # 兼容 v4

  suppressWarnings({seuratObj[["genomicScores"]] <- genomicAssay})
  seuratObj[["cnv_fraction"]] <- colMeans(abs(genomicScoresTrimmed) > 0)

  invisible(gc())
  return (seuratObj)
}
