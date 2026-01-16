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
#'
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
  # getting reference cells / spots
  if (is.null(referenceVar) || is.null(referenceLabel)){
    message(crayon::black,paste0("referenceVar and/or referenceLabel parameters not found. Computing the CNV without a reference."))
    # unable to scale the results on reference data if we don't know what the reference data is
    scaleOnReferenceLabel = FALSE
  } else {
    if (length(referenceLabel) == 1) {
      referenceCells <- Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == referenceLabel)]
      if (length(referenceCells) == 0) {
        message(crayon::black,paste0("CNVcalling : there is no annotation called ", referenceLabel," in the ", referenceVar," metadata slot of your seurat object.
Computing the CNV without a reference."))
        scaleOnReferenceLabel = FALSE
        rm(referenceCells)
      }
    } else {
      referenceCells <- list()
      for (i in referenceLabel){
        if (length(Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == i)]) >= 5) {
                  referenceCells[[i]] <- Seurat::Cells(seuratObj)[which(Seurat::FetchData(seuratObj, vars = referenceVar) == i)]
        }
      }
      if (length(referenceCells) == 0) {
        message(crayon::black,"Couldn't find any cells annotated as the referenceLabel in the ",referenceVar," metadata slot of the seurat object.
                Computing the CNV without a reference.")
        referenceLabel = NULL
        scaleOnReferenceLabel = FALSE
        rm(referenceCells)
      }
    }
  }

  # preparation of gene information
  geneMetadata <- geneMetadata[which(geneMetadata$gene_biotype %in% c("protein_coding","lncRNA") & geneMetadata$chromosome_name %in% c(1:22,"X") & geneMetadata$hgnc_symbol !=""),]
  geneMetadata$chromosome_num <- geneMetadata$chromosome_name
  geneMetadata$chromosome_num[which(geneMetadata$chromosome_num=="X")]<- 23
  geneMetadata$chromosome_num <- as.numeric(geneMetadata$chromosome_num)
  geneMetadata2 <- unique(geneMetadata[,c("hgnc_symbol","chromosome_num","start_position", "end_position", "chr_arm")])

  # internal functions
  funTrim <- function(normcounts,lo=-3,up=3){
    t(apply(normcounts,1, function(z) {
      z[which(z < lo)] <- lo ; z[which(z > up)]<-up;z } ))
  }


  # analysis
  if (is.null(assay)) {
    if ("AggregatedCounts" %in% Seurat::Assays(seuratObj)) {
      assay = "AggregatedCounts"
    } else {
      assay = Seurat::Assays(seuratObj)[1]
    }
  }

  if (dim(Seurat::GetAssay(seuratObj, assay = assay))[1] < topNGenes) {topNGenes = as.numeric(dim(Seurat::GetAssay(seuratObj, assay = assay))[1])}

  rawCounts <- as.matrix(Seurat::GetAssay(seuratObj, assay = assay)["counts"])
  invisible(gc())
  commonGenes <- intersect(rownames(rawCounts),geneMetadata2$hgnc_symbol)
  rawCounts <- rawCounts[commonGenes,]

  invisible(gc())

  if (!is.null(referenceVar) && !is.null(referenceLabel)){
    averageExpression <- rowMeans(rawCounts[,unlist(referenceCells)])
  } else {
    averageExpression <- rowMeans(rawCounts)
  }

  invisible(gc())

  topExprGenes <- commonGenes[order(averageExpression, decreasing = T)[1:topNGenes]]

  if(!is.null(genesToForce)) {
    topExprGenes <- union(topExprGenes, intersect(commonGenes, genesToForce))
  }

  if(!is.null(regionToForce)) {
    region_genes <- geneMetadata2 %>%
      filter(.data$chromosome_name == regionToForce[1], .data$start_position >= regionToForce[2], .data$end_position <= regionToForce[3]) %>%
      pull(.data$hgnc_symbol) %>%
      unique() %>%
      setdiff("")
    topExprGenes <- union(topExprGenes, intersect(commonGenes, region_genes))
  }

  topExprGenes_metadata <- geneMetadata2[geneMetadata2$hgnc_symbol %in% topExprGenes, ]
  topExprGenes_metadata$chr_arm_full <- paste0(topExprGenes_metadata$chromosome_num, topExprGenes_metadata$chr_arm)

  genes_by_arm <- split(
    topExprGenes_metadata$hgnc_symbol,
    topExprGenes_metadata$chr_arm_full
  )

  for (arm in unique(geneMetadata2$chr_arm)) {
    if (!(arm %in% names(genes_by_arm))) {
      genes_by_arm[[arm]] <- character(0)
    }
    if (length(genes_by_arm[[arm]]) < 200) {
      remaining_genes <- commonGenes[!commonGenes %in% genes_by_arm[[arm]]]
      top_arm_genes <- remaining_genes[order(averageExpression[commonGenes %in% remaining_genes], decreasing = TRUE)[1:200]]
      genes_by_arm[[arm]] <- unique(c(genes_by_arm[[arm]], top_arm_genes))
    }
  }

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

  final_selected_genes <- unlist(genes_by_arm)

  rawCounts <- rawCounts[final_selected_genes,]

  invisible(gc())

  normCounts <- log2(1+rawCounts)
  invisible(gc())
  normCounts <- scale(normCounts, scale = FALSE)

  rm(rawCounts)
  invisible(gc())

  if (scaleOnReferenceLabel) {
    if (length(referenceLabel) == 1) {
      scaleFactor <- rowMeans(normCounts[,referenceCells])
    } else {
      scaleFactor <- list()
      for (i in referenceLabel) {
        scaleFactor[[i]] <- rowMeans(normCounts[,referenceCells[[i]]])
      }
      scaleFactor <- do.call(rbind,scaleFactor)
      scaleFactor <- na.omit(scaleFactor)
      scaleFactor <- apply(scaleFactor, 2, median)
    }
  } else {
    scaleFactor <- rowMeans(normCounts)
  }

  invisible(gc())

  normCounts <- normCounts - scaleFactor
  invisible(gc())
  normCounts <- funTrim(normCounts, lo = -3, up = 3)
  invisible(gc())

  geneMetadata2 <- geneMetadata2[which(geneMetadata2$hgnc_symbol %in% topExprGenes),]
  geneMetadata2 <- geneMetadata2[order(geneMetadata2$chromosome_num,geneMetadata2$start_position),]

  # Preparation of genomic windows
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
          genesArm[(i - iter):(i + iter), "hgnc_symbol"] |> unlist() |> as.character()
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
  genomicScores <- sapply(genomicWindows, function(g) {
    if (length(g) == 1) {
      normCounts[g, ]
    } else {
      colMeans(normCounts[g, ])
    }
  })
  rm(normCounts)
  invisible(gc())

  if (scaleOnReferenceLabel) {
    if (length(referenceLabel) == 1) {
      genomicScoresReferenceLabel<- genomicScores[referenceCells,]
    } else {
      genomicScoresReferenceLabel <- list()
      for (i in referenceLabel){
        genomicScoresReferenceLabel[[i]] <- genomicScores[referenceCells[[i]],]
      }
      genomicScoresReferenceLabel <- do.call(rbind, genomicScoresReferenceLabel)
    }

    Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, "probs"=c(0+thresholdPercentile,1-thresholdPercentile))
    genomicScoresTrimmed <- apply(genomicScores, 1,function(v)
    {v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v })

  } else {
    if (is.null(referenceVar)) {
      Q01Q99 <- apply(genomicScores,2,stats::quantile,"probs"=c(0+thresholdPercentile,1-thresholdPercentile))
      genomicScoresTrimmed <- apply(genomicScores,1,function(v)
      {v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0;v})
    } else {
      cellLines <- split(Seurat::Cells(seuratObj), Seurat::FetchData(seuratObj, vars = referenceVar))

      high_threshold <- median(unlist(sapply(cellLines, function(z) apply (genomicScores[z, , drop = F], 2, function(v) quantile(v, probs = c(1-thresholdPercentile))))))
      low_threshold <- median(unlist(sapply(cellLines, function(z) apply (genomicScores[z, , drop = F], 2, function(v) quantile(v, probs = c(0+thresholdPercentile))))))
      genomicScoresTrimmed <- t(apply(genomicScores, 2, function(z){
        z[which( z>low_threshold & z<high_threshold)] = 0; z}))
    }
  }

  if (saveGenomicWindows){
    save(genomicWindows, file = paste0("genomicWindows_size",windowSize,"_step",windowStep,".RData"))
  }

  rawGenomicAssay <- Seurat::CreateAssayObject(counts = t(as.matrix(genomicScores)))
  suppressWarnings({seuratObj[["rawGenomicScores"]] <- rawGenomicAssay})
  genomicAssay <- Seurat::CreateAssayObject(counts = as.matrix(genomicScoresTrimmed))
  suppressWarnings({seuratObj[["genomicScores"]] <- genomicAssay})
  seuratObj[["cnv_fraction"]] <- colMeans(abs(genomicScoresTrimmed) > 0)

  invisible(gc())

  return (seuratObj)
}
