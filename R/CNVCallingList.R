#' CNVCalling for a List of Seurat Objects
#' Performs Copy Number Variation (CNV) analysis on a list of Seurat objects.
#'
#' @param seuratList A list of Seurat objects containing the data for CNV analysis.
#' Each object can be either **single-cell** or **spatial transcriptomics** data.
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
#' @return A list of Seurat objects, where each:
#' - Contains an **additional assay** with genomic scores per genomic window.
#' - Has a new **CNV fraction column** added to its metadata.
#'
#' @importFrom stats median
#' @import Seurat
#' @importFrom crayon black
#'
#' @export
#'


CNVCallingList <- function(seuratList,
                           assay = NULL,
                           referenceVar = NULL,
                           referenceLabel = NULL,
                           scaleOnReferenceLabel= TRUE,
                           thresholdPercentile = 0.01,
                           geneMetadata=getGenes(),
                           windowSize=150,
                           windowStep=10,
                           saveGenomicWindows = FALSE,
                           topNGenes=7000,
                           chrArmsToForce = NULL,
                           genesToForce = NULL,
                           regionToForce = NULL
){
  LrawcountsByPatient <- lapply(seuratList, function(x) {
    if (!is.null(assay)) {
      as.matrix(Seurat::GetAssay(x, assay = assay)["counts"])
    } else if ("AggregatedCounts" %in% Seurat::Assays(x)) {
      as.matrix(Seurat::GetAssay(x, assay = "AggregatedCounts")["counts"])
    } else {
      as.matrix(Seurat::GetAssay(x, assay = Seurat::Assays(x)[1])["counts"])
    } } )

  invisible(gc())

  names(LrawcountsByPatient) <- lapply(seuratList, function(x) Seurat::Project(x))
  Lannot <- lapply(seuratList, function(x) Seurat::FetchData(x, vars = referenceVar))
  names(Lannot) <- lapply(seuratList, function(x) Seurat::Project(x))
  invisible(gc())

  # getting reference items per patient
  if (is.null(referenceVar) || is.null(referenceLabel)){
    # unable to scale the results on reference data if we don't know what the reference data is
    message(crayon::black("referenceVar and/or referenceLabel parameters not found. Computing the CNV without a reference."))
    scaleOnReferenceLabel = FALSE
  } else if ( as.numeric(sum(sapply(lapply(Lannot, function(annot) rownames(annot)[which(annot == referenceLabel)]), function(x) length(x)))) == 0 ) {
    message(crayon::black("No observations found for the referenceLabel given. Computing the CNV without a reference."))
    scaleOnReferenceLabel = FALSE
  } else if (length(referenceLabel) == 1){
    LN <- lapply(Lannot, function(annot) rownames(annot)[which(annot == referenceLabel)])
    names(LN) <- names(Lannot)
    LN <- LN[which(sapply(LN,length) >= 5)]
    LN <- Filter(function(x) length(x) > 0, LN)
  } else {
    LN <- list()
    for (i in referenceLabel){
      LN[[i]] <- lapply(Lannot, function(annot) rownames(annot)[which(annot == i)])
      LN[[i]] <- LN[[i]][which(sapply(LN[[i]],length) >= 5)]
      LN <- Filter(function(x) length(x) > 0, LN)
    }
  }

  # preparation of gene information
  geneMetadata <- geneMetadata[which(geneMetadata$gene_biotype %in% c("protein_coding","lncRNA") & geneMetadata$chromosome_name %in% c(1:22,"X") & geneMetadata$hgnc_symbol !=""),]
  geneMetadata$chromosome_num <- geneMetadata$chromosome_name
  geneMetadata$chromosome_num[which(geneMetadata$chromosome_num=="X")]<- 23
  geneMetadata$chromosome_num <- as.numeric(geneMetadata$chromosome_num)
  geneMetadata2 <- unique(geneMetadata[,c("hgnc_symbol","chromosome_num","start_position", "chr_arm")])

  invisible(gc())

  # internal functions
  funTrim <- function(normcounts,lo=-3,up=3){
    t(apply(normcounts, 1, function(z) {z[which(z < lo)]<-lo ; z[which(z > up)]<-up ; z}))
  }

  funGenomicScore <- function(normcounts,GW=genomicWindows){
    res <- sapply(GW, function(g) colMeans(normcounts[g,]))
    res
  }

  # Analyses

  commonGenes <- Reduce(intersect, c(lapply(LrawcountsByPatient, rownames), list(geneMetadata2$hgnc_symbol)))
  LrawcountsByPatient <- lapply(LrawcountsByPatient, function(x) x[commonGenes,])
  invisible(gc())
  if (exists("LN")){
    aveExpr <- computeAverageExpression(LN, LrawcountsByPatient)
    aveExpr <- do.call(cbind,aveExpr)
    aveExpr <- rowMeans(aveExpr, na.rm = T)
  } else {
    aveExpr <- rowMeans(sapply(LrawcountsByPatient, rowMeans))
  }

  if (length(aveExpr) < topNGenes) {topNGenes = length(aveExpr)}
  topExprGenes <- commonGenes[order(aveExpr, decreasing = T)[1:topNGenes]]

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
      top_arm_genes <- remaining_genes[order(aveExpr[commonGenes %in% remaining_genes], decreasing = TRUE)[1:200]]
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

  LrawcountsByPatient <- lapply(LrawcountsByPatient, function(x) x[final_selected_genes,])
  invisible(gc())
  LnormcountsByPatient <- lapply(LrawcountsByPatient, function(d) log2(1+d))
  rm(LrawcountsByPatient) ; invisible(gc())
  LnormcountsByPatient <- lapply(LnormcountsByPatient, scale, "scale"=F)

  if(scaleOnReferenceLabel){
    if (length(referenceLabel) == 1) {
      SF <- rowMeans(sapply(names(LN), function(patient) rowMeans(LnormcountsByPatient[[patient]][,LN[[patient]]])))
    } else {
      if (length(LN) == 1) {
        SF <- rowMeans(sapply(names(LN[[1]]), function(patient) rowMeans(LnormcountsByPatient[[patient]][,LN[[1]][[patient]]])))
      } else {
        SF <- list()
        for (i in names(LN)) {
          SF[[i]] <- rowMeans(sapply(names(LN[[i]]), function(patient) rowMeans(LnormcountsByPatient[[patient]][,LN[[i]][[patient]]])))
        }
        SF <- do.call(rbind,SF)
        SF <- apply(SF, 2, median)
      }
    }
  }else{
    SF <- rowMeans(sapply(LnormcountsByPatient,rowMeans))
  }
  invisible(gc())
  LnormcountsByPatient <- lapply(LnormcountsByPatient, function(d) d-SF)
  invisible(gc())
  LnormcountsByPatient <- lapply(LnormcountsByPatient, function(d) funTrim(d, lo=-3, up=3))
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

  LgenomicScores <- lapply(LnormcountsByPatient,funGenomicScore,"GW"=genomicWindows)
  rm(LnormcountsByPatient) ; invisible(gc())


  if(scaleOnReferenceLabel){
    if (length(referenceLabel) == 1){
      genomicScoresReferenceLabel <- do.call(rbind, lapply(names(LN), function(patient) LgenomicScores[[patient]][LN[[patient]],]))
    } else {
      genomicScoresReferenceLabel <- list()
      for (i in referenceLabel){
        genomicScoresReferenceLabel[[i]] <- do.call(rbind, lapply(names(LN[[i]]), function(patient) LgenomicScores[[patient]][LN[[i]][[patient]],]))
      }
      genomicScoresReferenceLabel <- do.call(rbind, genomicScoresReferenceLabel)
    }
    Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, "probs"=c(0+thresholdPercentile,1-thresholdPercentile))
    LgenomicScoresTrimmed <- lapply(LgenomicScores, function(gs){
      apply(gs, 1, function(v) { v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v }) })

  } else {
    genomicScoresAll <- do.call(rbind, lapply(names(LgenomicScores), function(patient) LgenomicScores[[patient]]))
    Q01Q99 <- apply(genomicScoresAll,2,quantile,"probs"=c(0+thresholdPercentile,1-thresholdPercentile))
    LgenomicScoresTrimmed <- lapply(LgenomicScores,function(gs){
      apply(gs,1,function(v) { v[which(v >= Q01Q99[1,] & v <= Q01Q99[2,])] <- 0 ; v }) })
  }

  LgenomicScores <- lapply(LgenomicScores, function(x) t(x))


  if (saveGenomicWindows){
    save(genomicWindows, file = paste0("genomicWindows_size",windowSize,"_step",windowStep,".RData"))
  }

  for (i in 1:length(seuratList)) {
    rawGenomicAssay <- Seurat::CreateAssayObject(counts = LgenomicScores[i][[1]])
    suppressWarnings({seuratList[[i]][["rawGenomicScores"]] <- rawGenomicAssay})
    genomicAssay <- Seurat::CreateAssayObject(counts = LgenomicScoresTrimmed[i][[1]])
    suppressWarnings({seuratList[[i]][["genomicScores"]] <- genomicAssay})
    seuratList[[i]][["cnv_fraction"]] <- colMeans(abs(LgenomicScoresTrimmed[i][[1]]) > 0)
  }
  invisible(gc())
  return (seuratList)
}
