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
    if (!is.list(object)) {
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
      if (length(object) == 1) {
        object <- list(CNVCalling(object[1],
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
        if (pooledReference == TRUE) {
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

        # if (doRecapPlot == TRUE) {
        #   message("Plotting the CNV recap heatmap per category accross samples. This could take some time.")
        #   if (!is.null(referenceVar)) {
        #   palette <- as.character(paletteer::paletteer_d("ggsci::default_igv"))
        #   annot_colors <- setNames(palette[1:length(names(object))], names(object))
        #   uni <- list()
        #   for (i in names(object)){
        #     uni[[i]] <- unique(object[[i]]@meta.data[[referenceVar]])
        #   }
        #   cellTypes <- unique(unlist(uni))
        #   cellTypes <- cellTypes[nzchar(cellTypes)]
        #
        #   Lmat <- list()
        #   for (i in cellTypes){
        #     Lmat[[i]] <- lapply(object, function(x)
        #       Seurat::Cells(x)[which(x@meta.data[[referenceVar]] == i)])
        #     Lmat[[i]] <- Lmat[[i]][which(sapply(Lmat[[i]], length)>=5)]
        #     for (j in names(Lmat[[i]])) {
        #       Lmat[[i]][[j]] <- as.matrix(Seurat::GetAssay(object[[j]], assay = "genomicScores")["counts"])[,Lmat[[i]][[j]]]
        #     }
        #   }
        #   Lhm <- list()
        #   for (x in names(Lmat)) {
        #     modified <- list()
        #     sample_info <- data.frame(col_name = character(), sample = character(), stringsAsFactors = FALSE)
        #     for (sample in names(Lmat[[x]])) {
        #       matrix <- Lmat[[x]][[sample]]
        #       col_names <- colnames(matrix)
        #       new_col_names <- paste0(sample, "_", col_names)
        #       colnames(matrix) <- new_col_names
        #       modified[[sample]] <- rowMeans(matrix)
        #
        #       sample_info <- rbind(sample_info, data.frame(col_name = new_col_names, sample = sample))
        #     }
        #     rownames(sample_info) <- sample_info[,1]
        #     sample_info[,1] <- NULL
        #     combined <- t(do.call(cbind, modified))
        #
        #     sample_df <- as.data.frame(unique(sample_info))
        #     rownames(sample_df) <- sample_df[,1]
        #
        #     Lhm[[x]] <- ComplexHeatmap::pheatmap(combined,
        #                                          border=T,
        #                                          border_color = NA,
        #                                          use_raster = F,
        #                                          cluster_cols = F,
        #                                          show_rownames = F,
        #                                          show_colnames = F,
        #                                          clustering_distance_rows = "euclidean",
        #                                          clustering_method = "ward.D",
        #                                          column_split = as.numeric(sapply(strsplit(colnames(combined),".",fixed=T),function(z)z[1])),
        #                                          row_split = sample_df,
        #                                          annotation_row = sample_df,
        #                                          annotation_colors = list(sample = annot_colors),
        #                                          col=circlize::colorRamp2(c(-1, 0, 1), c("darkblue", "white", "darkred")))
        #
        #   }
        #   output_file <- "heatmaps_per_annotation.pdf"
        #
        #   pdf(output_file, width = 20, height = 10)
        #   for (i in seq_along(Lhm)){
        #     grid::grid.newpage()
        #     grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit.c(grid::unit(1, "cm"),grid::unit(1, "null")))))
        #     grid::pushViewport(grid::viewport(layout.pos.row = 1))
        #     grid::grid.text(paste0("CNV heatmap for annotation : ", names(Lhm)[i]), gp = grid::gpar(fontsize = 20))
        #     grid::popViewport()
        #     grid::pushViewport(grid::viewport(layout.pos.row = 2))
        #     ComplexHeatmap::draw(Lhm[[i]], newpage = FALSE)
        #     grid::popViewport()
        #
        #   }
        #   dev.off()
        # }
        # }

      }
    }
  message(crayon::green(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Done !")))
  return (object)
}
