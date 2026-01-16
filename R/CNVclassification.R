#' CNV Classification
#' Classifies the CNV results into loss, gain, or no alteration for each observation and chromosome arm.
#'
#' @param seuratObj A Seurat object containing the results of the CNV analysis (e.g., from `fastCNV`).
#' @param peaks A numeric vector containing the thresholds for classifying CNVs. The default is `c(-0.1, 0, 0.1)`, which defines:
#'   - Loss: CNV scores below `-0.1`
#'   - No alteration: CNV scores between `-0.1` and `0.1`
#'   - Gain: CNV scores above `0.1`
#'
#' @return The same Seurat object with an additional classification for each observation and chromosome arm in the metadata.
#'   The classification can be one of `"loss"`, `"gain"`, or `"no_alteration"`.
#'
#'
#' @export

CNVClassification <- function(seuratObj, peaks = c(-0.1,0,0.1)) {
  chromosomes <- c(paste0(1:22, ".p_CNV"), paste0(1:22, ".q_CNV"))
  chromosomes <- c(chromosomes, "X.p_CNV", "X.q_CNV")

  metadata <- Seurat::FetchData(seuratObj, vars = colnames(seuratObj@meta.data))


  for (chrom in chromosomes) {
    cnvVector <- metadata[[chrom]]
    classification <- classifyCNV(cnvVector, peaks)
    metadata[[paste0(chrom, "_classification")]] <- classification
  }

  seuratObj@meta.data <- metadata

  return (seuratObj)
}
