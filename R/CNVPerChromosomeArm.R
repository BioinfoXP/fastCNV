#' CNV Per Chromosome Arm
#' Computes the CNV fraction of each spot/cell per chromosome arm, then stores the results into the metadata.
#'
#' @param seuratObj A Seurat object, typically the output from the `fastCNV()` function, containing genomic scores for CNV analysis.
#'
#' @return The function returns the same Seurat object with the CNV fraction for each chromosome arm added to the metadata.
#'
#'
#' @export


CNVPerChromosomeArm <- function(seuratObj) {
  genomicScores <- as.matrix(Seurat::GetAssay(seuratObj, assay = "genomicScores")["data"])
  window_names <- rownames(genomicScores)

  extract_chrom_arm <- function(window_name) {
    chrom_arm <- sub("(\\d+\\.\\w).*", "\\1", window_name)
    chrom_arm <- sub("^23\\.(p|q)$", "X.\\1", chrom_arm)
    return(chrom_arm)
  }

  window_info <- data.frame(
    window = window_names,
    chrom_arm = sapply(window_names, extract_chrom_arm)
  )

  chrom_arms_all <- unique(c(
    paste0(1:22, ".p"), paste0(1:22, ".q"),
    "X.p", "X.q"
  ))

  arm_averages <- list()
  for (chrom_arm in chrom_arms_all) {
    if (chrom_arm %in% unique(window_info$chrom_arm)) {
      windows <- window_info$window[window_info$chrom_arm == chrom_arm]
      subset_scores <- genomicScores[windows, , drop = FALSE]
      avg_values <- colMeans(subset_scores, na.rm = TRUE)
    } else {
      avg_values <- rep(0, ncol(genomicScores))
    }
    arm_averages[[chrom_arm]] <- avg_values
  }

  arm_averages_df <- do.call(rbind, lapply(names(arm_averages), function(chrom_arm) {
    data.frame(
      chrom_arm = chrom_arm,
      value = arm_averages[[chrom_arm]],
      barcode = colnames(genomicScores),
      check.names = FALSE
    )
  }))


  meta <- seuratObj@meta.data
  meta$barcode <- rownames(meta)

  for (chrom_arm in unique(arm_averages_df$chrom_arm)) {
    arm_df <- arm_averages_df[arm_averages_df$chrom_arm == chrom_arm, ]
    lookup <- setNames(arm_df$value, arm_df$barcode)
    new_col <- paste0(chrom_arm, "_CNV")
    meta[[new_col]] <- NA
    idx <- meta$barcode %in% arm_df$barcode
    meta[[new_col]][idx] <- lookup[meta$barcode[idx]]
  }

  meta$barcode <- NULL
  seuratObj@meta.data <- meta
  return(seuratObj)
}
