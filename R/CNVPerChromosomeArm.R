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
  # --- 1. Safe Data Access (v4/v5 Compatible) ---
  genomicScores <- tryCatch({
    as.matrix(Seurat::GetAssayData(seuratObj, assay = "genomicScores", slot = "data"))
  }, error = function(e) {
    as.matrix(Seurat::GetAssayData(seuratObj, assay = "genomicScores", layer = "data"))
  })

  window_names <- rownames(genomicScores)

  extract_chrom_arm <- function(window_name) {
    chrom_arm <- sub("(\\d+\\.\\w).*", "\\1", window_name)
    chrom_arm <- sub("^23\\.(p|q)$", "X.\\1", chrom_arm)
    return(chrom_arm)
  }

  window_info <- data.frame(
    window = window_names,
    chrom_arm = sapply(window_names, extract_chrom_arm),
    stringsAsFactors = FALSE
  )

  # --- 2. Define Standard Arms (46 total) ---
  chrom_arms_standard <- c(
    paste0(rep(1:22, each=2), c(".p", ".q")),
    "X.p", "X.q"
  )

  # --- 3. Compute Averages (Fill missing with 0) ---
  arm_averages <- list()
  present_arms <- unique(window_info$chrom_arm)

  for (chrom_arm in chrom_arms_standard) {
    if (chrom_arm %in% present_arms) {
      windows <- window_info$window[window_info$chrom_arm == chrom_arm]
      valid_windows <- intersect(windows, rownames(genomicScores))

      if(length(valid_windows) > 0) {
        subset_scores <- genomicScores[valid_windows, , drop = FALSE]
        arm_averages[[chrom_arm]] <- colMeans(subset_scores, na.rm = TRUE)
      } else {
        arm_averages[[chrom_arm]] <- rep(0, ncol(genomicScores))
      }
    } else {
      # CRITICAL: Fill missing arms with 0 to maintain column structure
      arm_averages[[chrom_arm]] <- rep(0, ncol(genomicScores))
    }
  }

  # --- 4. Write to Metadata ---
  meta <- seuratObj@meta.data

  # Clean up old columns to prevent duplication issues
  old_cols <- grep("_CNV$", colnames(meta), value = TRUE)
  if(length(old_cols) > 0) {
    meta <- meta[, !colnames(meta) %in% old_cols]
  }

  # Add columns in standard order
  for (chrom_arm in chrom_arms_standard) {
    col_name <- paste0(chrom_arm, "_CNV")
    meta[[col_name]] <- arm_averages[[chrom_arm]]
  }

  seuratObj@meta.data <- meta
  return(seuratObj)
}
