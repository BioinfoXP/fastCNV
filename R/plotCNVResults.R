#' Plot CNV Results into a Heatmap
#' Builds a heatmap to visualize the CNV results based on genomic scores.
#'
#' @param seuratObj A Seurat object containing the genomic scores computed previously.
#' @param referenceVar The name of the metadata column in the Seurat object containing reference annotations.
#' @param clustersVar The name of the metadata column containing cluster information (default = `"cnv_clusters"`).
#' @param splitPlotOnVar The name of the metadata column used to split the heatmap rows (e.g., cell type or cluster) (default = `clustersVar`).
#' @param denoise If `TRUE`, the denoised data will be used in the heatmap (default = `TRUE`).
#' @param savePath The path where the heatmap will be saved. If `NULL`, the plot will not be saved (default = `"."`).
#' @param printPlot Logical. If `TRUE`, prints the heatmap to the console.
#' @param referencePalette A color palette for `referenceVar`.
#' You can provide a custom palette as a vector of color codes (e.g., `c("#FF0000", "#00FF00")`).
#' @param clusters_palette A color palette for `clustersVar`.
#' You can provide a custom palette as a vector of color codes (e.g., `c("#F8766D", "#A3A500", "#00BF7D")`).
#' @param outputType Character. Specifies the file format for saving the plot, either `"png"` or `"pdf"`.
#'
#' @importFrom Seurat GetAssay FetchData
#' @importFrom ComplexHeatmap Heatmap rowAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit grid.newpage pushViewport viewport grid.layout grid.text popViewport
#' @importFrom paletteer paletteer_d
#' @importFrom scales hue_pal
#' @importFrom crayon black
#' @import stats
#' @import grDevices
#'
#' @return This function generates a heatmap and saves it as a `.pdf` or `.png` file in the specified path (default = working directory).
#'
#' @export


plotCNVResults <- function(seuratObj,
                           referenceVar = NULL,
                           clustersVar = "cnv_clusters",
                           splitPlotOnVar = clustersVar,
                           denoise = TRUE,
                           savePath = ".",
                           printPlot = FALSE,
                           referencePalette = "default",
                           clusters_palette = "default",
                           outputType = "png"){
  if (outputType != "png" && outputType != "pdf"){
    message("Warning : outputType not valid, should be 'pdf' or 'png'. Setting outputType to 'png'")
    outputType = "png"
  }

  if (!is.null(clustersVar)){
    if (clustersVar == "cnv_clusters"){
      if (!(clustersVar %in% names(seuratObj@meta.data))) {
        if(!is.null(splitPlotOnVar)) {
          if(splitPlotOnVar == clustersVar){
            splitPlotOnVar = referenceVar
          }
        }
        clustersVar = NULL
      }
    }
  }

  if (denoise == TRUE) {
    mat <- as.matrix(Seurat::GetAssay(seuratObj, "genomicScores")["data"])
    arms <- rownames(mat)

    # define grouping variable (e.g. 1.p, 1.q, 13.q, etc.)
    arms_group <- stringr::str_extract(arms, "^[0-9XY]+\\.[pq]")

    # unique arms in the order they appear
    unique_arms <- unique(arms_group)

    # build new matrix
    M <- do.call(rbind, lapply(unique_arms, function(arm) {
      idx <- which(arms_group == arm)
      rows <- mat[idx, , drop = FALSE]

      # duplicate if fewer than 3
      if (nrow(rows) < 3) {
        need <- 3 - nrow(rows)
        dup_rows <- rows[rep(1:nrow(rows), length.out = need), , drop = FALSE]

        # give unique rownames to avoid collisions
        rownames(dup_rows) <- paste0(rownames(rows)[rep(1:nrow(rows), length.out = need)], "_dup")

        rows <- rbind(rows, dup_rows)
      }
      rows
    }))
    M <- t(M)

  } else {
    mat <- as.matrix(Seurat::GetAssay(seuratObj, "rawGenomicScores")["data"])
    arms <- rownames(mat)

    # define grouping variable (e.g. 1.p, 1.q, 13.q, etc.)
    arms_group <- stringr::str_extract(arms, "^[0-9XY]+\\.[pq]")

    # unique arms in the order they appear
    unique_arms <- unique(arms_group)

    # build new matrix
    M <- do.call(rbind, lapply(unique_arms, function(arm) {
      idx <- which(arms_group == arm)
      rows <- mat[idx, , drop = FALSE]

      # duplicate if fewer than 3
      if (nrow(rows) < 3) {
        need <- 3 - nrow(rows)
        dup_rows <- rows[rep(1:nrow(rows), length.out = need), , drop = FALSE]

        # give unique rownames to avoid collisions
        rownames(dup_rows) <- paste0(rownames(rows)[rep(1:nrow(rows), length.out = need)], "_dup")

        rows <- rbind(rows, dup_rows)
      }
      rows
    }))
    M <- t(M)
  }
  if (any(referencePalette == "default")) {
    referencePalette = as.character(paletteer::paletteer_d("pals::glasbey"))
  }
  if (!is.null(referenceVar)) {
    annotation_df <- as.data.frame(seuratObj@meta.data[[referenceVar]])
    colnames(annotation_df) <- "Annotations"
    annot_colors <- setNames(referencePalette[1:length(unique(annotation_df$Annotations))], unique(annotation_df$Annotations))
    if (!is.null(splitPlotOnVar)) {
      split_df <- as.data.frame(Seurat::FetchData(seuratObj, vars = splitPlotOnVar))
      colnames(split_df) <- "Split"
    } else {
      split_df <- NULL
    }
  } else if (is.null(splitPlotOnVar)) {
    split_df <- NULL
  }
  if (!is.null(clustersVar)) {
    clusters_df <- as.data.frame(seuratObj@meta.data[[clustersVar]])
    colnames(clusters_df) <- "Clusters"
    if (any(clusters_palette == "default")){
      clusters_palette = scales::hue_pal()(length(unique(clusters_df$Clusters)))
    }
    clusters_colors <- setNames(clusters_palette, sort(unique(clusters_df$Clusters)))
  }

  if (!is.null(referenceVar) && is.null(clustersVar)) {
    annotation_heatmap <- ComplexHeatmap::rowAnnotation(
      Annotations = annotation_df$Annotations,
      col = list(Annotations = annot_colors),
      annotation_legend_param = list(
        title = "Annotations",
        title_gp = grid::gpar(fontsize = 11),
        labels_gp = grid::gpar(fontsize = 8),
        legend_height = grid::unit(3, "cm"),
        legend_width = grid::unit(1.5, "cm"),
        grid_height = grid::unit(0.6, "cm"),
        grid_width = grid::unit(0.6, "cm")
      )
    )
  }

  if (is.null(referenceVar) && !is.null(clustersVar)) {
    annotation_heatmap <- ComplexHeatmap::rowAnnotation(
      Clusters = clusters_df$Clusters,
      col = list(Clusters = clusters_colors),
      annotation_legend_param = list(
        title = "Clusters",
        title_gp = grid::gpar(fontsize = 11),
        labels_gp = grid::gpar(fontsize = 8),
        legend_height = grid::unit(3, "cm"),
        legend_width = grid::unit(1.5, "cm"),
        grid_height = grid::unit(0.6, "cm"),
        grid_width = grid::unit(0.6, "cm")
      )
    )
    if (!is.null(splitPlotOnVar)) {
      split_df <- as.data.frame(Seurat::FetchData(seuratObj, vars = splitPlotOnVar))
      colnames(split_df) <- "Split"
    } else {
      split_df <- NULL
    }
  }

  if (is.null(referenceVar) && is.null(clustersVar)) {
    annotation_heatmap <- NULL
  }

  if (!is.null(referenceVar) && !is.null(clustersVar)) {
    annotation_heatmap <- ComplexHeatmap::rowAnnotation(
      Annotations = annotation_df$Annotations,
      Clusters = clusters_df$Clusters,
      col = list(
        Annotations = annot_colors,
        Clusters = clusters_colors
      ),
      annotation_legend_param = list(
        Annotations = list(
          title = "Annotations",
          title_gp = grid::gpar(fontsize = 11),
          labels_gp = grid::gpar(fontsize = 8)
        ),
        Clusters = list(
          title = "CNV cluster",
          title_gp = grid::gpar(fontsize = 11),
          labels_gp = grid::gpar(fontsize = 8)
        )
      )
    )
  }

  if (is.null(split_df)) {
    splitting = NULL
  } else {
    splitting = as.factor(split_df[[1]])
  }

  hm <-  ComplexHeatmap::Heatmap(
    M,
    right_annotation = annotation_heatmap,
    border = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D",
    column_split = as.numeric(sapply(strsplit(colnames(M), ".", fixed = TRUE), function(z) z[1])),
    column_title_gp = grid::gpar(fontsize = 8),
    column_title = c(1:22,"X"),
    row_split = splitting,
    row_title = NULL,
    col = circlize::colorRamp2(c(-1, -0.6, -0.3, 0, 0.3, 0.6, 1), c("#0B2F7EFF", "#2A4D9EFF", "#A0A0FFFF", "white", "#E3807D", "#A4161A","#7A0A0D")),
    heatmap_legend_param = list(
      title = "CNV Score",
      title_gp = grid::gpar(fontsize = 11),
      labels_gp = grid::gpar(fontsize = 8),
      grid_height = grid::unit(1, "cm"),
      grid_width = grid::unit(0.6, "cm"))
    )

  if(printPlot == TRUE) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit.c(grid::unit(1, "cm"), grid::unit(1, "null")))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.text(paste0("CNV heatmap for sample ", seuratObj@project.name), gp = grid::gpar(fontsize = 20))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    ComplexHeatmap::draw(hm, newpage = FALSE)
    grid::popViewport()
  }

  if(!is.null(savePath)) {
    if (outputType == "png") {
      fname <- file.path(savePath, paste0("heatmap.fastCNV_",seuratObj@project.name,".png"))
      png(width = 3500, height = 2200, filename = fname, res = 300)
    }
    if (outputType == "pdf"){
      fname <- file.path(savePath, paste0("heatmap.fastCNV_",seuratObj@project.name,".pdf"))
      pdf(width = 12, height = 7, file = fname)
    }
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit.c(grid::unit(1, "cm"),grid::unit(1, "null")))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1,gp = grid::gpar(fill = "white")))
    grid::grid.text(paste0("CNV heatmap for sample ", seuratObj@project.name), gp = grid::gpar(fontsize = 20))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    ComplexHeatmap::draw(hm, newpage = FALSE)
    grid::popViewport()
    Sys.sleep(3)
    dev.off()
    message(crayon::black,"CNV plot for sample ",seuratObj@project.name, " saved at ", fname)
  }
  return(hm)
}
