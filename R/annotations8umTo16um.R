#' Project 8µm Spatial Annotation onto 16µm Spots
#'
#' This function projects annotations from a high-resolution (8µm) spatial assay onto a lower-resolution (16µm) spatial assay
#' by finding the nearest 8µm spot to each 16µm spot based on spatial coordinates.
#'
#' @param HDobj A \code{Seurat} object containing both 8µm and 16µm spatial assays (named \code{Spatial.008um} and \code{Spatial.016um}).
#' @param referenceVar A character string specifying the name of the metadata column in the 8µm assay to project (e.g., a clustering or annotation label).
#'
#' @return A modified \code{Seurat} object with a new metadata column named \code{projected_<referenceVar>} containing the projected annotation on 16µm spots.
#'
#' @details
#' The function uses \code{FNN::get.knnx()} to find the nearest 8µm spot for each 16µm spot based on tissue coordinates.
#' It assigns the annotation from the closest 8µm spot to each 16µm spot. The new annotation column is added to the metadata of \code{HDobj}.
#'
#' @import Seurat
#' @importFrom FNN get.knnx
#' @importFrom stats setNames
#'
#'
#' @export


annotations8umTo16um <- function(HDobj,
                                    referenceVar){

  Seurat::DefaultAssay(HDobj) = "Spatial.008um"
  visium8um <- suppressWarnings(suppressMessages(subset(HDobj, cells = Seurat::Cells(HDobj))))
  invisible(gc())

  Seurat::DefaultAssay(HDobj) = "Spatial.016um"
  visium16um <- suppressWarnings(suppressMessages(subset(HDobj, cells = Seurat::Cells(HDobj))))
  invisible(gc())

  coords_8um <- Seurat::GetTissueCoordinates(visium8um)
  coords_16um <- Seurat::GetTissueCoordinates(visium16um)

  mat_8um <- as.matrix(coords_8um[, c("x", "y")])
  mat_16um <- as.matrix(coords_16um[, c("x", "y")])

  nn <- FNN::get.knnx(mat_8um, mat_16um, k = 1)

  nearest_8um_index <- nn$nn.index[, 1]

  annotations_8um <- Seurat::FetchData(visium8um, vars = referenceVar)[,1]
  annotations_8um <- setNames(visium8um@meta.data[[referenceVar]], Seurat::Cells(visium8um))


  nearest_8um_names <- rownames(coords_8um)[nearest_8um_index]

  projected_annotations <- annotations_8um[nearest_8um_names]

  names(projected_annotations) <- rownames(coords_16um)

  visium16um$projected_annotation <- projected_annotations


  HDobj[[paste0("projected_",referenceVar)]] <- NA
  HDobj@meta.data[Seurat::Cells(visium16um),paste0("projected_",referenceVar)] = visium16um$projected_annotation
  rm(coords_16um,coords_8um,mat_16um,mat_8um,nn,annotations_8um,nearest_8um_index,nearest_8um_names,projected_annotations,visium16um,visium8um) ; invisible(gc())

  print(paste0("New annotation column (new referenceVar) is named : projected_",referenceVar,"."))

  return (HDobj)
}
