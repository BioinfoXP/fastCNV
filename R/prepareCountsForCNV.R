#' Aggregate Observations by Cell Type for CNV Analysis
#' Aggregates observations with the same cell types to increase counts per observation,
#' improving Copy Number Variation (CNV) computation.
#'
#' @param seuratObj A Seurat object containing the data.
#' @param sampleName A character string specifying the sample name.
#' @param referenceVar The name of the metadata column in the Seurat object that contains reference annotations.
#' @param aggregateByVar Logical. If `TRUE` (default), aggregates observations based on `referenceVar` annotations.
#' @param aggregFactor Integer. The target number of counts per observation (default = `15 000`).
#' @param seuratClusterResolution Numeric. The resolution used for Seurat clustering (default = `0.8`).
#' @param reClusterSeurat Logical. If `TRUE`, re-runs clustering on the Seurat object.
#'
#' @return A Seurat object with:
#' - A new assay called **"AggregatedCounts"** containing the modified count matrix.
#' - **Seurat clusters** stored in the metadata.
#'
#' @import Seurat
#' @importFrom crayon black
#'
#' @export
#'
prepareCountsForCNVAnalysis <- function(seuratObj,
                                        sampleName = NULL,
                                        referenceVar = NULL,
                                        aggregateByVar = T,
                                        aggregFactor=15000,
                                        seuratClusterResolution = 0.8,
                                        reClusterSeurat = F ){

    assay <- Seurat::Assays(seuratObj)[1]
      if (is.null(referenceVar) || aggregateByVar == F) {
        if (!"seurat_clusters" %in% colnames(seuratObj[[]]) | reClusterSeurat) {
          if (!is.null(sampleName)){message(crayon::black,paste0("Running Seurat SCTransform and clustering for sample ", sampleName,". This could take some time."))
          }else{message(crayon::black,"Running Seurat SCTransform and clustering. This could take some time.")}
          seuratObj <- Seurat::SCTransform(seuratObj, assay = assay, verbose = F)
          seuratObj <- Seurat::RunPCA(seuratObj, assay = "SCT", verbose = F)
          seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = "pca", dims = 1:10, verbose = F)
          seuratObj <- Seurat::FindClusters(seuratObj, resolution = seuratClusterResolution, verbose = F)
          if (!is.null(sampleName)){message(crayon::black,paste0("Seurat SCTransform and clustering done for sample ", sampleName))}
          else {message(crayon::black,"Seurat SCTransform and clustering done.")}
        }

        countsMat <- as.matrix(Seurat::GetAssay(seuratObj, assay = assay)$counts)
        LC = split(Seurat::Cells(seuratObj), Seurat::FetchData(seuratObj, vars = "seurat_clusters"))
        LC = lapply(lapply(LC, function(cells) {nc <- Seurat::FetchData(seuratObj, vars = paste0("nCount_",assay))[cells,] ; names(nc) <- cells ; nc}), sort)

        seuratObj <- Seurat::AddMetaData(seuratObj, metadata = 0, col.name = "metaSpots")

        nbMetaSpots <- 1

        for (X in LC) {
          if (length(X) != 0) {
            barcodes = c()
            taille = 0

            for(i in 1:length(X)) {
              taille = taille + as.numeric(X[i])
              barcodes = append(barcodes, names(X[i]))

              if (taille >= aggregFactor) {
                if (length(barcodes)>1) {
                  countsMat[,barcodes] <- rowSums(countsMat[,barcodes])
                }
                seuratObj@meta.data[barcodes,"metaSpots"] <- nbMetaSpots
                nbMetaSpots = nbMetaSpots + 1
                taille = 0
                barcodes = c()
              }
              if (i == length(X)){
                if (length(barcodes)>1) {
                  countsMat[,barcodes] <- rowSums(countsMat[,barcodes])
                }
                seuratObj@meta.data[barcodes,"metaSpots"] <- nbMetaSpots
                nbMetaSpots = nbMetaSpots + 1
              }
            }
          }
        }
        aggregAssay <- Seurat::CreateAssayObject(counts = countsMat)
        seuratObj[["AggregatedCounts"]] <- aggregAssay

      } else {

        if (!"seurat_clusters" %in% colnames(seuratObj[[]]) | reClusterSeurat) {
          if (!is.null(sampleName)){message(crayon::black,paste0("Running Seurat SCTransform and clustering for sample ", sampleName,". This could take some time."))
          }else{message(crayon::black,"Running Seurat SCTransform and clustering. This could take some time.")}
          seuratObj <- Seurat::SCTransform(seuratObj, assay = assay, verbose = F)
          seuratObj <- Seurat::RunPCA(seuratObj, assay = "SCT", verbose = F)
          seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = "pca", dims = 1:10, verbose = F)
          seuratObj <- Seurat::FindClusters(seuratObj, resolution = seuratClusterResolution, verbose = F)
          if (!is.null(sampleName)){message(crayon::black,paste0("Seurat SCTransform and clustering done for sample ", sampleName,"."))}
          else {message(crayon::black,"Seurat SCTransform and clustering done.")}
        }

        countsMat <- as.matrix(Seurat::GetAssay(seuratObj, assay = Seurat::Assays(seuratObj)[1])$counts)
        spotsCateg <- split(Seurat::Cells(seuratObj),Seurat::FetchData(seuratObj, vars = referenceVar))
        LC <- lapply(spotsCateg, function(x) split(x, Seurat::FetchData(seuratObj, vars = "seurat_clusters", cells = x)))

        remove_empty_lists <- function(list_of_lists) {
          lapply(list_of_lists, function(sub_list) {
            Filter(function(x) length(x) > 0, sub_list)
          })
        }

        LC <- remove_empty_lists(LC)
        LC = lapply(LC, function(x) {lapply(
                lapply(x, function(cells) {
                  nc <- Seurat::FetchData(seuratObj, vars = paste0("nCount_",assay))[cells,]
                  names(nc) <- cells
                  nc}
                ), sort)})

        seuratObj <- Seurat::AddMetaData(seuratObj, metadata = 0, col.name = "metaSpots")

        nbMetaSpots <- 1

        for (i in LC) {
          for (X in i) {
            if (length(X) != 0) {
              barcodes = c()
              taille = 0

              for(i in 1:length(X)) {
                taille = taille + as.numeric(X[i])
                barcodes = append(barcodes, names(X[i]))

                if (taille >= aggregFactor) {
                  if (length(barcodes)>1) {
                    countsMat[,barcodes] <- rowSums(countsMat[,barcodes])
                  }
                  seuratObj@meta.data[barcodes,"metaSpots"] <- nbMetaSpots
                  nbMetaSpots = nbMetaSpots + 1
                  taille = 0
                  barcodes = c()
                }
                if (i == length(X)){
                  if (length(barcodes)>1) {
                    countsMat[,barcodes] <- rowSums(countsMat[,barcodes])
                  }
                  seuratObj@meta.data[barcodes,"metaSpots"] <- nbMetaSpots
                  nbMetaSpots = nbMetaSpots + 1
                }
              }
            }
          }
        }
        aggregAssay <- Seurat::CreateAssayObject(counts = countsMat)
        seuratObj[["AggregatedCounts"]] <- aggregAssay

      }
  invisible(gc())

   if (!is.null(sampleName)) {
    if (Seurat::Project(seuratObj) == "SeuratProject") {
       Seurat::Project(seuratObj) <- sampleName
    }}
  return(seuratObj)
}
