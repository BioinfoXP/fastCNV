#' Genes Data from Ensembl Version 113
#'
#' Data downloaded from the Ensembl website (version 113), containing detailed gene information for approximately 76,000 genes.
#' The dataset includes Ensembl gene IDs, HUGO nomenclature (HGNC symbol), Entrez gene IDs, chromosome locations,
#' gene biotype, and gene length for each gene.
#'
#' @docType data
#'
#' @usage data(geneMetadata)
#'
#' @format An object of class `list`, containing gene information as described above.
#'
#' @keywords datasets
#'
#' @source Ensembl Genome Browser, Version 113: https://www.ensembl.org/index.html
#'
#' @examples
#' data(geneMetadata)
#' hgnc <- geneMetadata$hgnc_symbol
#' entrez <- geneMetadata$entrezgene_id
"geneMetadata"
