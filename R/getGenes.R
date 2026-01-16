#' Download Gene Information from Ensembl
#'
#' This function retrieves gene information from the Ensembl database using the specified filters.
#' It can either fetch the latest data or use cached data if available.
#'
#' @param filters A character vector of filters to be applied in the query. These filters determine
#' which genes and their associated information are returned from the Ensembl database.
#' @param cache Logical. If `TRUE`, the function will use cached data if available. If `FALSE`, it will
#' download the latest version of the gene data from Ensembl.
#'
#' @return A list containing gene information retrieved from Ensembl, with each element representing
#' data for a specific gene (e.g., gene IDs, descriptions, associated attributes).
#'
#' @importFrom biomaRt useMart useDataset getBM
#' @importFrom utils data
#'
#' @export

getGenes <- function(filters=NULL, cache=TRUE){
  if(cache) {
    message(crayon::black,"Using genes data from ensembl version 113.")
    data(geneMetadata,envir = environment())
  } else {
    message("Retrieving data from ensembl current version...")
    message("This could take some time.")
    ensembl = biomaRt::useMart("ensembl")
    ensembl = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)
    columnstokeep <- c('ensembl_gene_id',   'hgnc_symbol',   'entrezgene_id',"chromosome_name", 'start_position','end_position','gene_biotype')

    if(is.null(filters)){
      geneMetadata <- biomaRt::getBM(attributes=columnstokeep,values=T,mart = ensembl)
    }else{
      geneMetadata <- biomaRt::getBM(attributes=columnstokeep,filters = filters,values=T,mart = ensembl)
    }
    geneMetadata$length <- 1+geneMetadata$end_position - geneMetadata$start_position
    geneMetadata$length.kb=geneMetadata$length/1000

    centromere_positions <- list(
      "1"  = 121535434,
      "2"  = 93216878,
      "3"  = 91189767,
      "4"  = 50738735,
      "5"  = 48268472,
      "6"  = 60335031,
      "7"  = 58972616,
      "8"  = 45610351,
      "9"  = 49418102,
      "10" = 40000000,
      "11" = 53000000,
      "12" = 35000000,
      "13" = 17000000,
      "14" = 17000000,
      "15" = 19000000,
      "16" = 36000000,
      "17" = 25000000,
      "18" = 17000000,
      "19" = 26000000,
      "20" = 28000000,
      "21" = 13000000,
      "22" = 14000000,
      "X"  = 58000000,
      "Y"  = 10000000
    )

    determine_chr_arm <- function(chromosome, position, centromere_list) {
      if (chromosome %in% names(centromere_list)) {
        cent_pos <- centromere_list[[chromosome]]
        if (position < cent_pos) {
          return("p")
        } else {
          return("q")
        }
      } else {
        return(NA)
      }
    }

    geneMetadata <- geneMetadata |>
      dplyr::rowwise() |>
      dplyr::mutate(
        chr_arm = determine_chr_arm(.data$chromosome_name, .data$start_position, centromere_positions)
      ) |>
      dplyr::ungroup()
  }

  return(geneMetadata)
}
