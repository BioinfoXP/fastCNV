#' CNV Classification for a Given Vector
#' Classifies each CNV value in a vector as "loss", "gain", or "no alteration" based on specified thresholds.
#'
#' @param cnv_vector A numeric vector of CNV scores for a given chromosome arm.
#' @param peaks A numeric vector containing the thresholds for classification.
#'   The default is `c(-0.1, 0, 0.1)`:
#'   - Loss: CNV values below the first threshold (`peaks[1]`).
#'   - Gain: CNV values above the third threshold (`peaks[3]`).
#'   - No alteration: CNV values between the first and third thresholds.
#'
#' @return A character vector of classifications for each CNV value in the input vector,
#'   with possible values of `"loss"`, `"gain"`, or `"no_alteration"`.
#'
#' @keywords internal
#'
#' @export

classifyCNV <- function(cnv_vector, peaks) {
  cnv_classification <- c()  # Initialize a vector to store the classification

  for (cnv in cnv_vector){
    if (cnv < peaks[1]) {
      classification = "loss"
    } else if (cnv > peaks[3]) {
      classification = "gain"
    } else {
      classification = "no_alteration"
    }
    cnv_classification <- c(cnv_classification, classification)
  }

  return(cnv_classification)
}
