#' printMedPathCoeff
#'
#' @description Prints out mediation path coefficients in a formatted table
#'
#' @param pathCoeff List of mediation path coefficients for each PDM
#'
#' @importFrom pracma fprintf
#'
#' @return A printed, formatted table of the mediation path coefficient values
#' for a, b, ab, c' for each calculated PDM
#'
#' @export


printMedPathCoeff <- function(pathCoeff){

  if ( 'redFeatWeights' %in% names(pathCoeff)){
    stop('This function only accepts $pathCoeff list.
         Do not input entire output list from getDirectionsofMed()')
  }

  # Number of PDMs calculated
  num <- length(pathCoeff)

  # Unlist pathCoeff
  pathCoeff <- matrix(unlist(pathCoeff), nrow = 5, ncol = num)

  # Print table
  dashes <- '_____________________________________________________________________'
  cat(sprintf("\nPDM path coefficients \n%s \n \tpath a\t\tpath b\t\tpath ab\t\tpath c'\n",dashes))
  for (k in 1:num){
    pracma::fprintf('PDM%d \t%0.4f\t\t%0.4f\t\t%0.4f\t\t%0.4f\n',k,
                    pathCoeff[3,k],pathCoeff[4,k],pathCoeff[5,k],pathCoeff[2,k])
  }
  cat(sprintf('%s\n',dashes))
}
