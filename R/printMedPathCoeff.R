#' printMedPathCoeff
#'
#' @description Prints out mediation path coefficients in a formatted table
#'
#' @param theta List of mediation path coefficients for each PDM
#'
#' @return A printed, formatted table of the mediation path coefficient values
#' for a, b, ab, c' for each calculated PDM
#'
#' @export


printMedPathCoeff <- function(theta){

  if ('w_k' %in% names(theta)){
    stop('This function only accepts $Theta list.
         Do not input entire output list from getDirectionsofMed()')
  }

  # Number of PDMs calculated
  num <- length(theta)

  # Unlist theta
  theta <- matrix(unlist(theta), nrow = 5, ncol = num)

  # Print table
  dashes <- '_____________________________________________________________________'
  cat(sprintf("\nPDM path coefficients \n%s \n \tpath a\t\tpath b\t\tpath ab\t\tpath c'\n",dashes))
  for (k in 1:num){
    pracma::fprintf('PDM%d \t%0.4f\t\t%0.4f\t\t%0.4f\t\t%0.4f\n',k,
                    theta[3,k],theta[4,k],theta[5,k],theta[2,k])
  }
  cat(sprintf('%s\n',dashes))
}
