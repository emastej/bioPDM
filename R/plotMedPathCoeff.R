#' plotMedPathCoeff
#'
#' @description Plots mediation path coefficients of a, b, and |ab|
#'
#' @param theta List of mediation path coefficients for each PDM
#'
#' @importFrom graphics par
#'
#' @return A figure with 3 subplots showing how Path A (red), Path B (blue),
#' and |Path AB| (green) change in value for each PDM
#'
#' @export
#'


plotMedPathCoeff <- function(theta){

  if ('w_k' %in% names(theta)){
    stop('This function only accepts $Theta list.
         Do not input entire output list from getDirectionsofMed()')
  }

  # Number of PDMs
  num <- length(theta)

  # Unlist theta
  theta <- matrix(unlist(theta), nrow = 5, ncol = num)

  # Subplot
  graphics::par(mfrow=c(1,3))

  # Plot Path A
  plot(theta[3,],
       type = 'o',
       col = 'red',
       lwd = 4,
       main = 'path a',
       xlab = 'PDM #',
       ylab = 'Coefficients',
       cex.lab = 1.5,
       cex.main = 2,
       cex.axis = 1.5)

  # Plot Path B
  plot(theta[4,],
       type = 'o',
       col = 'blue',
       lwd = 4,
       main = 'path b',
       xlab = 'PDM #',
       ylab = 'Coefficients',
       cex.lab = 1.5,
       cex.main = 2,
       cex.axis = 1.5)

  # Plot Path AB
  plot(abs(theta[5,]),
       type = 'o',
       col = 'green',
       lwd = 4,
       main = '|path ab|',
       xlab = 'PDM #',
       ylab = 'Coefficients',
       cex.lab = 1.5,
       cex.main = 2,
       cex.axis = 1.5)

  # Reset the figure subplot configuration
  graphics::par(mfrow=c(1,1))
}
