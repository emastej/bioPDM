#' plotMedPathCoeff
#'
#' @description Plots mediation path coefficients of a, b, and |ab|
#'
#' @param pathCoeff List of mediation path coefficients for each PDM
#'
#' @importFrom graphics par
#'
#' @return A figure with 3 subplots showing how Path A (red), Path B (blue),
#' and |Path AB| (green) change in value for each PDM
#'
#' @export
#'


plotMedPathCoeff <- function(pathCoeff){

  if ('redFeatWeights' %in% names(pathCoeff)){
    stop('This function only accepts $pathCoeff list.
         Do not input entire output list from getDirectionsofMed()')
  }

  # Number of PDMs
  num <- length(pathCoeff)

  # Unlist pathCoeff
  pathCoeff <- matrix(unlist(pathCoeff), nrow = 5, ncol = num)

  # Subplot
  graphics::par(mfrow=c(1,3))

  # Plot Path A
  plot(pathCoeff[3,],
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
  plot(pathCoeff[4,],
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
  plot(abs(pathCoeff[5,]),
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
