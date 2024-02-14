#' .calculateSparseEN
#'
#' @description Internal function to find sparse feature weights and calculate 
#' the resulting PDMs and mediation pathway coefficients using the elastic net 
#' method
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param m matrix (N x b) Mediator matrix with reduced dimension space of b
#' features
#' @param W Nonsparse feature weights 
#'
#' @importFrom pracma pinv
#'
#' @return A list
#' \itemize{
#'     \item featWeights: sparse weight vectors for each PDM
#'     \item pathCoeff: resulting mediation pathway coefficients after sparse 
#'     weights calculated
#'     \item PDM: resulting PDM values after sparse weights are calculated
#'     }
#'
#' @noRd

.calculateSparseEN <- function(x, y, m, W){
  
}