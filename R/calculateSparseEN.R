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
#' @param PDMs Original PDM values; Calculated from non sparse features  
#' @param nPDM number of PDMs
#'
#' @importFrom glmnet cv.glmnet
#'
#' @return A list
#' \itemize{
#'     \item spfeatWeights: sparse weight vectors for each PDM
#'     \item spPathCoeff: resulting mediation pathway coefficients after sparse 
#'     weights calculated
#'     \item spPDM: resulting PDM values after sparse weights are calculated
#'     }
#'
#' @noRd

.calculateSparseEN <- function(x, y, m, PDMs, nPDM){
  
  # Define variables
  B_EN <- NULL
  PDM_EN <- NULL
  pathCoeff_EN <- NULL
  
  # Unlist M if neccessary
  if (is.list(m)){
    M <- matrix(unlist(m), ncol = length(m))
  }
  
  for (i in 1:nPDM){
    
    # Do a cross validation elastic net fit with PDM ~ M
    cv_fit <- cv.glmnet(M,PDMs[[i]])
    # What are the coefficients that correspond to the lambda with the lowest MSE
    B <- coef(cvfit, s = "lambda.min")
    # Remove the Intercept
    B <- B[-1,]
    B_EN[[i]] <- B
    
    # Calculate the new PDM 
    PDM <-  M %*% B
    PDM_EN[[sprintf('PDM%d',i)]] <- PDM
    
    # Calculate the new path coefficient 
    c <- summary(lm(y ~ x))$coefficients[2,1]
    a <- summary(lm(PDM ~ x))$coefficients[2,1]
    mod2_coefficients <- summary(lm(y ~ x + PDM))$coefficients
    c_prime <- mod2_coefficients[2,1]
    b <- mod2_coefficients[3,1]
    
    pathCoeff_EN[[i]] <- c(c, c_prime, a, b, a*b )
    
    
  }
  
  return(list('spFeatWeights' = B_EN,
              'spPathCoeff' = pathCoeff_EN,
              'spPDM' = PDM_EN))
  
}