#' .calculateSparseThresh
#'
#' @description Internal function to find sparse feature weights and calculate 
#' the resulting PDMs and mediation pathway coefficients using the weight threshold
#' method
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param m matrix (N x b) Mediator matrix with reduced dimension space of b
#' features
#' @param W Nonsparse feature weights 
#' @param nPDM number of PDMs
#'
#'
#' @return A list
#' \itemize{
#'     \item spFeatWeights: sparse weight vectors for each PDM
#'     \item spPathCoeff: resulting mediation pathway coefficients after sparse 
#'     weights calculated
#'     \item spPDM: resulting PDM values after sparse weights are calculated
#'     }
#'
#' @noRd

.calculateSparseThresh <- function(x, y, m, W, nPDM){
  
  # Define variables
  B_thresh <- W
  PDM_thresh <- NULL
  pathCoeff_thresh <- NULL
  
  # Define the weight threshold value 
  threshold <- 1/sqrt(length(m))
  
  # Unlist M if neccessary
  if (is.list(m)){
    M <- matrix(unlist(m), ncol = length(m))
  }
  
  for (i in 1:nPDM){
    
    # Set any weights whos absolute value is less than threshold to zero
    B_thresh_i <- B_thresh[[i]]
    B_thresh_i[base::abs(B_thresh_i) < threshold] <- 0
    
    # Calculate the new PDM 
    PDM <-  M %*% B_thresh_i
    PDM_thresh[[sprintf('PDM%d',i)]] <- PDM
    
    # Calculate the new path coefficient 
    c <- summary(lm(y ~ x))$coefficients[2,1]
    a <- summary(lm(PDM ~ x))$coefficients[2,1]
    mod2_coefficients <- summary(lm(y ~ x + PDM))$coefficients
    c_prime <- mod2_coefficients[2,1]
    b <- mod2_coefficients[3,1]
    
    pathCoeff_thresh[[i]] <- c(c, c_prime, a, b, a*b )
    
    B_thresh[[i]] <- B_thresh_i
    
  }
  
  return(list('spFeatWeights' = B_thresh,
              'spPathCoeff' = pathCoeff_thresh,
              'spPDM' = PDM_thresh))
  
}



