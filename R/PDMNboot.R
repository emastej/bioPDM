#' .PDMNboot
#'
#' @description Objective function for radPDM. Optimization algorithm that finds PDM weights that maximize |ab|
#' while being orthogonal to previous PDM weights. Uses initial weights from PDMN
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param m matrix (N x b) Mediator matrix with reduced dimension space of b features
#' @param W Weights of previously calculated directions of mediation
#' @param initialValues matrix (b x q) Initial weights used for PDM calculation optimization
#'
#' @return A list
#' \itemize{
#'     \item weights: 95% CI, mean, and std of bootstrap results for each feature weight for each PDM
#'     \item theta: pathway effects from each bootstrapped sample
#'     }
#'
#' @noRd
#'

.PDMNboot <- function(x, y, m, W, initialValues){

  if (typeof(initialValues)=='list'){
    initialValues = unlist(initialValues)
  }

  ## Set Up the initial variables

  set.seed(12345)

  len <- dim(m)[2] # Number of mediators (or SVD components)
  N1 <- length(W)  # Number of PDMs previously calculated
  J <- rep(1, dim(m)[1])  # Design matrix for first regression
  X1 <- cbind(J, x)

  ## Create orthogonality constraints

  # A and b will be empty on the first PDM
  A <- NULL
  b <- NULL

  # N1 is the number of PDMs already calculated
  # Not relevant for first PDM calculation
  if (N1 > 0){

    #Create a matrix of previoulsy calculated weights
    A <- matrix(0, len, N1)
    for (i in 1:N1){
      A[,i] <- W[[i]]
    }
    A  <- t(A)

    # Create a vector of zeroes for the orthogonality constraint
    b <- rep(0,N1)
  }


  ## Find weights that maximize |ab| given the initial values of the optimization algorithm are the same from PDMN.R


  optim_result <- pracma::fmincon(initialValues, objfun, m = m, y = y, X1 = X1, Aeq = A, beq = b, heq = ceq)

  # Weights
  weights <- optim_result$par

  # Find Gamma
  gamma <- pracma::pinv(X1) %*% y

  # Find Alpha: [n x b] x [b x 1]
  mw_N <- m %*% weights
  # Solve for alpha since M = X * alpha
  alpha <- pracma::pinv(X1) %*% mw_N

  # Find Beta and Gamma'
  X2 <- cbind(X1, mw_N)

  beta <- pracma::pinv(X2) %*% y

  ab <- abs(alpha[2]*beta[3])

  # theta = gamma, gamma' , alpha, beta
  theta <- c(gamma[2], beta[2], alpha[2], beta[3])
  w_N <- weights

  # Return the final PDM weights and theta values
  return(list('weights' = w_N, 'theta' = theta))

}

# Functions for the optimization algorithm
###############################################################################

# Objective Function for Optimization of max of |ab|
objfun <- function(initialValues,m,y,X1){

  mw_N <- m %*% initialValues
  alpha <- pracma::pinv(X1) %*% mw_N
  r <- mw_N - X1 %*% alpha
  beta <- (sum(r*y)) /(sum(r*r))
  mval <- -abs(alpha[2]*beta)

  return(mval)
}

# Constraints of Object function (Weights must be normalized to 1)
ceq <- function(initialValues,m,y,X1){
  return( t(initialValues) %*% initialValues - 1 )

}

