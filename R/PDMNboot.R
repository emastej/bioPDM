#' .PDMNboot
#'
#' @description Objective function for radPDM. Optimization algorithm that finds
#' PDM weights that maximize |ab|
#' while being orthogonal to previous PDM weights. Uses initial weights from PDMN
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param m matrix (N x b) Mediator matrix with reduced dimension space of b
#' features
#' @param W Weights of previously calculated directions of mediation
#' @param initialValues matrix (b x k) Initial weights used for PDM calculation
#' optimization
#'
#' @importFrom pracma fmincon
#' @importFrom pracma pinv
#'
#' @return A list
#' \itemize{
#'     \item weights: 95% CI, mean, and std of bootstrap results for each
#'     feature weight for each PDM
#'     \item path_coeff: pathway effects from each bootstrapped sample
#'     }
#'
#' @noRd
#'

.PDMNboot <- function(x, y, m, W, initialValues, timeout) {

  if (typeof(initialValues)=='list'){
    initialValues = unlist(initialValues)
  }

  ## Set Up the initial variables

  num_calc_pdms <- length(W)  # Number of PDMs previously calculated
  vector_ones <- rep(1, dim(m)[1])  # Design matrix for first regression
  X1 <- cbind(vector_ones, x)

  ## Create orthogonality constraints

  # A and b will be empty on the first PDM
  A <- NULL
  b <- NULL

  # num_calc_pdms is the number of PDMs already calculated
  # Not relevant for first PDM calculation
  if (num_calc_pdms > 0){

    # Unlist weight list
    A <- t(matrix(unlist(W), ncol = num_calc_pdms))

    # Create a vector of zeroes for the orthogonality constraint
    b <- rep(0,num_calc_pdms)
  }


  ## Find weights that maximize |ab| given the initial values of the
  ## optimization algorithm are the same from PDMN.R

  optim_result <- tryCatch({

      # Limit time spent on this - if it fails just try again
      result <- R.utils::withTimeout({
        pracma::fmincon(initialValues,
                        objfun,
                        m = m,
                        y = y,
                        X1 = X1,
                        Aeq = A,
                        beq = b,
                        heq = ceq)
        }, timeout = timeout)

      # Weights
      weights <- result$par

      # Find total_effect
      total_effect <- pracma::pinv(X1) %*% y

      # Find alpha_indirect_effect: [n x b] x [b x 1]
      calc_pdm <- m %*% weights

      # Solve for alpha since M = X * alpha_indirect_effect
      alpha_indirect_effect <- pracma::pinv(X1) %*% calc_pdm

      # Find Beta and total_effect
      X2 <- cbind(X1, calc_pdm)

      eq_coeff <- pracma::pinv(X2) %*% y
      beta_indirect_effect <- eq_coeff[3]
      direct_effect <- eq_coeff[2]

      ab <- abs(alpha_indirect_effect[2]*beta_indirect_effect)

      # theta = total_effect, direct_effect , alpha_indirect_effect, beta
      path_coeff <- c(total_effect[2], direct_effect,
                      alpha_indirect_effect[2], beta_indirect_effect)

      # Return relevant values
      list("success" = TRUE, 'weights' = weights, 'path_coeff' = path_coeff)

    }, error = function(e) {

      # If the above errors out (takes too much time) set success to FALSE
      # and return nothing
      list("success" = FALSE, 'weights' = NA, 'path_coeff' = NA)

    })


  # optim_result <- pracma::fmincon(initialValues,
  #                                 objfun,
  #                                 m = m,
  #                                 y = y,
  #                                 X1 = X1,
  #                                 Aeq = A,
  #                                 beq = b,
  #                                 heq = ceq)
  #
  # # Weights
  # weights <- optim_result$par

  # # Find total_effect
  # total_effect <- pracma::pinv(X1) %*% y

  # # Find alpha_indirect_effect: [n x b] x [b x 1]
  # calc_pdm <- m %*% weights
  # # Solve for alpha since M = X * alpha_indirect_effect
  # alpha_indirect_effect <- pracma::pinv(X1) %*% calc_pdm

  # # Find Beta and total_effect'
  # X2 <- cbind(X1, calc_pdm)
  #
  # eq_coeff <- pracma::pinv(X2) %*% y
  # beta_indirect_effect <- eq_coeff[3]
  # direct_effect <- eq_coeff[2]
  #
  # ab <- abs(alpha_indirect_effect[2]*beta_indirect_effect)
  #
  # # theta = total_effect, direct_effect , alpha_indirect_effect, beta
  # path_coeff <- c(total_effect[2], direct_effect, alpha_indirect_effect[2], beta_indirect_effect)

  # Return the final PDM weights and theta values
  return(optim_result)

}

# Functions for the optimization algorithm
###############################################################################

# Objective Function for Optimization of max of |ab|
objfun <- function(initialValues,m,y,X1){

  calc_pdm <- m %*% initialValues
  alpha <- pracma::pinv(X1) %*% calc_pdm
  r <- calc_pdm - X1 %*% alpha
  beta <- (sum(r*y)) /(sum(r*r))
  mval <- -abs(alpha[2]*beta)

  return(mval)
}

# Constraints of Object function (Weights must be normalized to 1)
ceq <- function(initialValues,m,y,X1){
  return( t(initialValues) %*% initialValues - 1 )

}

