#' .PDMN
#'
#' @description Objective function for radPDM. Optimization algorithm that finds
#'  PDM weights that maximize |ab|
#' while being orthogonal to previous PDM weights.
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param m matrix (N x b) Mediator matrix with reduced dimension space of b
#' features
#' @param W Weights of previously calculated directions of mediation
#'
#' @importFrom R.utils withTimeout
#' @importFrom stats rnorm
#' @importFrom pracma fmincon
#' @importFrom pracma pinv
#'
#' @return A list
#' \itemize{
#'     \item weights -  Weights for the Nth direction of mediation
#'     \item theta - Pathway coefficients for Nth direction of mediation
#'     (c, c', a, b)
#'     \item init_weights - Initial weights used for PDM calculation optimization
#'     }
#'
#' @noRd
#'


.PDMN <- function(x, y, m, W){

  ## Set Up the initial variables

  len <- dim(m)[2] # Number of mediators (SVD components)
  N1 <- length(W)  # Number of PDMs previously calculated
  K <- 50          # Number of times to repeat optimization
  J <- rep(1, dim(m)[1])  # Design matrix for first regression
  X1 <- cbind(J, x)


  ## Create orthogonality constraints


  # A and b will be empty on the first PDM
  A <- NULL
  b <- NULL

  # N1 is the number of PDMs already calculated
  # Not relevant for first PDM calculation
  if (N1 > 0){

    # Create a matrix of the previously calculated PDM weights
    A <- matrix(0, len, N1)
    for (i in 1:N1){
      A[,i] <- W[[i]]
    }
    A  <- t(A)

    # Create a vector of zeroes for the orthogonality constraint
    b <- rep(0,N1)
  }


  ## Find optimal PDM weight values


  # Set crtmax to 0. If |ab| > 0, then crtmax will be replaced
  crtmax <- 0

  # Find the the PDM weights that max |ab|
  for (i in 1:K){

    # Error Handling (bypass iterations where the WM cause a complex eigenvalue error)
    tryCatch({

      # Error Handling (Bypass iterations that take way too long)
      # I assume this is from being caught in a local minimum
      tryCatch(
        expr = {
          R.utils::withTimeout({

            # Generate random initial values where mu = 0 and std = 1
            WM <- stats::rnorm(len, 0, 1)
            WM <- WM/sqrt(sum(WM^2))

            WMtmp <- WM

            # Find optimum weights to maximize ab
            optim_result <- pracma::fmincon(WM, objfun,
                                            m = m,
                                            y = y,
                                            X1 = X1,
                                            Aeq = A,
                                            beq = b,
                                            heq = ceq)

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

            # Quick Proof
            # y = beta[1] + beta[2]*x + beta[3]*m
            # where beta[1] = beta_0, beta[2] = c', beta[3] = b

            ab <- abs(alpha[2]*beta[3])

            # Repeat k times. If ab > than previous ab, Keep that iteration info
            if (ab > crtmax){
              crtmax <- ab
              # theta = gamma, gamma' , alpha, beta
              theta <- c(gamma[2], beta[2], alpha[2], beta[3])
              w_N <- weights
              WMi <- WMtmp
            }

          },
          timeout = 2.3)},
        TimeoutException = function(ex) {})}, error=function(e){})

  }

  # Return final PDM weights, theta values,  and initial values
  return(list('weights' = w_N, 'theta' = theta, 'init_weights' = WMi))

}


# Functions for the optimization algorithm
################################################################################

# Objective Function for Optimization of max of |ab|
objfun <- function(WM,m,y,X1){

  mw_N <- m %*% WM
  alpha <- pracma::pinv(X1) %*% mw_N
  r <- mw_N - X1 %*% alpha
  beta <- (sum(r*y)) /(sum(r*r))
  mval <- -abs(alpha[2]*beta)

  return(mval)
}

# Constraints of Objective function (Weights must be normalized to 1)
ceq <- function(WM,m,y,X1){
  return( t(WM) %*% WM - 1 )

}

