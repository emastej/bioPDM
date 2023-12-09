#' .PDMN
#'
#' @description Objective function for radPDM. Optimization algorithm that finds
#' PDM weights that maximize |ab|
#' while being orthogonal to previous PDM weights.
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param m matrix (N x b) Mediator matrix with reduced dimension space of b
#' features
#' @param W Weights of previously calculated directions of mediation
#' @param numCores Number of cores to use for parallel execution. Default is
#' a single core.
#' @param timeout NEED TO ADD DOCUMENTATION
#'
#' @importFrom R.utils withTimeout
#' @importFrom stats rnorm
#' @importFrom pracma fmincon
#' @importFrom pracma pinv
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @return A list
#' \itemize{
#'     \item weights -  Weights for the Kth direction of mediation
#'     \item path_coeff - Pathway coefficients for Nth direction of mediation
#'     (c, c', a, b)
#'     \item init_weights - Initial weights used for PDM calculation optimization
#'     }
#'
#' @noRd
#'


.PDMN <- function(x, y, m, W, numCores, timeout) {

  ## Set Up the initial variables

  num_M_feat <- dim(m)[2] # Number of mediator features (SVD components)
  num_calc_pdms <- length(W)  # Number of PDMs previously calculated
  Z <- 50          # Number of times to repeat optimization
  vec_of_ones <- rep(1, dim(m)[1])  # Design matrix for first regression
  X1 <- cbind(vec_of_ones, x)


  ## Create orthogonality constraints


  # A and b will be empty on the first PDM
  # These are for orthogonality constraints
  A <- NULL
  b <- NULL

  # num_calc_pdms is the number of PDMs already calculated
  # Not relevant for first PDM calculation
  if (num_calc_pdms > 0){

    # Unlist matrix of the previously calculated PDM weights
    A <- t(matrix(unlist(W), ncol = num_calc_pdms))

    # Create a vector of zeroes for the orthogonality constraint
    b <- rep(0, num_calc_pdms)
  }

  ## Find optimal PDM weight values

  # Set max_indirect_effect to 0. If |ab| > 0, then max_indirect_effect will be replaced
  max_indirect_effect <- 0

  # Initialize clusters
  cluster <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cluster)

  # Prepare cluster environment
  parallel::clusterExport(cluster, c("num_M_feat", "m", "y", "X1", "A", "b",
                                     "max_indirect_effect", "objfun", "ceq",
                                     "fmincon", "pinv"),
                          envir = environment())

  # Parallelized loop
  results <- parallel::parLapply(cluster, 1:Z, fun = function (k, timeout) {
  # results <- lapply(1:Z, FUN = function (k, timeout) {

    # Generate random values until fmincon() returns non-imaginary
    # results
    start_user_time <- proc.time()[["elapsed"]]
    success <- FALSE

    while (!success &
           ((proc.time()[["elapsed"]] - start_user_time) <= 100)) {

      # Generate random initial values where mu = 0 and std = 1
      rand_init_vars <- stats::rnorm(num_M_feat, 0, 1)
      rand_init_vars <- rand_init_vars / sqrt(sum(rand_init_vars^2))

      WM <- rand_init_vars

      optim_result <- tryCatch({

          # Limit time spent on this - if it fails just try again
          result <- R.utils::withTimeout({
            pracma::fmincon(WM, objfun, m = m, y = y, X1 = X1,
                            Aeq = A, beq = b, heq = ceq)
            }, timeout = timeout)

          list("success" = TRUE, result = result)

        }, error = function(e) {

          list("success" = FALSE, result = NA)

        })


      if (optim_result[["success"]]) {
        success <- TRUE
      }

    }

    # If no non-complex results are found within timeout length, return NULL
    if (!success) {
      warning("Some random number generation for optimization failed, results are fewer than expected")
      return(NULL)
    } else {
      optim_result <- optim_result[["result"]]
    }

    # Weights
    feat_weights <- optim_result$par

    # Solve for Total Effect aka Gamma (y = gamma * x)
    total_effect <- pracma::pinv(X1) %*% y


    # calc_pdm is the calculated PDM from the calculated weights
    calc_pdm <- m %*% feat_weights

    # Solve for Alpha: [n x b] x [b x 1] (M = alpha * X)
    alpha_indirect_effect <- pracma::pinv(X1) %*% calc_pdm

    # Solve for Beta and Direct Effect aka Gamma' (y = beta * M + gamma' * X)
    X2 <- cbind(X1, calc_pdm)
    eq_coeff <- pracma::pinv(X2) %*% y
    beta_indirect_effect <- eq_coeff[3]
    direct_effect <- eq_coeff[2]

    # Quick Proof
    # y = beta[1] + beta[2]*x + beta[3]*m
    # where beta[1] = beta_0, beta[2] = c', beta[3] = b

    # Solve for Indirect Effect
    indirect_effect <- abs(alpha_indirect_effect[2]*beta_indirect_effect)

    # If indirect_effect > max_indirect_effect, Keep that iteration info
    if (indirect_effect > max_indirect_effect){
      max_indirect_effect <- indirect_effect
      # path_coeff = gamma, gamma' , alpha, beta
      path_coeff <- c(total_effect[2], direct_effect,
                 alpha_indirect_effect[2], beta_indirect_effect)
      w_k <- feat_weights
      WM_init <- rand_init_vars
    }

    # Return final PDM weights, theta values,  and initial values
    return(list('weights' = w_k,
                'path_coeff' = path_coeff,
                'init_weights' = WM_init,
                'max_indirect_effect' = max_indirect_effect))

  }, timeout = timeout)

  # Stop cluster
  parallel::stopCluster(cluster)

  # Choose final weights with largest indirect effect
  not_nullIndex <- unlist(lapply(results, function (x) { !is.null(x) }))
  results <- results[not_nullIndex]
  max_indirect_effect <- unlist(lapply(results, function (x) { x$max_indirect_effect }))
  max_indirect_effectIndex <- which(max_indirect_effect == max(max_indirect_effect))
  final_results <- results[[max_indirect_effectIndex]]

  # Return final PDM weights, theta values,  and initial values
  final_results$max_indirect_effect <- NULL
  return(final_results)

}


# Functions for the optimization algorithm
################################################################################

# Objective Function for Optimization of max of |ab|
objfun <- function(WM,m,y,X1){

  calc_pdm <- m %*% WM
  alpha <- pracma::pinv(X1) %*% calc_pdm
  r <- calc_pdm - X1 %*% alpha
  beta <- (sum(r*y)) /(sum(r*r))
  mval <- -abs(alpha[2]*beta)

  return(mval)
}

# Constraints of Objective function (Weights must be normalized to 1)
ceq <- function(WM,m,y,X1){
  return( t(WM) %*% WM - 1 )

}

