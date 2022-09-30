#' .runPDM
#'
#' @description  Calculates PDM weights for each number of desired PDMs and
#' organizes results
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param M_tilde matrix (N x b) Mediator matrix with reduced dimension space of
#'  b features
#' @param Dt matrix (p x b)' Eigenvector space from SVD that maps M_tilde back
#' to the original Mediator matrix
#' @param nPDM Number of PDMs to calculate
#' @param doJointPDM Calculate the JointPDM. Can be True or False
#'
#' @return A list
#' \itemize{
#'     \item W - matrix (b x q) Weights for each b feature calculated for each
#'     of the q PDMs
#'     \item Theta - matrix (5 x q) Mediation Path Coefficients for each of the
#'     q PDMs (gamma, gamma', alpha, beta, ab)
#'     \item Wfull - matrix (p x q) Weights mapped back onto original
#'     Mediator Matrix
#'     \item WMinit - matrix (b x q) Initial weights used for PDM calculation
#'     optimization
#'     \item WfullJoint - vector (p x 1) Weights for the JointPDM. Each Mediator
#'     feature has one weight (if applicable)
#'     }
#'
#' @noRd
#'


.runPDM <- function(x, y, M_tilde, Dt, nPDM, doJointPDM){


  ## Calculate the Principal Directions of Mediation


  cat('Computing ',nPDM, ' PDMs...')

  # Initialize variables

  # W is a matrix of weights from previous PDM calculated. Will be empty for
  # first PDM
  W <- NULL
  # gamma, alpha, beta values
  Theta <- NULL
  # Full weights mapped back onto original Mediators
  # (not just the mediator matrix with reduced dimension space)
  Wfull <- NULL
  # Initial Weights randomly generated in the optimization process
  WMinit <- NULL

  # Transpose Dt
  D <-  t(Dt)

  # Loop through number of PDMs that are to be calculated
  for (k in 1:nPDM){

    # Print the number of PDM currently being calculated
    cat(' ',k)

    # Use PDMN function to calculate each PDM. Output is a list
    pdmk_results <- .PDMN(x, y, M_tilde, W)

    # Weights for kth direction
    w_k <- pdmk_results[['weights']]
    # Theta values for kth PDM. theta = gamma, gamma' , alpha, beta
    theta_k <- pdmk_results[['theta']]
    # Initial weight values input into optim alg
    wmi <- pdmk_results[['init_weights']]

    # Check that alpha is positive
    if (sign(theta_k[3]) == -1){
      theta_k[3:4] <- -1 * theta_k[3:4]
      w_k <- -1 * w_k
    }

    # Add Kth PDM Weights to W matrix
    W[[k]] <-  w_k

    # Add a*b to the theta value (gamma, gamma', alpha, beta, ab)
    Theta[[k]] <- c(theta_k, theta_k[3]*theta_k[4])

    # Weights mapped back onto original mediators
    Wfull[[k]] <- D %*% w_k

    # Initial weights
    WMinit[[k]] <- wmi
  }


  ## Compute JointPDM (if applicable)
  if (doJointPDM){

    cat('.\nComputing jointPDM')

    WAll <- matrix(unlist(W), nrow = dim(Dt)[1], ncol = nPDM )
    tf <- matrix(unlist(Theta), nrow = 5, ncol = nPDM )
    a <- tf[3,]
    b <- tf[4,]
    WJoint <- D %*% (t(((a*b) %*% t(WAll))))

    # Quick Linear Algebra Proof
    # [p x b] x ([1 x q] x [q x b])'  = [p x b] x [b x 1] = p x 1

  } else {
    WJoint <- NULL
  }

  # Return a list with final weights, theta,
  # Weights mapped back onto the original mediators, initial weights,
  # and JointPDM weights
  return(list('w_k' = W, 'w_k*Dt' = Wfull, 'JointW' = WJoint,
              'Theta' = Theta, 'WMinit' = WMinit ))

  pracma::fprintf(' - Done.')
}
