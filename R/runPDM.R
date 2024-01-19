#' .runPDM
#'
#' @description  Calculates PDM weights for each number of desired PDMs and
#' organizes results
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param M matrix (N x b) Mediator matrix with reduced dimension space of
#'  b features
#' @param tLoadingMatrix matrix (p x b)' Eigenvector space from SVD that maps M back
#' to the original Mediator matrix
#' @param nPDM Number of PDMs to calculate
#' @param doJointPDM Calculate the JointPDM. Can be True or False
#' @param numCores Number of cores to use for parallel execution. Default is
#' a single core
#'
#' @importFrom pracma fprintf
#' @importFrom MASS ginv
#'
#' @return A list
#' \itemize{
#'     \item redFeatWeights - matrix (b x k) Weights for each b feature calculated for each
#'     of the k PDMs
#'     \item pathCoeff - matrix (5 x k) Mediation Path Coefficients for each of the
#'     k PDMs (gamma, gamma', alpha, beta, ab)
#'     \item featWeights - matrix (p x k) Weights mapped back onto original
#'     Mediator Matrix
#'     \item initValues - matrix (b x k) Initial weights used for PDM calculation
#'     optimization
#'     \item jointWeights - vector (p x 1) Weights for the JointPDM. Each Mediator
#'     feature has one weight (if applicable)
#'     }
#'
#' @noRd
#'


.runPDM <- function(x, y, M, tLoadingMatrix, nPDM, doJointPDM, numCores,
                    timeout){


  ## Calculate the Principal Directions of Mediation


  cat('Computing ',nPDM, ' PDMs...')

  # Initialize variables

  # Feat_weights is a matrix of weights from previous PDM calculated.
  # Will be empty for first PDM
  Feat_weights <- NULL
  # Path Coeff: gamma, alpha, beta values (Total, Indirect, and Direct Effects)
  Path_coeff <- NULL
  # Full weights mapped back onto original Mediators
  # (not just the mediator matrix with reduced dimension space)
  Full_feat_weights <- NULL
  # Initial Weights randomly generated in the optimization process
  Init_weight_values <- NULL

  # Inverse tLoadingMatrix
  loading_matrix <-  MASS::ginv(tLoadingMatrix)

  # Loop through number of PDMs that are to be calculated
  for (k in 1:nPDM){

    # Print the number of PDM currently being calculated
    cat(' ', k)

    # Use PDMN function to calculate each PDM. Output is a list
    pdmk_results <- .PDMN(x, y, m = M, W = Feat_weights, numCores, timeout)

    # Weights for kth direction
    weights_k <- pdmk_results[['weights']]
    # Path_coeff values for kth PDM. Path_coeff = gamma, gamma' , alpha, beta
    path_coeff_k <- pdmk_results[['path_coeff']]
    # Initial weight values input into optim alg
    init_weights_k <- pdmk_results[['init_weights']]

    # Check that alpha is positive
    if (sign(path_coeff_k[3]) == -1){
      path_coeff_k[3:4] <- -1 * path_coeff_k[3:4]
      weights_k <- -1 * weights_k
    }

    # Add Kth PDM Weights to W matrix
    Feat_weights[[k]] <-  weights_k

    # Add a*b to the path_coeff value (gamma, gamma', alpha, beta, ab)
    Path_coeff[[k]] <- c(path_coeff_k, path_coeff_k[3]*path_coeff_k[4])

    # Weights mapped back onto original mediators
    Full_feat_weights[[k]] <- loading_matrix %*% weights_k

    # Initial weights
    Init_weight_values[[k]] <- init_weights_k
  }


  ## Compute JointPDM (if applicable)
  if (doJointPDM){

    cat('.\nComputing jointPDM')

    all_pdm_weights <- matrix(unlist(Feat_weights), nrow = dim(tLoadingMatrix)[1], ncol = nPDM )
    path_coeff_matrix <- matrix(unlist(Path_coeff), nrow = 5, ncol = nPDM )
    alpha <- path_coeff_matrix[3,]
    beta <- path_coeff_matrix[4,]
    joint_weight_values <- loading_matrix %*% (t(((alpha*beta) %*% t(all_pdm_weights))))

    # Quick Linear Algebra Proof
    # [p x b] x ([1 x k] x [k x b])'  = [p x b] x [b x 1] = p x 1

  } else {
    joint_weight_values <- NULL
  }

  # Return a list with final weights, path_coeff,
  # Weights mapped back onto the original mediators, initial weights,
  # and JointPDM weights
  return(list('redFeatWeights' = Feat_weights,
              'featWeights' = Full_feat_weights,
              'jointWeights' = joint_weight_values,
              'pathCoeff' = Path_coeff,
              'initValues' = Init_weight_values ))

  pracma::fprintf(' - Done.')
}
