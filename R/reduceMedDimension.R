#' reduceMedDimension
#'
#' @description Reduces dimensions of radiomic mediator space using Singular Value Decomposition (SVD).
#' Moves treatment, outcome, and radiomic mediator data into a list to be used as an input in later radPDM steps.
#'
#' @param treatment vector (N x 1) Treatment vector for N Subjects
#' @param outcome vector (N x 1) Outcome vector for N Subjects
#' @param mediator matrix (N x p) Mediator matrix for N subjects and p features (Assumed p >> N)
#' @param exVar Percentage of Mediator matrix variance user wants SVD to explain
#'
#' @return A list
#' \itemize{
#'     \item X - vector (N x 1) Treatment vector for N subjects
#'     \item Y - vector (N x 1) Outcome vector for N subjects
#'     \item M - matrix (N x p) Original Mediator matrix for N subjects and p features
#'     \item M_tilde - matrix (N x b) Mediator matrix with reduced dimension space of b features where N > b
#'     \item Dt - matrix (p x b)' Eigenvector space from SVD that maps M_tilde (reduced mediator matrix) back to the original Mediator matrix
#'     \item exVar - Percentage of Mediator matrix variance user wants SVD to explain
#'     }
#'
#' @export
#'


reduceMedDimension <- function(treatment, outcome, mediator, exVar = 0.9){


  ## Check the data type for X, Y, and M

  # X is numeric, M is numeric OR a list, Y is numeric
  if (is.numeric(treatment) && (is.numeric(mediator) | is.list(mediator)) && is.numeric(outcome)){
  } else {
    stop('Check X, Y, and M input structure. X and Y should be numeric. M should be numeric or list ')
  }

  # Check X is a matrix
  if (!(is.matrix(treatment))){
    X <- as.matrix(treatment)
  } else{
    X <- treatment
  }

  # Check Y is a matrix
  if (!(is.matrix(outcome))){
    Y <- as.matrix(outcome)
  } else{
    Y <- outcome
  }

  ## Reduce dimensions using SVD

  pracma::fprintf('SVD for %d subjects...\n', length(treatment))

  # Calculate the SVD of Mediator matrix
  svd_results <- svd(t(mediator))
  U <- svd_results$u
  diag_S <- svd_results$d
  V <- svd_results$v

  # Make sure the user input for exvar parameter is between 0 and 1
  if (abs(exVar) > 1){
    stop('exVar input (percent variance explained) needs to be between 0 and 1')
  } else {
    # Calculate the number of components that explains 90% variance
    nComps <- (which(cumsum(diag_S)/sum(diag_S) >= exVar))[1]
  }

  M_tilde <- V[,1:nComps] %*% diag(diag_S)[1:nComps,1:nComps]
  Dt <- t(U[,1:nComps])

  # Constrain the Dimensions to at most 50 (What is smaller -> 50 components or 90%)
  svddim <- min(50, dim(M_tilde)[2])
  M_tilde <- M_tilde[,1:svddim]
  Dt <- Dt[1:svddim,]

  # Print statement to convey how many features are being kept
  pracma::fprintf('Keeping %d components.\n',svddim)

  # Return a list
  return(list('X'= X, "Y" = Y, "M" = mediator, "M_tilde" = M_tilde, "Dt" = Dt, 'exVar' = exVar))

}
