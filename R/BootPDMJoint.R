#' .BootPDMJoint
#'
#' @description Perform bootstrap process to calculate stats (95% CI, mean, std)
#'  on JointPDM weights
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param M_reduced matrix (N x b) Mediator matrix with reduced dimension space of
#' b features
#' @param redFeatWeights matrix (b x k) Weights for each b feature calculated for each of the
#' k PDMs
#' @param tLoadingMatrix matrix (p x b)' Eigenvector space from SVD that maps M_reduced back
#' to the original Mediator matrix
#' @param bootSamp Number of bootstrap samples
#' @param initValues matrix (b x k) Initial weights used for PDM calculation optimization
#'
#' @importFrom stats t.test
#' @importFrom stats sd
#'
#' @return A list
#' \itemize{
#'     \item weightStats - 95% CI, mean, and std of bootstrap results
#'     \item weightBootSamples - Resulting weights from bootstrap process
#'     \item pathBootSamples - Resulting mediation path coefficients from
#'     bootstrap process
#'     }
#'
#' @noRd
#'

.BootPDMJoint <- function(x,
                          y,
                          M_reduced,
                          redFeatWeights,
                          tLoadingMatrix,
                          bootSamp,
                          initValues){

  # Number of features
  num_feat <- dim(tLoadingMatrix)[2]
  # Number of subjects
  num_sub <- length(y)
  # Number of Caluclated PDMs
  nPDM <- length(redFeatWeights)
  # JointPDM weights calculated after bootstrapping
  jointPDM_weight_samples <- matrix(0,num_feat,bootSamp)
  # jointPDM_path_coeff_samples dimensions are 4 x bootstrap sample size x nPDM
  jointPDM_path_coeff_samples <- array(0, dim = c(4,bootSamp, nPDM))
  jointPDM_weight_stats <- vector("list", num_feat)

  # Create matrix of randomly generated samples
  rand_samp <- matrix(0,num_sub,bootSamp)

  loading_matrix <- t(tLoadingMatrix)
  cat('\n')

  # Bootstrap samples of the data by randomly sample subjects with replacement
  for (i in 1:bootSamp){
    ransamp <- sample(1:num_sub, size = num_sub, replace = TRUE)
    rand_samp[,i] = ransamp
  }

  # Get the Bootstrapped Samples
  for (i in 1:bootSamp){

    # Print current bootstrapping iteration
    cat('\r', 'Bootstrap Sample:', i, '/', bootSamp)

    tryCatch(expr = {

      # Randomly sampled subjects
      subject_index <- rand_samp[,i]

      # Index data using bootstrap samples
      x_indexed <- x[subject_index]
      y_indexed <- y[subject_index]
      M_indexed <- M_reduced[subject_index,]

      # Create matrices for PDM weights and Thetas
      weights_from_boot_sample <- matrix(0, dim(M_reduced)[2], nPDM)
      path_coeff_from_boot_sample <- matrix(0, 4, nPDM)

      ## Loop through the PDMs for each sample
      ## Bootstrap and calculate EACH PDM from resampled Data

      for (k in 1:nPDM){
        if (k==1){
          boot_jpdmn_result <-  .PDMNboot(x = x_indexed,
                                          y = y_indexed,
                                          m = M_indexed,
                                          W = NULL,
                                          initialValues = initValues[k])
          weights_from_boot_sample[,k] <- boot_jpdmn_result[['weights']]
          path_coeff_from_boot_sample[,k] <- boot_jpdmn_result[['path_coeff']]
        } else {
          boot_jpdmn_result <-  .PDMNboot(x  = x_indexed,
                                          y = y_indexed,
                                          m = M_indexed,
                                          W = redFeatWeights[1:k-1],
                                          initialValues = initValues[k])
          weights_from_boot_sample[,k] <- boot_jpdmn_result[['weights']]
          path_coeff_from_boot_sample[,k] <- boot_jpdmn_result[['path_coeff']]
        }
      }

      # get path coefficients
      a_n <- path_coeff_from_boot_sample[3,]
      b_n <- path_coeff_from_boot_sample[4,]

      # Calculate the Full Feature Weight Vector
      # [p x b] x ( [b x npdm] x [npdm x 1] ) = [p x 1]
      jointPDM_weight_samples[,i] <- (loading_matrix %*% (weights_from_boot_sample %*% (a_n * b_n)))
      jointPDM_path_coeff_samples[,i,] <- path_coeff_from_boot_sample

    },
    error=function(e){})
  }

  # Remove any columns of zeroes that may have happened due to iteration skipping
  jointPDM_weight_samples = jointPDM_weight_samples[,-(which(colSums(jointPDM_weight_samples)==0))]
  jointPDM_path_coeff_samples = jointPDM_path_coeff_samples[,-(which(colSums(jointPDM_weight_samples)==0)),]


  ## Calculate the 95% CI, mean and std for the bootstrapped joint feature weights


  for (z in 1:num_feat){
    ci <- stats::t.test(jointPDM_weight_samples[z,])
    jointPDM_weight_stats[[z]][['95CI_lb']] <- ci$conf.int[1]
    jointPDM_weight_stats[[z]][['95CI_ub']] <- ci$conf.int[2]
    jointPDM_weight_stats[[z]][['mean']] <- mean(jointPDM_weight_samples[z,])
    jointPDM_weight_stats[[z]][['sd']] <- stats::sd(jointPDM_weight_samples[z,])
    names(jointPDM_weight_stats)[z] <- sprintf('feat_%d',z)
  }

  return(list('weightStats' = jointPDM_weight_stats,
              'weightBootSamples' = jointPDM_weight_samples,
              'pathBootSamples'= jointPDM_path_coeff_samples))

}
