#' .BootPDM
#'
#' @description Perform bootstrap process to calculate stats (95% CI, mean, std)
#' on PDM weights
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
#'     \item pathBootSamples - Resulting mediation path coefficients
#'     from bootstrap process
#'     }
#'
#' @noRd
#'

.BootPDM <- function(x,
                     y,
                     M_reduced,
                     redFeatWeights,
                     tLoadingMatrix,
                     bootSamp,
                     initValues){

  # Number of mediation features
  num_feat <- dim(tLoadingMatrix)[2]
  # Number of subjects
  num_sub <- length(y)
  # Matrix for the Weight bootstrapped samples for the kth pdm
  weight_samples_k <- matrix(0,num_feat,bootSamp)
  # Matrix for the path coefficients calculated from the boot samples for kth pdm
  path_samples_k <- matrix(0,4,bootSamp)
  # List for resulting stats (CI, mean, std) from bootstrapped samples for kth pdm
  weight_stats_k <- vector("list", num_feat)
  # Empty matrix of randomly generated samples
  rand_samp <- matrix(0,num_sub,bootSamp)

  loading_matrix <- t(tLoadingMatrix)
  cat('\n')

  # Bootstrap samples of the data by randomly sampling subjects with replacement
  for (i in 1:bootSamp){
    ransamp <- sample(1:num_sub, size = num_sub, replace = TRUE)
    rand_samp[,i] <- ransamp
  }

  # Get the bootstrapped samples
  for (i in 1:bootSamp){

    # Print current bootstrapping iteration
    cat('\r', 'Bootstrap Sample:', i, '/', bootSamp)

    tryCatch(expr={

      # Randomly sampled subjects
      subject_index <- rand_samp[,i]

      # Index data using bootstrap samples
      x_indexed <- x[subject_index]
      y_indexed <- y[subject_index]
      M_indexed <- M_reduced[subject_index,]

      # Estimate the nth PDM using PDMN
      # Skip iterations that caused an eigenvalue error in fmincon

      boot_pdmn_results <-  .PDMNboot(x = x_indexed,
                                      y = y_indexed,
                                      m = M_indexed,
                                      W = redFeatWeights,
                                      initialValues = initValues)

      # Save new weights and path_coeff values after calculating
      # PDM on this bootstrapped sample
      weights_from_boot_sample <- boot_pdmn_results[['weights']]
      path_coeff_from_boot_sample <- boot_pdmn_results[['path_coeff']]

      # Save this iteration of bootstrap info
      weight_samples_k[,i] <- loading_matrix %*% weights_from_boot_sample
      path_samples_k[,i] <- path_coeff_from_boot_sample
    },
    error=function(e){})

  }

  # Remove any columns of zeros that may have happened due to iteration skipping
  weight_samples_k <- weight_samples_k[, !(colSums(weight_samples_k) == 0)]
  path_samples_k <- path_samples_k[, !(colSums(path_samples_k) == 0)]

  # Calculate the 95% CI, mean and standard dev for the bootstrapped feature weights
  for (z in 1:num_feat){
    ci <- stats::t.test(weight_samples_k[z,])
    weight_stats_k[[z]][['95CI_lb']] <- ci$conf.int[1]
    weight_stats_k[[z]][['95CI_ub']] <- ci$conf.int[2]
    weight_stats_k[[z]][['mean']] <- mean(weight_samples_k[z,])
    weight_stats_k[[z]][['sd']] <- stats::sd(weight_samples_k[z,])
    names(weight_stats_k)[z] <- sprintf('feat_%d',z)
  }

  return(list('weightStats' = weight_stats_k,
              'weightBootSamples' = weight_samples_k,
              'pathBootSamples' = path_samples_k))

}
