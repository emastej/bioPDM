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
#' @param timeout Time allowed for optimization for each set of starting values.
#' Reducing timeout may reduce time required for calculating PDMs, but reducing
#' it too much may lead to reduced accuracy of results.
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
                     initValues,
                     numCores,
                     timeout){

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

  # Initialize clusters
  cluster <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cluster)

  # Prepare cluster environment
  parallel::clusterExport(cluster, c('x', 'y', 'M_reduced', 'redFeatWeights',
                                     'initValues', 'num_sub'),
                          envir = environment())

  # Parallelized loop
  results <- parallel::parLapply(cluster, 1:bootSamp,
                                 fun = function (i, timeout) {
  # results <- lapply(1:bootSamp, FUN = function (i, timeout) {

    # Generate random values until fmincon() returns non-imaginary
    # results
    start_user_time <- proc.time()[["elapsed"]]

    # This timer is an "overall" in case a given data set is "weird" and
    # no amount of resampling generates values which yield non-imaginary
    # results. If this fails, we warn users that they are receiving fewer BS
    # samples than they expected
    while (((proc.time()[["elapsed"]] - start_user_time) <= 300)) {

      # Bootstrap samples of the data by randomly sampling subjects with
      # replacement
      subject_index <- sample(1:num_sub, size = num_sub, replace = TRUE)

      # Index data using bootstrap samples
      x_indexed <- x[subject_index]
      y_indexed <- y[subject_index]
      M_indexed <- M_reduced[subject_index,]

      # Estimate the nth PDM using PDMN
      boot_pdmn_results <-  .PDMNboot(x = x_indexed,
                                      y = y_indexed,
                                      m = M_indexed,
                                      W = redFeatWeights,
                                      initialValues = initValues,
                                      timeout = timeout)

      if (boot_pdmn_results[['success']] == TRUE) {

        # Save new weights and path_coeff values after calculating
        # PDM on this bootstrapped sample
        weights_from_boot_sample <- boot_pdmn_results[['weights']]
        path_coeff_from_boot_sample <- boot_pdmn_results[['path_coeff']]

        # Return successful results
        return(list('weight_samples_k' = loading_matrix %*%
                      weights_from_boot_sample,
                    'path_samples_k' = path_coeff_from_boot_sample,
                    'subject_index' = subject_index,
                    'iteration' = i,
                    'success' = TRUE))

      }

    }

    # If we reach this point, it means none of the attempts within the
    # timeout limit were successful, so we return empty values along with
    # a flag indicating uncuccessful iteration
    return(list('weight_samples_k' = NA,
                'path_samples_k' = NA,
                'subject_index' = NA,
                'iteration' = i,
                'success' = FALSE))

  }, timeout = timeout)

  # Stop cluster
  parallel::stopCluster(cluster)

  # Set weight samples, path samples
  for (i in 1:bootSamp) {
    if (results[[i]][['success']] == TRUE) {
      weight_samples_k[, i] <- results[[i]][['weight_samples_k']]
      path_samples_k[, i] <- results[[i]][['path_samples_k']]
    } else {
      weight_samples_k[, i] <- rep(NA, length(weight_samples_k[,i]))
      path_samples_k[, i] <- rep(NA, length(path_samples_k[,i]))
    }
  }

  # Store vector of which iterations were successful
  successes <- lapply(results, function(item) {
    item[['success']]
  }) |>
    unlist()

  # Warn users about unsuccessful iterations
  if (sum(!successes) > 0) {
    warning(sum(!successes), " bootstrap samples were unsuccessful")
  }

  # Remove any empty that may exist due to iteration skipping
  weight_samples_k <- weight_samples_k[, successes]
  path_samples_k <- path_samples_k[, successes]

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
