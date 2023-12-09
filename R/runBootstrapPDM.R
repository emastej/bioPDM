#' .runBootstapPDM
#'
#' @description Performs boostrap process for PDM weights (All PDMs and JointPDM).
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param M_reduced matrix (N x b) Mediator matrix with reduced dimension space of
#'  b features
#' @param redFeatWeights matrix (b x k) Weights for each b feature calculated for each of the
#'  k PDMs
#' @param tLoadingMatrix matrix (p x b)' Eigenvector space from SVD that maps
#' M_reduced back to the original Mediator matrix
#' @param initValues matrix (b x k) Initial weights used for PDM calculation optimization
#' @param bootSamp Number of bootstrap samples
#' @param whichPDM Vector of indices to determine which PDM to bootstrap OR 'JointPDM'
#'
#' @return A list
#' \itemize{
#'     \item weightStats - 95% CI, mean, and std of bootstrap results for each
#'     feature weight for each PDM
#'     \item weightBootSamples - Bootstrapped weight samples
#'     \item pathBootSamples - Bootstrapped mediation path coefficients
#'     }
#'
#' @noRd
#'

.runBootstrapPDM <- function(x,
                             y,
                             M_reduced,
                             redFeatWeights,
                             tLoadingMatrix,
                             initValues,
                             bootSamp,
                             whichPDM,
                             num,
                             numCores,
                             timeout){

  ## Bootstrap individual PDMs


  if (is.numeric(whichPDM)){
    cat("\nPDMs selected for bootstrap: ", whichPDM,"\n",sep="\t")

    # Set empty lists to save results
    weight_stats = weight_samples = path_coeff_samples = vector('list', length(whichPDM))

    # Loop through the PDMs to bootstrap.
    for (k in whichPDM){

      cat('\n')
      cat('\n Bootstrapping PDM', k, 'with', bootSamp, 'samples')
      cat('\n')

      # First PDM
      if (k == 1){

        boot_pdm_results <- .BootPDM(x = x,
                                     y = y,
                                     M_reduced = M_reduced,
                                     redFeatWeights = NULL,
                                     tLoadingMatrix = tLoadingMatrix,
                                     bootSamp = bootSamp,
                                     initValues = initValues[1],
                                     numCores = numCores,
                                     timeout = timeout)

        # All other PDMs
      } else {

        boot_pdm_results <- .BootPDM(x = x,
                                     y = y,
                                     M_reduced = M_reduced,
                                     redFeatWeights = redFeatWeights[1:k-1],
                                     tLoadingMatrix = tLoadingMatrix,
                                     bootSamp = bootSamp,
                                     initValues = initValues[k],
                                     numCores = numCores,
                                     timeout = timeout)
      }

      weight_stats[[k]] <- boot_pdm_results[['weightStats']]
      weight_samples[[k]] <- boot_pdm_results[['weightBootSamples']]
      path_coeff_samples[[k]] <- boot_pdm_results[['pathBootSamples']]

      # Set element names
      names(weight_stats)[k] <- sprintf('pdm_%d',k)
      names(weight_samples)[k] <- sprintf('pdm_%d',k)
      names(path_coeff_samples)[k] <- sprintf('pdm_%d',k)
    }


    ## Joint PDM


  } else if (is.character(whichPDM) && whichPDM == 'jointPDM'){

    # Set variable lists
    weight_stats = weight_samples = path_coeff_samples = list()

    cat('\n')
    cat('\nPDM selected for bootstrap: JointPDM')

    # Bootstrap JointPDM
    boot_joint_results <- .BootPDMJoint(x = x,
                                        y = y,
                                        M_reduced = M_reduced,
                                        redFeatWeights = redFeatWeights,
                                        tLoadingMatrix = tLoadingMatrix,
                                        bootSamp = bootSamp,
                                        initValues = initValues)

    weight_stats[['jpdm']] = boot_joint_results[['weightStats']]
    weight_samples[['jpdm']] = boot_joint_results[['weightBootSamples']]
    path_coeff_samples[['jpdm']] = boot_joint_results[['pathBootSamples']]

  }

  return(list('weightStats' = weight_stats, 'weightBootSamples' = weight_samples,
              'pathBootSamples' = path_coeff_samples))

}


