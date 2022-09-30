#' .runBootstapPDM
#'
#' @description Performs boostrap process for PDM weights (All PDMs and JointPDM).
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param M_tilde matrix (N x b) Mediator matrix with reduced dimension space of
#'  b features
#' @param W matrix (b x q) Weights for each b feature calculated for each of the
#'  q PDMs
#' @param Dt matrix (p x b)' Eigenvector space from SVD that maps M_tilde back
#' to the original Mediator matrix
#' @param WMi matrix (b x q) Initial weights used for PDM calculation optimization
#' @param Bsamp Number of bootstrap samples
#' @param whPDM Vector of indices to determine which PDM to bootstrap OR 'JointPDM'
#'
#' @return A list
#' \itemize{
#'     \item stats - 95% CI, mean, and std of bootstrap results for each
#'     feature weight for each PDM
#'     \item Wboot - Bootstrapped weights
#'     \item Tboot - Bootstrapped parameters
#'     }
#'
#' @noRd
#'

.runBootstrapPDM <- function(x, y, M_tilde, W , Dt, WMi, Bsamp, whPDM){

  ## Bootstrap individual PDMs


  if (is.numeric(whPDM)){
    cat("\nPDMs selected for bootstrap: ", whPDM,"\n",sep="\t")

    # Set empty lists to save results
    w_stats = Wboot = Tboot = vector('list', length(whPDM))

    # Loop through the PDMs to bootstrap.
    for (k in whPDM){

      cat('\n')
      cat('\n Bootstrapping PDM', k, 'with', Bsamp, 'samples')
      cat('\n')

      # First PDM
      if (k == 1){

        boot_pdm_results <- .BootPDM(x = x,
                                     y = y,
                                     M_tilde = M_tilde,
                                     W = NULL,
                                     Dt = Dt,
                                     Bsamp = Bsamp,
                                     WMi = WMi[1])

        # All other PDMs
      } else {

        boot_pdm_results <- .BootPDM(x = x,
                                     y = y,
                                     M_tilde = M_tilde,
                                     W = W[1:k-1],
                                     Dt = Dt,
                                     Bsamp = Bsamp,
                                     WMi = WMi[k])
      }

      w_stats[[k]] <- boot_pdm_results[['stats']]
      Wboot[[k]] <- boot_pdm_results[['Wboot']]
      Tboot[[k]] <- boot_pdm_results[['Tboot']]

      # Set element names
      names(w_stats)[k] <- sprintf('pdm_%d',k)
      names(Wboot)[k] <- sprintf('pdm_%d',k)
      names(Tboot)[k] <- sprintf('pdm_%d',k)
    }


    ## Joint PDM


  } else if (is.character(whPDM) && whPDM == 'jointPDM'){

    # Set variable lists
    w_stats = Wboot = Tboot = list()

    cat('\n')
    cat('\nPDM selected for bootstrap: JointPDM')

    # Bootstrap JointPDM
    boot_joint_results <- .BootPDMJoint(x = x,
                                        y = y,
                                        M_tilde = M_tilde,
                                        W = W,
                                        Dt = Dt,
                                        Bsamp = Bsamp,
                                        WMi = WMi)

    w_stats[['jpdm']] = boot_joint_results[['stats']]
    Wboot[['jpdm']] = boot_joint_results[['Wboot']]
    Tboot[['jpdm']] = boot_joint_results[['Tboot']]

  }

  return(list('stats' = w_stats, 'Wboot' = Wboot, 'Tboot' = Tboot))

}


