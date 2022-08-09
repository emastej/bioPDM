#' .BootPDM
#'
#' @description Perform bootstrap process to calculate stats (95% CI, mean, std) on PDM weights
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param M_tilde matrix (N x b) Mediator matrix with reduced dimension space of b features
#' @param W matrix (b x q) Weights for each b feature calculated for each of the q PDMs
#' @param Dt matrix (p x b)' Eigenvector space from SVD that maps M_tilde back to the original Mediator matrix
#' @param Bsamp Number of bootstrap samples
#' @param WMi matrix (b x q) Initial weights used for PDM calculation optimization
#'
#' @return A list
#' \itemize{
#'     \item stats - 95% CI, mean, and std of bootstrap results
#'     \item Wboot - Resulting weights from bootstrap process
#'     \item Tboot - Resulting theta values (mediation path coefficients) from bootstrap process
#'     }
#'
#' @noRd
#'

.BootPDM <- function(x, y, M_tilde, W, Dt, Bsamp, WMi){

  p <- dim(Dt)[2]
  nsub <- length(y)
  # Matrix for the Weights
  kWboot <- matrix(0,p,Bsamp)
  # Matrix for the Thetas
  kTboot <- matrix(0,4,Bsamp)
  # List for resulting stats (CI, mean, std)
  kweight_stats <- vector("list", p)
  # Empty matrix of randomly generated samples
  rand_samp <- matrix(0,nsub,Bsamp)

  D <- t(Dt)
  cat('\n')

  # Bootstrap samples of the data by randomly sampling subjects with replacement
  for (i in 1:Bsamp){
    ransamp <- sample(1:nsub, size = nsub, replace = TRUE)
    rand_samp[,i] <- ransamp
  }

  # Get the bootstrapped samples
  for (i in 1:Bsamp){

    # Print current bootstrapping iteration
    cat('\r', 'Bootstrap Sample:', i, Bsamp)

    tryCatch(expr={

      # Randomly sampled subjects
      subind <- rand_samp[,i]

      # Index data using bootstrap samples
      xB <- x[subind]
      yB <- y[subind]
      MB <- M_tilde[subind,]

      # Estimate the nth PDM using PDMN
      # Skip iterations that caused an eigenvalue error in fmincon

      boot_pdmn_results <-  .PDMNboot(x = xB, y = yB, m = MB, W = W, initialValues = WMi)

      # Save new weights and theta values
      w_n <- boot_pdmn_results[['weights']]
      theta_n <- boot_pdmn_results[['theta']]

      # Save this iteration of bootstrap info to the Wboot and Tboot matrices
      kWboot[,i] <- D %*% w_n
      kTboot[,i] <- theta_n
    },
    error=function(e){})

  }

  # Remove any columns of zeros that may have happened due to iteration skipping
  kWboot = kWboot[,-(which(colSums(kWboot)==0))]
  kTboot = kTboot[,-(which(colSums(kTboot)==0))]

  # Calculate the 95% CI, mean and standard dev for the bootstrapped feature weights
  for (z in 1:p){
    ci <- stats::t.test(kWboot[z,])
    kweight_stats[[z]][['95CI_lb']] <- ci$conf.int[1]
    kweight_stats[[z]][['95CI_ub']] <- ci$conf.int[2]
    kweight_stats[[z]][['mean']] <- mean(kWboot[z,])
    kweight_stats[[z]][['sd']] <- stats::sd(kWboot[z,])
    names(kweight_stats)[z] <- sprintf('feat_%d',z)
  }

  return(list('stats' = kweight_stats, 'Wboot' = kWboot, 'Tboot' = kTboot))

}
