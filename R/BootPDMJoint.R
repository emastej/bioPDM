#' .BootPDMJoint
#'
#' @description Perform bootstrap process to calculate stats (95% CI, mean, std)
#'  on JointPDM weights
#'
#' @param x vector (N x 1) Treatment vector for N Subjects
#' @param y vector (N x 1) Outcome vector for N Subjects
#' @param M_tilde matrix (N x b) Mediator matrix with reduced dimension space of
#' b features
#' @param W matrix (b x q) Weights for each b feature calculated for each of the
#' q PDMs
#' @param Dt matrix (p x b)' Eigenvector space from SVD that maps M_tilde back
#' to the original Mediator matrix
#' @param Bsamp Number of bootstrap samples
#' @param WMi matrix (b x q) Initial weights used for PDM calculation optimization
#'
#' @return A list
#' \itemize{
#'     \item stats - 95% CI, mean, and std of bootstrap results
#'     \item Wboot - Resulting weights from bootstrap process
#'     \item Tboot - Resulting theta values (mediation path coefficients) from
#'     bootstrap process
#'     }
#'
#' @noRd
#'

.BootPDMJoint <- function(x, y, M_tilde, W, Dt, Bsamp, WMi){

  p <- dim(Dt)[2]
  nsub <- length(y)
  nPDM <- length(W)
  jWboot <- matrix(0,p,Bsamp)

  # jTboot dimensions are 4 x bootstrap sample size x nPDM
  jTboot <- array(0, dim = c(4,Bsamp, nPDM))
  jweight_stats <- vector("list", p)

  # Create matrix of randomly generated samples
  rand_samp <- matrix(0,nsub,Bsamp)

  D <- t(Dt)
  cat('\n')

  # Bootstrap samples of the data by randomly sample subjects with replacement
  for (i in 1:Bsamp){
    ransamp <- sample(1:nsub, size = nsub, replace = TRUE)
    rand_samp[,i] = ransamp
  }

  # Get the Bootstrapped Samples
  for (i in 1:Bsamp){

    # Print current bootstrapping iteration
    cat('\r', 'Bootstrap Sample:', i, '/', Bsamp)

    tryCatch(expr = {

      # Randomly sampled subjects
      subind <- rand_samp[,i]

      # Index data using bootstrap samples
      xB <- x[subind]
      yB <- y[subind]
      MB <- M_tilde[subind,]

      # Create matrices for PDM weights and Thetas
      w_n <- matrix(0, dim(M_tilde)[2], nPDM)
      theta_n <- matrix(0, 4, nPDM)

      ## Loop through the PDMs for each sample
      ## Bootstrap and calculate EACH PDM from resampled Data

      for (k in 1:nPDM){
        if (k==1){
          boot_jpdmn_result <-  .PDMNboot(x = xB,
                                          y = yB,
                                          m = MB,
                                          W = NULL,
                                          initialValues = WMi[k])
          w_n[,k] <- boot_jpdmn_result[['weights']]
          theta_n[,k] <- boot_jpdmn_result[['theta']]
        } else {
          boot_jpdmn_result <-  .PDMNboot(x  = xB,
                                          y = yB,
                                          m = MB,
                                          W = W[1:k-1],
                                          initialValues = WMi[k])
          w_n[,k] <- boot_jpdmn_result[['weights']]
          theta_n[,k] <- boot_jpdmn_result[['theta']]
        }
      }

      # Compute jointPDM
      a_n <- theta_n[3,]
      b_n <- theta_n[4,]

      # Calculate the Full Feature Weight Vector
      # [p x b] x ( [b x npdm] x [npdm x 1] ) = [p x 1]
      jWboot[,i] <- (D %*% (w_n %*% (a_n * b_n)))
      jTboot[,i,] <- theta_n

    },
    error=function(e){})
  }

  # Remove any columns of zeroes that may have happened due to iteration skipping
  jWboot = jWboot[,-(which(colSums(jWboot)==0))]
  jTboot = jTboot[,-(which(colSums(jWboot)==0)),]


  ## Calculate the 95% CI, mean and std for the bootstrapped joint feature weights


  for (z in 1:p){
    ci <- stats::t.test(jWboot[z,])
    jweight_stats[[z]][['95CI_lb']] <- ci$conf.int[1]
    jweight_stats[[z]][['95CI_ub']] <- ci$conf.int[2]
    jweight_stats[[z]][['mean']] <- mean(jWboot[z,])
    jweight_stats[[z]][['sd']] <- stats::sd(jWboot[z,])
    names(jweight_stats)[z] <- sprintf('feat_%d',z)
  }

  return(list('stats' = jweight_stats, 'Wboot' = jWboot, 'Tboot'= jTboot))

}
