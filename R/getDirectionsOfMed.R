#' getDirectionsOfMed
#'
#' @description Calculates Principal Directions of Mediation
#' (PDMs) that mediate the effect of an exposure on an outcome.
#'
#' @param data_list List containing X, Y, M, M_reduced, tLoadingMatrix,
#' and exVar elements
#' @param nPDM Number of PDMs to calculate
#' @param doJointPDM Calculate the JointPDM. Can be True or False
#' @param doBootPDM Bootstrap samples to calculate feature weight stats
#' (95% CI, mean, standard deviation). Can be True or False
#' @param doBootJPDM Bootstrap samples to calculate feature weight stats
#' (95% CI, mean, standard deviation) for JointPDM. Can be True or False
#' @param bootSamp Number of bootstrap samples
#' @param returnBootSamples Return the PDM weights and path coefficients
#'  calculated from bootstrapped samples. Can be True or False
#' @param saveResults Save final list and a .txt file with parameters and any
#' free text notes
#' @param saveDir Directory for which the user wants to save results to
#' @param notes Parameter to write a free text note that will be saved in .txt
#' file (if saveResults = TRUE)
#' @param numCores Number of cores to use for parallel execution. Default is
#' a single core
#' @param timeout Time allowed for optimization for each set of starting values.
#' Reducing timeout may reduce time required for calculating PDMs, but reducing
#' it too much may lead to reduced accuracy of results.
#'
#' @return A list
#' \itemize{
#'     \item redFeatWeights - matrix (b x k) Weights for each reduced b feature
#'     calculated for each of the k PDMs
#'     \item featWeights - matrix (p x k) Weights mapped back onto original Mediator
#'     Matrix
#'     \item jointWeights - vector (p x 1) Weights for the JointPDM. Each Mediator
#'     feature has one weight (if applicable)
#'     \item pathCoeff - matrix (5 x k) Mediation Path Coefficients for each of the
#'     k PDMs (c, c', a, b, ab)
#'     \item PDMk - vector (p x 1) PDM vector for each k PDM
#'     \item JointPDM - vector (p x 1) JointPDM
#'     \item Boot - List of bootstrap results (if applicable)
#'     \itemize{
#'         \item featWeightStats -  95% CI, mean, and standard deviation
#'         (for each feature for each PDM) calculated from the PDM weights
#'         generated from bootstrapped samples
#'         \item featWeightSamples -  The PDM weights generated from bootstrapped
#'         samples (from which the featWeightStats were calculated)
#'         \item pathCoeffSamples - Resulting mediation pathway coefficients from
#'         each bootstrapped sample
#'         \item jointWeightStats - 95% CI, mean, and standard deviation
#'         (for each feature) calculated from JointW weights generated from
#'         bootstrapped results
#'         \item jointWeightSamples - The JointPDM weights generated from bootstrapped
#'         samples (from which the jointWeightStats were calculated)
#'     }
#'    }
#'
#' @export
#'


getDirectionsOfMed <- function(data_list = NULL,
                               nPDM = 5,
                               doJointPDM = TRUE,
                               doBootPDM = FALSE,
                               doBootJPDM = FALSE,
                               bootSamp = 1000,
                               returnBootSamples = FALSE,
                               saveResults = FALSE,
                               saveDir = NULL,
                               notes = NULL,
                               numCores = 1,
                               timeout = 5) {


  ## Check Input Data


  if (is.list(data_list)){
  } else {
      stop('Make sure your data input is a list')
  }

  # Check the structure of the data
  if (!(.checkDataStruct(data_list))){
    stop('Check the structure of input list elements')
  }

  if (saveResults == TRUE){
    if (is.null(saveDir)){
      stop("saveDir parameter is empty.
           Set the directory for which results get saved to")
    }
  }


  ## Compute PDM weights


  # Calculate the PDM weights using the runPDM function
  pdm_list <- .runPDM(x = data_list[['X']],
                      y = data_list[['Y']],
                      M = data_list[['M_reduced']],
                      tLoadingMatrix = data_list[['tLoadingMatrix']],
                      nPDM = nPDM,
                      doJointPDM = doJointPDM,
                      numCores = numCores,
                      timeout = timeout)

  # Remove Joint W from the list if JointPDM is not requested to be calculated
  if (doJointPDM == FALSE){
    pdm_list[['jointWeights']] = NULL
  }


  ## Compute PDMs


  # Calculate PDMs using the previously calculated weights
  pdm_list <- .calculatePDM(list = pdm_list,
                            M = data_list[['M']],
                            npdm = nPDM,
                            doJPDM = doJointPDM)


  ## Bootstrap Results


  # Bootstrap individual PDMs
  if (doBootPDM){

    if (all(c('initValues', 'redFeatWeights') %in% names(pdm_list)) &&
        all(c('X', 'Y', 'M_reduced', 'tLoadingMatrix') %in% names(data_list))){

      # Calculate bootstrap results using runBootstrapPDM.R
      boot_results <- .runBootstrapPDM(x = data_list[['X']],
                                       y = data_list[['Y']],
                                       M_reduced = data_list[['M_reduced']],
                                       redFeatWeights = pdm_list[['redFeatWeights']],
                                       tLoadingMatrix = data_list[['tLoadingMatrix']],
                                       initValues = pdm_list[['initValues']],
                                       bootSamp = bootSamp,
                                       whichPDM = 1:nPDM,
                                       numCores = numCores,
                                       timeout = timeout)

    } else {
      stop('Missing data to perform PDM bootstrapping')
    }

    # Collect bootstrap outputs
    pdm_list[['boot']][['featWeightStats']] <- boot_results[['weightStats']]

    if (returnBootSamples){
      pdm_list[['boot']][['featWeightSamples']] <- boot_results[['weightBootSamples']]
      pdm_list[['boot']][['pathCoeffSamples']] <- boot_results[['pathBootSamples']]
    }
  }

  # Bootstrap Joint PDM
  if (doBootJPDM){

    if (all(c('initValues', 'redFeatWeights') %in% names(pdm_list)) &&
        all(c('X', 'Y', 'M_reduced', 'tLoadingMatrix') %in% names(data_list))){

      # Calculate JointPDM boostrap results using runBootstrapPDM.R
      jboot_results <- .runBootstrapPDM(x = data_list[['X']],
                                        y = data_list[['Y']],
                                        M_reduced = data_list[['M_reduced']],
                                        redFeatWeights = pdm_list[['redFeatWeights']],
                                        tLoadingMatrix = data_list[['tLoadingMatrix']],
                                        initValues = pdm_list[['initValues']],
                                        bootSamp = bootSamp,
                                        whichPDM = 'jointPDM')

    } else {
      stop('Missing data to perform JointPDM bootstrapping')
    }

    # Collect bootstrapped outputs
    pdm_list[['boot']][['jointWeightStats']] <- jboot_results[['stats']]

    if (returnBootSamples){
      pdm_list[['boot']][['jointWeightSamples']] <- jboot_results[['weightBootSamples']]
    }
  }


  ## Save results (if saveResults = TRUE)


  if (saveResults){
    .saveOutputs(result_list = pdm_list,
                 nPDM = nPDM,
                 doBootPDM = doBootPDM,
                 doBootJPDM = doBootJPDM,
                 Bsamp = bootSamp,
                 notes = notes,
                 dir_to_save = saveDir)
  }


  # Remove initValues element in the pdm_list before returning the pdm_list
  # This element was only required for the bootstrapping step
  pdm_list$initValues <- NULL

  return(pdm_list)

}


## Helper functions
###############################################################################
# Function to calculate PDMs using the weights calculated in .runPDM() step
# Input: pdm_list, Mediator, and the nPDM
# Output: pdm_list with added elements of PDMs and JointPDM

.calculatePDM <- function(list, M, npdm, doJPDM){

  # Unlist M if neccessary
  if (is.list(M)){
    M <- matrix(unlist(M), ncol = length(M))
  }

  # Define weights
  weights <- list[['featWeights']]

  # Calculate PDM
  for (i in 1:npdm){
    pdm <- M %*% weights[[i]]
    list[[sprintf('PDM%d',i)]] <- pdm
  }

  # Calculate JointPDM
  if (doJPDM == TRUE){
    list[['JointPDM']] <- M %*% list[['jointWeights']]
  }

  return(list)
}

############################################################################
# Function to check the Data list structure before PDM calculation
# Input: Data list (should contain X, Y, M_reduced, tLoadingMatrix, and exVar)
# Output: True/ False... is the list ready for PDM calculation

.checkDataStruct <- function(DAT){

  element_names <- c('X', 'Y', 'M_reduced', 'tLoadingMatrix', 'exVar')

  # Check list element names and element type
  if (all(element_names %in% names(DAT)) && is.numeric(DAT[['X']]) &&
      is.numeric(DAT[['Y']]) && is.numeric(DAT[['M_reduced']]) &&
      is.numeric(DAT[['tLoadingMatrix']])){
  } else {
    stop("Missing fields in input structure")
  }

  # Check the dimensions of X and Y match
  if (dim(DAT[['X']])[1] == dim(DAT [['Y']])[1]){
  } else {
    stop("Data dimension mismatch between X and Y")
  }

  # Check the dimensions of X and M_reduced match
  if (dim(DAT[['X']])[1] == dim(DAT[['M_reduced']])[1]){
  } else {
    stop("Data dimension mismatch between X and M_reduced")
  }

  # Check the dimensions of tLoadingMatrix and M_reduced match
  if (dim(DAT[['tLoadingMatrix']])[1] == dim(DAT[['M_reduced']])[2]){
  } else {
    stop("Data dimension mismatch between tLoadingMatrix and M_reduced")
  }

  return(TRUE)

}

###############################################################################

# Function to save results from getDirectionsOfMed()
# Saves the final output list and a notes .txt file

.saveOutputs <- function(result_list, nPDM, doBootPDM, doBootJPDM, Bsamp,
                         notes, dir_to_save){

  # Create a folder to save results
  out_Dir <- "radPDM_results"

  # Create the radPDM_results folder if it does not already exist in the directory
  if (!file.exists(paste0(dir_to_save,'/',out_Dir))){
    dir.create(paste0(dir_to_save,'/',out_Dir))
  }

  # Create a new folder with current Time and Date for every iteration
  current_folder <- .currentTime()
  current_dir <- paste0(dir_to_save,'/',out_Dir)
  dir.create(file.path(current_dir, current_folder))

  # Save the final output list
  save(result_list, file = paste0(current_dir,'/',current_folder,'/', 'pdm_results.rda'))

  # Check if bootstrapping occurred
  if (!(doBootJPDM | doBootPDM)){
    Bsamp <- NULL
  }

  # Save out .text file that includes time, nPDM, Bsamp (if applicable),
  # and notes (if applicable)
  txt <- c('Time:', current_folder, '\nnPDM:', nPDM, '\nBootstrapSamples:',
           Bsamp, '\nNotes:', notes)
  writeLines(txt, paste0(current_dir,'/',current_folder,'/',
                         'Parameters and Notes.txt'), sep = '\t')

  }

###############################################################################

# Function to set the Time and Date in order to define a new time stamped dir
# for saving results to

.currentTime <- function(time = Sys.time()){
  formatted_time = format(time, "%F %H-%M")
  return(formatted_time)
}




