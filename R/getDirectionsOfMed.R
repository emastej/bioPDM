#' getDirectionsOfMed
#'
#' @description Calculates radiomic-based Principal Directions of Mediation
#' (PDMs) that mediate the effect of an exposure on an outcome.
#'
#' @param data_list List containing X, Y, M, M_tilde, Dt, and exVar elements
#' @param nPDM Number of PDMs to calculate
#' @param doJointPDM Calculate the JointPDM. Can be True or False
#' @param doBootPDM Bootstrap samples to calculate radiomic feature weight stats
#' (95% CI, mean, standard deviation). Can be True or False.
#' @param doBootJPDM Bootstrap samples to calculate feature weight stats
#' (95% CI, mean, standard deviation) for JointPDM. Can be True or False
#' @param BootSamp Number of bootstrap samples
#' @param returnBootsamples Return the PDM weights and theta calculated from
#' bootstrapped samples. Can be True or False.
#' @param saveResults Save final list and a .txt file with parameters and any
#' free text notes
#' @param saveDir Directory for which the user wants to save results to
#' @param notes Parameter to write a free text note that will be saved in .txt
#' file (if saveResults = TRUE)
#'
#' @return A list
#' \itemize{
#'     \item w_k - matrix (b x q) Weights for each b feature calculated for each
#'      of the q PDMs
#'     \item w_k*Dt - matrix (p x q) Weights mapped back onto original Mediator
#'     Matrix
#'     \item JointW - vector (p x 1) Weights for the JointPDM. Each Mediator
#'     feature has one weight (if applicable)
#'     \item Theta - matrix (5 x q) Mediation Path Coefficients for each of the
#'     q PDMs (c, c', a, b, ab)
#'     \item PDMk - vector (p x 1) PDM vector for each k PDM
#'     \item JointPDM - vector (p x 1) JointPDM
#'     \item Boot - List of bootstrap results (if applicable)
#'     \itemize{
#'         \item w_kStats -  95% CI, mean, and standard deviation
#'         (for each feature for each PDM) calculated from the PDM weights
#'         generated from bootstrapped samples
#'         \item w_kSamples -  The PDM weights generated from bootstrapped
#'         samples (from which the w_kStats were calculated)
#'         \item ThetaSamples - Resulting mediation pathway coefficients from
#'         each bootstrapped sample
#'         \item JointWStats - 95% CI, mean, and standard deviation
#'         (for each feature) calculated from JointW weights generated from
#'         bootstrapped results
#'         \item JointWSamples - The JointPDM weights generated from bootstrapped
#'         samples (from which the JointWStats were calculated)
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
                                  BootSamp = 1000,
                                  returnBootsamples = FALSE,
                                  saveResults = FALSE,
                                  saveDir = NULL,
                                  notes = NULL){


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
                      M_tilde = data_list[['M_tilde']],
                      Dt = data_list[['Dt']],
                      nPDM = nPDM,
                      doJointPDM = doJointPDM)

  # Remove Joint W from the list if JointPDM is not requested to be calculated
  if (doJointPDM == FALSE){
    pdm_list[['JointW']] = NULL
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

    if (all(c('WMinit', 'w_k') %in% names(pdm_list)) &&
        all(c('X', 'Y', 'M_tilde', 'Dt') %in% names(data_list))){

      # Calculate bootstrap results using runBootstrapPDM.R
      boot_results <- .runBootstrapPDM(x = data_list[['X']],
                                       y = data_list[['Y']],
                                       M_tilde = data_list[['M_tilde']],
                                       W = pdm_list[['w_k']],
                                       Dt = data_list[['Dt']],
                                       WMi = pdm_list[['WMinit']],
                                       Bsamp = BootSamp,
                                       whPDM = 1:nPDM)

    } else {
      stop('Missing data to perform PDM bootstrapping')
    }

    # Collect bootstrap outputs
    pdm_list[['boot']][['w_kStats']] <- boot_results[['stats']]

    if (returnBootsamples){
      pdm_list[['boot']][['w_kSamples']] <- boot_results[['Wboot']]
      pdm_list[['boot']][['ThetaSamples']] <- boot_results[['Tboot']]
    }
  }

  # Bootstrap Joint PDM
  if (doBootJPDM){

    if (all(c('WMinit', 'w_k') %in% names(pdm_list)) &&
        all(c('X', 'Y', 'M_tilde', 'Dt') %in% names(data_list))){

      # Calculate JointPDM boostrap results using runBootstrapPDM.R
      jboot_results <- .runBootstrapPDM(x = data_list[['X']],
                                        y = data_list[['Y']],
                                        M_tilde = data_list[['M_tilde']],
                                        W = pdm_list[['w_k']],
                                        Dt = data_list[['Dt']],
                                        WMi = pdm_list[['WMinit']],
                                        Bsamp = BootSamp,
                                        whPDM = 'jointPDM')

    } else {
      stop('Missing data to perform JointPDM bootstrapping')
    }

    # Collect bootstrapped outputs
    pdm_list[['boot']][['JointWStats']] <- jboot_results[['stats']]

    if (returnBootsamples){
      pdm_list[['boot']][['JointWSamples']] <- jboot_results[['Wboot']]
    }
  }


  ## Save results (if saveResults = TRUE)


  if (saveResults){
    .saveOutputs(result_list = pdm_list,
                 nPDM = nPDM,
                 doBootPDM = doBootPDM,
                 doBootJPDM = doBootJPDM,
                 Bsamp = BootSamp,
                 notes = notes,
                 dir_to_save = saveDir)
  }


  # Remove WMinit element in the pdm_list before returning the pdm_list
  # This element was only required for the bootstrapping step
  pdm_list$WMinit <- NULL

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
  weights <- list[['w_k*Dt']]

  # Calculate PDM
  for (i in 1:npdm){
    pdm <- M %*% weights[[i]]
    list[[sprintf('PDM%d',i)]] <- pdm
  }

  # Calculate JointPDM
  if (doJPDM == TRUE){
    list[['JointPDM']] <- M %*% list[['JointW']]
  }

  return(list)
}

############################################################################
# Function to check the Data list structure before PDM calculation
# Input: Data list (should contain X, Y, M_tilde, Dt, and exvar)
# Output: True/ False... is the list ready for PDM calculation

.checkDataStruct <- function(DAT){

  element_names <- c('X', 'Y', 'M_tilde', 'Dt', 'exVar')

  # Check list element names and element type
  if (all(element_names %in% names(DAT)) && is.numeric(DAT[['X']]) &&
      is.numeric(DAT[['Y']]) && is.numeric(DAT[['M_tilde']]) &&
      is.numeric(DAT[['Dt']])){
  } else {
    stop("Missing fields in input structure")
  }

  # Check the dimensions of X and Y match
  if (dim(DAT[['X']])[1] == dim(DAT [['Y']])[1]){
  } else {
    stop("Data dimension mismatch between X and Y")
  }

  # Check the dimensions of X and M_tilde match
  if (dim(DAT[['X']])[1] == dim(DAT[['M_tilde']])[1]){
  } else {
    stop("Data dimension mismatch between X and M_tilde")
  }

  # Check the dimensions of Dt and M_tilde match
  if (dim(DAT[['Dt']])[1] == dim(DAT[['M_tilde']])[2]){
  } else {
    stop("Data dimension mismatch between Dt and M_tilde")
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




