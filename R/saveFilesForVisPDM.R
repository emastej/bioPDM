#' saveFilesForVisPDM
#'
#' @description Saves .csv files that can be uploaded into visPDM app. 
#'
#' @param data_list List containing X, Y, M, M_reduced, tLoadingMatrix,
#' and exVar elements. Output from reduceMedDimension()
#' @param pdm_result_list List output from getDirectionsOfMed()
#' @param saveDir Directory for which the user wants to save results to
#' 
#'
#' @return .csv files. *PDMs_X_Y.csv files contain the PDM, X, and Y values for 
#' each subject. *PDM_feature_loading.csv files contain the feature weights for 
#' each PDM
#' 
#' @importFrom utils write.csv
#'
#' @export


saveFilesForVisPDM <- function(data_list, pdm_result_list, saveDir = NULL){
  
  # Check that user defined a directory
  if (is.null(saveDir)){
    stop('Please define directory of which to save files')
  }
  
  
  # Save the PDM files - non sparse 
  
  
  # Save PDM, X, Y
  PDM_df <- data.frame(pdm_result_list$PDMs)
  PDM_df$X <- data_list$X[,1]
  PDM_df$Y <- data_list$Y[,1]
  utils::write.csv(PDM_df, paste0(saveDir, "/PDMs_X_Y.csv"), row.names = FALSE)

  # Save the feature loadings for each PDM
  weights_df <- data.frame(pdm_result_list$featWeights)
  colnames(weights_df) <- names(pdm_result_list$PDMs)
  rownames(weights_df) <- names(data_list$M)
  utils::write.csv(weights_df, paste0(saveDir, "/PDM_feature_loadings.csv"), row.names = TRUE)
  
  
  # Save sparse PDM files (if sparsity was calculated)
  
  
  if ('sparseThresh' %in% names(pdm_result_list)){
    # Save PDM, X, Y
    spPDM_df <- data.frame(pdm_result_list$sparseThresh$spPDM)
    spPDM_df$X <- data_list$X[,1]
    spPDM_df$Y <- data_list$Y[,1]
    utils::write.csv(spPDM_df, paste0(saveDir, "/spPDMThresh_X_Y.csv"), row.names = FALSE)
    
    # Save the feature loadings for each PDM
    spweights_df <- data.frame(pdm_result_list$sparseThresh$spFeatWeights)
    spweights_df[spweights_df == 0] <- NA
    colnames(spweights_df) <- names(pdm_result_list$PDMs)
    rownames(spweights_df) <- names(data_list$M)
    utils::write.csv(spweights_df, paste0(saveDir, "/spPDMThresh_feature_loadings.csv"), row.names = TRUE)
 
  }
  
  if ('sparseEN' %in% names(pdm_result_list)){
    # Save PDM, X, Y
    spePDM_df <- data.frame(pdm_result_list$sparseEN$spPDM)
    spePDM_df$X <- data_list$X[,1]
    spePDM_df$Y <- data_list$Y[,1]
    utils::write.csv(spePDM_df, paste0(saveDir, "/spPDMElasticNet_X_Y.csv"), row.names = FALSE)
    
    # Save the feature loadings for each PDM
    speweights_df <- data.frame(pdm_result_list$sparseEN$spFeatWeights)
    speweights_df[speweights_df == 0] <- NA
    colnames(speweights_df) <- names(pdm_result_list$PDMs)
    rownames(speweights_df) <- names(data_list$M)
    utils::write.csv(speweights_df, paste0(saveDir, "/spPDMElasticNet_feature_loadings.csv"), row.names = TRUE)
    
  }
  
}


