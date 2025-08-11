#' Model extraction function for MMGSEM
#'
#' @description
#' Function used to extract an specific model out of a Model Selection (mmgsem) fitted object.
#'
#' @usage extract(object)
#'
#' @param object A resulting fitted object from the ModelSelection function.
#' @param nclus The number of clusters of the desired model
#'
#' @return selected_model: the selected model
#'
#' @export
extract <- function(object, nclus){
  # First extract all models from the general object
  Models <- modelSelection_fit$Models

  # Check the number of clusters of each model and extract the desired model
  for(k in 1:length(Models)){
    nclus_model <- ncol(Models[[k]]$posteriors) - 1 # Number of clusters
    if(nclus_model == nclus){
      selected_model <- Models[[k]]
      class(selected_model) <- "mmgsem"
      return(selected_model)
    }
  }
}
