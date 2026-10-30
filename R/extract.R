#' Extract a model from model selection
#'
#' @description
#' Extract a model with a specified number of clusters from a model-selection
#' result created by [ModelSelection()].
#'
#' @param object A fitted object returned by [ModelSelection()].
#' @param nclus Number of clusters in the desired model.
#'
#' @return The selected `mmgsem` or `lmmgsem` model.
#'
#' @export
extract <- function(object, nclus){
  # First extract all models from the general object
  Models <- object$Models

  # Check the number of clusters of each model and extract the desired model
  for (k in seq_along(Models)) {
    nclus_model <- if (inherits(Models[[k]], "lmmgsem")) {
      Models[[k]]$nclus
    } else {
      ncol(Models[[k]]$posteriors) - 1
    }
    if(nclus_model == nclus){
      selected_model <- Models[[k]]
      class(selected_model) <- class(Models[[k]])
      return(selected_model)
    }
  }
  stop("No model with the requested number of clusters was found.",
       call. = FALSE)
}
