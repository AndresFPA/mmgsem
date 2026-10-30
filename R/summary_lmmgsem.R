#' Summary for longitudinal MMG-SEM models
#'
#' @param object A fitted `lmmgsem` model or longitudinal model-selection
#'   result.
#' @param se Standard errors are not currently available for `lmmgsem` models.
#' @param model_selection If `TRUE`, print the model-selection overview.
#' @param logL Select observed or factor likelihoods for model selection.
#' @param ... Additional arguments passed to [summary.mmgsem()].
#'
#' @method summary lmmgsem
#' @export
summary.lmmgsem <- function(object, se = NULL, model_selection = FALSE,
                            logL = "observed", ...) {
  if (isTRUE(model_selection)) {
    class(object) <- "mmgsem"
    return(summary.mmgsem(
      object,
      se = se,
      model_selection = TRUE,
      logL = logL,
      ...
    ))
  }
  if (!is.null(se)) {
    stop("Standard errors are not currently implemented for lmmgsem models.",
         call. = FALSE)
  }

  beta_ks <- object$param$beta_ks
  if (!is.list(beta_ks)) {
    beta_ks <- list(beta_ks)
  }

  cat("\nLongitudinal Mixture Multigroup SEM Summary\n")
  cat("-------------------------------------------\n")
  cat("Number of groups: ", object$ngroups, "\n", sep = "")
  cat("Number of time points: ", object$ntimes, "\n", sep = "")
  cat("Number of states: ", object$nclus, "\n", sep = "")
  cat("Number of random starts: ", object$nstarts, "\n", sep = "")
  cat(
    if (isTRUE(object$converged)) "Best start converged after "
    else "Best start did not meet the convergence tolerance after ",
    object$iterations, " iterations\n",
    sep = ""
  )
  cat("Sequence log-likelihood: ", object$logLik$loglik, "\n", sep = "")

  cat("\nInitial state probabilities:\n")
  print(object$initial_probabilities)
  cat("\nTransition probabilities:\n")
  print(object$transition)

  for (k in seq_along(beta_ks)) {
    beta_df <- as.data.frame(as.table(beta_ks[[k]]))
    names(beta_df) <- c("Response", "Predictor", "estimate")
    beta_df <- beta_df[beta_df$estimate != 0, , drop = FALSE]
    beta_df$parameter <- paste(beta_df$Response, "~", beta_df$Predictor)
    beta_df <- beta_df[, c("parameter", "estimate"), drop = FALSE]
    beta_df$estimate <- round(beta_df$estimate, 3)
    cat("\nState ", k, " regression coefficients:\n", sep = "")
    print(beta_df, row.names = FALSE)
  }

  cat("\nModal state by group and time:\n")
  print(object$modal_post)
  invisible(object)
}
