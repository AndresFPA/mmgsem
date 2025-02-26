#' @export
summary_MMGSEM <- function(model, model_selection = F) {
  if (!is.list(model)) {
    stop("Input must be a list, likely the output of MMGSEM().")
  }

  if(isFALSE(model_selection)){
    # Extract key components
    beta_ks     <- model$param$beta_ks     # Regression parameters
    nstarts     <- model$nstarts           # Random starts
    iterations  <- model$iterations        # Random starts
    loglik      <- model$logLik$loglik     # Factors' loglikelihood
    obs_loglik  <- model$logLik$obs_loglik # Observed loglikelihood
    BIC         <- model$model_sel$BIC     # BIC
    AIC         <- model$model_sel$AIC     # AIC
    # entropy <- model$model_sel$R2_entropy

    # Display results in a readable format
    cat("\nMixture Multi-Group Structural Equation Modelling (MMGSEM) Summary\n")
    cat("---------------------------------------------------------------------\n")
    cat("Number of random starts: ", nstarts, "\n")
    cat("Best start converged normally after", iterations, "iterations\n")

    cat("\nIf 'local' sam is used, MMG-SEM converges using the loglikelihood from the\nsecond step (only structural model parameters), here called 'Step 2'. 'Observed'\nloglikelihood is reconstructed using the measurement model parameters.\n")
    cat("\nFinal Log-Likelihood (Step 2): ", loglik, "\n")
    cat("Final Log-Likelihood (Observed): ", obs_loglik, "\n")
    cat("BIC (Step 2): ", BIC$Factors$BIC_N, " (N-based), ", BIC$Factors$BIC_G, " (G-based)\n")
    cat("BIC (Observed): ", BIC$observed$BIC_N, " (N-based), ", BIC$observed$BIC_G, " (G-based)\n")
    cat("AIC (Step 2): ", AIC$Factors, "\n")
    cat("AIC (Observed): ", AIC$observed, "\n")
    # cat("Entropy R²: ", entropy, "\n")

    # Display regression coefficients per cluster
    cat("\nRegression Coefficients (Beta Matrices):\n")
    cat("\nIn the matrix, rows are response variables and columns are predictors.\n")
    if (is.list(beta_ks)) {
      for (k in seq_along(beta_ks)) {
        cat("\nCluster", k, ":\n")
        print(beta_ks[[k]])
      }
    } else {
      print(beta_ks)
    }
  } else if (isTRUE(model_selection)){
    Overview <- model$Overview
    lower_limit <- Overview$Clusters[1]
    upper_limit <- Overview$Clusters[length(Overview$Clusters)]
    cat("\nModel Selection Summary for Mixture Multi-Group Structural Equation Modelling (MMGSEM)\n")
    cat("---------------------------------------------------------------------\n")
    cat("MMG-SEM models from", lower_limit, "to", upper_limit, "clusters were run.\n")
    cat("-----------------------------------------\n")
    cat("Convex Hull selected the model with", which.max(Overview$Chull), "clusters\n")
    cat("BIC_G selected the model with", which.min(Overview$BIC_G), "clusters\n")
    cat("BIC_N selected the model with", which.min(Overview$BIC_N), "clusters\n")
    cat("AIC selected the model with",   which.min(Overview$AIC), "clusters\n")
    cat("AIC3 selected the model with",  which.min(Overview$AIC3), "clusters\n")
    cat("ICL selected the model with",   which.min(Overview$ICL), "clusters\n")
    cat("\n")

    idx <- which(colnames(Overview) == "ICL")
    print(Overview[, 1:idx])

    cat("\n")
    cat("To check the scree plots of each model selection measure, please use the plot_mmgsem() function\n")
  }

}
