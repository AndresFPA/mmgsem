#' @export
summary_MMGSEM <- function(model) {
  if (!is.list(model)) {
    stop("Input must be a list, likely the output of MMGSEM().")
  }
  
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
}