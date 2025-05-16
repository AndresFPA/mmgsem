#' Summary function for mmgsem
#'
#' Summarizes the most important information of the results from a MMG-SEM model
#'
#' @param model Resulting object from a mmgsem model (must be a list). It can also be the resulting object from the ModelSelection() function.
#' @param se Must be the resulting object from the se() function of the mmgsem package. If included, the summary function will return the hypothesis testing of the relevant parameters (regressions).

#' @export
test_MMGSEM <- function(model, se) {
  if (!is.list(model)) {
    stop("Input must be a list, likely the output of MMGSEM() or ModelSelection().")
  }

  # Extract key components
  beta_ks     <- model$param$beta_ks     # Regression parameters
  nstarts     <- model$nstarts           # Random starts
  iterations  <- model$iterations        # Random starts
  loglik      <- model$logLik$loglik     # Factors' loglikelihood
  obs_loglik  <- model$logLik$obs_loglik # Observed loglikelihood
  BIC         <- model$model_sel$BIC     # BIC
  AIC         <- model$model_sel$AIC     # AIC
  posteriors  <- model$posteriors        # Posterior matrix
  K           <- length(beta_ks)         # Number of clusters

  # Prepare table for presentation of the hypothesis testing
  clusters <- paste0("K", seq_along(1:K))       # Get a string vector containing the clusters
  combinations <- t(combn(x = clusters, m = 2)) # Get all possible combinations of cluster comparisons
  combinations <- paste0(combinations[, 1], "=", combinations[, 2])

  test_df  <- data.frame(
    Cluster_comparison = combinations,
    Difference = NA,
    z_value    = NA,
    p_value    = NA
  )

  # Get standard errors and betas
  # Extract standard errors
  betas_se           <- se$SE_vector$betas_se
  betas_se_corrected <- se$SE_vector$betas_se_corrected

  # Extract betas as vector
  betas_vector       <- se$estimates_vector$betas_est

  # Identify how many regression parameters we have
  regressions        <- betas_se_corrected[grepl(paste0("k1$"), names(betas_se_corrected))]
  n_reg              <- length(regressions)

  # Get a "clean" object of regression parameters (without the .k)
  regressions        <- sub(pattern = "\\..*", replacement = "", x = names(regressions))

  # Display results in a readable format
  cat("\nMixture Multi-Group Structural Equation Modelling (MMGSEM)\n")
  cat("Hypothesis Testing: Comparing regression parameters across clusters\n")
  cat("----------------------------------------------------------------------------\n")

  # Display regression coefficients per cluster
  cat("\nTested Regression Coefficient (cluster comparison):\n")
  for (b in 1:n_reg) {
    # Print relevant regression parameter
    cat("\nParameter", regressions[b], ":\n")

    # Extract the relevant parameters
    betas_relevant        <- betas_vector[grepl(paste0("^", regressions[b]), names(betas_vector))]
    se_corrected_relevant <- betas_se_corrected[grepl(paste0("^", regressions[b]), names(betas_se_corrected))]

    test_df$Cluster_comparison <- combinations

    for(i in 1:K){
      for(j in 1:K){
        if(i == j | i > j){next}

        # Which comparison?
        current_comp <- paste0("K", i, "=K", j)

        # Get correct betas
        beta_k_i <- betas_relevant[grepl(paste0("k", i, "$"), names(betas_relevant))]
        beta_k_j <- betas_relevant[grepl(paste0("k", j, "$"), names(betas_relevant))]
        diff_beta <- beta_k_i - beta_k_j

        # Get correct SE
        se_k_i <- se_corrected_relevant[grepl(paste0("k", i, "$"), names(se_corrected_relevant))]
        se_k_j <- se_corrected_relevant[grepl(paste0("k", j, "$"), names(se_corrected_relevant))]
        pooled_se <- sqrt((se_k_i^2) + (se_k_j^2))

        # Z-score
        z_score <- diff_beta/pooled_se

        #P-value
        p_value <- 2*(stats::pnorm(q = abs(z_score), lower.tail = F))

        # Add to the table (rounded for cleaner presentation)
        test_df[test_df$Cluster_comparison == current_comp, ]$Difference <- round(diff_beta, 3)
        test_df[test_df$Cluster_comparison == current_comp, ]$z_value    <- round(z_score, 3)
        test_df[test_df$Cluster_comparison == current_comp, ]$p_value    <- round(p_value, 3)
      }

    }
    test_df$Cluster_comparison <- gsub(pattern = "K", replacement = "Cluster ", x = test_df$Cluster_comparison)
    print(test_df, row.names = FALSE)
  }
}





















