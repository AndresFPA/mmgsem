#' Summary function for mmgsem
#'
#' Summarizes the most important information of the results from a MMG-SEM model
#'
#' @param model Resulting object from a mmgsem model (must be a list). It can also be the resulting object from the ModelSelection() function.
#' @param se Must be the resulting object from the se() function of the mmgsem package. If included, the summary function will return the hypothesis testing of the relevant parameters (regressions).

#' @export
summary_MMGSEM <- function(model, se = NULL, model_selection = F) {
  if (!is.list(model)) {
    stop("Input must be a list, likely the output of MMGSEM() or ModelSelection().")
  }

  if (!is.null(se) & isTRUE(model_selection)){
    stop("se argument can only be used with a single resulting model from MMGSEM(). Thus, is not compatible model_selection = T")
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
    posteriors  <- model$posteriors        # Posterior matrix

    if(!is.null(se)){
      # Extract standard errors
      betas_se           <- se$SE_vector$betas_se
      betas_se_corrected <- se$SE_vector$betas_se_corrected
      # Extract betas as vector
      betas_vector       <- se$estimates_vector$betas_est

      # Prepare z-scores and p-values
      betas_z           <- betas_vector/betas_se
      betas_z_corrected <- betas_vector/betas_se_corrected

      # For the p-values:
      # Use the pnorm() function with absolute z scores and lower.tail = F to obtain correct TWO-TAILED p-values.
      # They must be multiplied by 2.
      betas_pvalue           <- 2*(stats::pnorm(q = abs(betas_z), lower.tail = F))
      betas_pvalue_corrected <- 2*(stats::pnorm(q = abs(betas_z_corrected), lower.tail = F))
    }

    # entropy <- model$model_sel$R2_entropy

    # Display results in a readable format
    cat("\nMixture Multi-Group Structural Equation Modelling (MMGSEM) Summary\n")
    cat("---------------------------------------------------------------------\n")
    cat("Number of clusters: ", (ncol(posteriors)-1), "\n")
    cat("Number of random starts: ", nstarts, "\n")
    cat("Best start converged normally after", iterations, "iterations\n")

    cat("\nIf 'local' sam is used, MMG-SEM converges using the loglikelihood from the\nsecond step (only structural model parameters), here called 'Step 2' Log-Likelihood.\n 'Observed'loglikelihood is reconstructed using the measurement model parameters.\n")
    cat("\nFinal Log-Likelihood (Step 2): ", loglik, "\n")
    cat("Final Log-Likelihood (Observed): ", obs_loglik, "\n")
    cat("BIC (Step 2): ", BIC$Factors$BIC_N, " (N-based), ", BIC$Factors$BIC_G, " (G-based)\n")
    cat("BIC (Observed): ", BIC$observed$BIC_N, " (N-based), ", BIC$observed$BIC_G, " (G-based)\n")
    cat("AIC (Step 2): ", AIC$Factors, "\n")
    cat("AIC (Observed): ", AIC$observed, "\n")
    # cat("Entropy R²: ", entropy, "\n")

    # Display regression coefficients per cluster
    cat("\nRegression Coefficients:\n")
    for (k in seq_along(beta_ks)) {
      # Prepare the output
      # Convert matrix to a more suitable data frame
      if (is.list(beta_ks)) {
        beta_df <- as.data.frame(as.table(model$param$beta_ks[[k]]))
      } else { # Which means we only have one cluster
        beta_df <- as.data.frame(as.table(model$param$beta_ks))
      }

      # Rename columns
      colnames(beta_df) <- c("Response", "Predictor", "Estimate")

      # Remove rows with zero estimates
      beta_df <- beta_df[beta_df$Estimate != 0, ]

      # Convert to formula-like representation
      beta_df$Parameter <- paste0(beta_df$Response, " ~ ", beta_df$Predictor)

      # Select relevant columns
      beta_df <- beta_df[, c("Parameter", "Estimate")]

      # Round beta results for cleaner presentation
      beta_df$Estimate <- round(beta_df$Estimate, 3)

      # Add standard errors if requested
      if(!is.null(se)){
        # Extract the relevant values for this cluster
        pattern <- paste0("k", k, "$")
        betas_se_k               <- betas_se[grepl(pattern, names(betas_se))]
        betas_se_corrected_k     <- betas_se_corrected[grepl(pattern, names(betas_se_corrected))]
        betas_z_k                <- betas_z[grepl(pattern, names(betas_z))]
        betas_pvalue_k           <- betas_pvalue[grepl(pattern, names(betas_pvalue))]
        betas_z_corrected_k      <- betas_z_corrected[grepl(pattern, names(betas_z_corrected))]
        betas_pvalue_corrected_k <- betas_pvalue_corrected[grepl(pattern, names(betas_pvalue_corrected))]

        # Round for cleaner presentation
        betas_se_k               <- round(betas_se_k, 3)
        betas_se_corrected_k     <- round(betas_se_corrected_k, 3)
        betas_z_k                <- round(betas_z_k, 3)
        betas_pvalue_k           <- round(betas_pvalue_k, 3)
        betas_z_corrected_k      <- round(betas_z_corrected_k, 3)
        betas_pvalue_corrected_k <- round(betas_pvalue_corrected_k, 3)

        # Original hypothesis testing (without correction) is not printed
        # beta_df$se      <- betas_se_k
        # beta_df$z_score <- betas_z_k
        # beta_df$p_value <- betas_pvalue_k

        beta_df$se      <- betas_se_corrected_k
        beta_df$z_score <- betas_z_corrected_k
        beta_df$p_value <- betas_pvalue_corrected_k
      }

      # Print results
      cat("\nCluster", k, ":\n")
      print(beta_df, row.names = FALSE)
    }

    # Print the (modal) clustering results. Posterior matrix is avoided due to possible large results
    cat("\nClustering results:\n")
    # Transform soft clustering to hard clustering
    post_no_group <- posteriors[, 1:(ncol(posteriors) - 1)]
    Modal_posteriors <- t(apply(post_no_group, 1, function(x) as.numeric(x == max(x))))
    Modal_posteriors <- as.data.frame(Modal_posteriors)
    Modal_posteriors$Group <- posteriors$Group

    for(k in 1:(ncol(posteriors)-1)){
      # browser()
      members_k <- Modal_posteriors$Group[Modal_posteriors[, k] == 1]
      cat("\nCluster", k, ":\n")
      cat(members_k)
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
