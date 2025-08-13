#' Summary function for mmgsem
#'
#' Summarizes the most important information of the results from a MMG-SEM model
#'
#' @param model Resulting object from a mmgsem model (must be a list). It can also be the resulting object from the ModelSelection() function.
#' @param se Must be the resulting object from the se() function of the mmgsem package. If included, the summary function will return the hypothesis testing of the relevant parameters (regressions).

#' @export
test.mmgsem <- function(model, se, multiple_comparison = FALSE) {
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
  clusters    <- paste0("K", seq_along(1:K)) # Get a string vector containing the clusters

  # Extract and define some important objects related to the SE
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

  # Perform hypothesis testing
  if(isFALSE(multiple_comparison)){ # Wald test ------------------------------------------------------------------------
    # First, do a WALD test to compare the parameters across all clusters simultaneously
    # Extract the vcov matrix of our relevant parameters
    vcov_betas <- se$vcov_betas

    # Define the restriction/contrast matrix
    n_contrast <- K - 1
    contrast <- rep(x = c(1, -1, rep(0, times = (n_contrast))),
                    times = n_contrast)
    contrast <- contrast[1:(K*n_contrast)]

    # Matrix of contrasts/ restrictions (R)
    R <- matrix(data = contrast,
                ncol = K,
                nrow = n_contrast,
                byrow = TRUE)

    # Vector containing our null value (r)
    # Null value is 0. We test if beta_1 = beta_2 as beta_1 - beta_2 = 0
    r <- rep(0, n_contrast)

    # Display results in a readable format
    cat("\nMixture Multi-Group Structural Equation Modelling (MMGSEM)\n")
    cat("Hypothesis Testing: Comparing regression parameters across clusters\n")
    cat("Wald test\n")
    cat("----------------------------------------------------------------------------\n")
    cat("\nTested Regression Coefficient (across all clusters):\n")
    # Compute the Wald statistic PER regression parameter
    for (b in 1:n_reg) {
      cat("\nParameter", regressions[b], ":\n")
      # Extract the relevant regression parameters
      betas_relevant <- betas_vector[grepl(paste0("^", regressions[b]), names(betas_vector))]
      vcov_relevant  <- vcov_betas[grepl(paste0("^", regressions[b]), colnames(vcov_betas)), grepl(paste0("^", regressions[b]), colnames(vcov_betas))]

      # Wald stat is defined as t(R %*% betas - r) %*% (R %*% cov_betas %*% t(R))^-1 %*% (R %*% betas - r)
      # Computation broken down in smaller steps

      # Differences between betas and reference value
      diff <- R %*% matrix(data = betas_relevant) - r

      # Covariance of the parameters
      V    <- MASS::ginv(R %*% vcov_relevant %*% t(R))

      # Compute the wald stat
      wald_stat <- t(diff) %*% V %*% diff
      wald_stat <- as.numeric(wald_stat) # From class matrix to numeric

      # Get the p-value
      # Degrees of freedom = number of restrictions/contrasts = 3
      df <- nrow(R)
      p_value <- pchisq(wald_stat, df, lower.tail = FALSE)

      # Create a matrix for the results
      result_df <- data.frame(
        wald_stat = round(wald_stat, 3),
        p_value   = round(p_value, 3)
      )

      # Print it
      print(result_df)
    }
  } else if(isTRUE(multiple_comparison)){ # Multiple comparison testing ------------------------------------------------
    # Prepare table for presentation of the hypothesis testing
    combinations <- t(combn(x = clusters, m = 2)) # Get all possible combinations of cluster comparisons
    combinations <- paste0(combinations[, 1], " = ", combinations[, 2])

    test_df  <- data.frame(
      Cluster_comparison = combinations,
      Difference = NA,
      z_stat     = NA,
      p_value    = NA,
      sig        = NA
    )

    # Compute the total number of comparisons
    n_tests   <- length(combinations) * n_reg
    new_alpha <- (0.05/n_tests)

    # Display results in a readable format
    cat("\nMixture Multi-Group Structural Equation Modelling (MMGSEM)\n")
    cat("Hypothesis Testing: Comparing regression parameters across clusters\n")
    cat("----------------------------------------------------------------------------\n")
    cat("\n")

    cat(n_tests, "multiple tests were computed.")
    cat("\nIf the base alpha level is 0.05, we recommend using an alpha level of", new_alpha, "(Bonferroni correction).")
    cat("\nAn * is added when the p-value is below the recommeded threshold.\n")

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
          current_comp <- paste0("K", i, " = K", j)

          # Get correct betas
          beta_k_i <- betas_relevant[grepl(paste0("k", i, "$"), names(betas_relevant))]
          beta_k_j <- betas_relevant[grepl(paste0("k", j, "$"), names(betas_relevant))]
          diff_beta <- beta_k_i - beta_k_j

          # Get correct SE
          se_k_i <- se_corrected_relevant[grepl(paste0("k", i, "$"), names(se_corrected_relevant))]
          se_k_j <- se_corrected_relevant[grepl(paste0("k", j, "$"), names(se_corrected_relevant))]
          pooled_se <- sqrt((se_k_i^2) + (se_k_j^2))

          # Z-score
          z_stat <- diff_beta/pooled_se

          #P-value
          p_value <- 2*(stats::pnorm(q = abs(z_stat), lower.tail = F))

          # Add to the table (rounded for cleaner presentation)
          test_df[test_df$Cluster_comparison == current_comp, ]$Difference <- round(diff_beta, 3)
          test_df[test_df$Cluster_comparison == current_comp, ]$z_stat     <- round(z_stat, 3)
          test_df[test_df$Cluster_comparison == current_comp, ]$p_value    <- round(p_value, 4)
          test_df[test_df$Cluster_comparison == current_comp, ]$sig        <- ifelse(p_value < new_alpha, "*", "")
        }

      }
      test_df$Cluster_comparison <- gsub(pattern = "K", replacement = "Cluster ", x = test_df$Cluster_comparison)
      print(test_df, row.names = FALSE)
    }
  }
}





















