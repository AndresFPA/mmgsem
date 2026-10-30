sem_est_fast_long <- function(S2,
                         post,
                         cov_eta,
                         nstates,
                         ntimes,
                         ngroups,
                         N_gs,
                         ite,
                         sem_par){
  # For consistency
  i <- ite

  # Get necessary parameters for iterations where i > 1
  if(i > 1){
    beta_ks <- sem_par$beta_ks
    psi_gks <- sem_par$psi_gks
  } else {
    # Initialize psi_gks
    psi_gks <- matrix(data = list(NA), nrow = (ngroups*ntimes), ncol = nstates)
  }

  # Swap back the dimensions of the posterior matrices for easier manipulation
  post <- aperm(a = post, perm = c(1,3,2))

  # Necessary objects
  nclus <- nstates
  gro_clu <- nclus*ngroups*ntimes
  # browser()
  # ----------------------------------------------------------------------------
  # PREPARE TRICK TO SPEED-UP THE ESTIMATION
  # ----------------------------------------------------------------------------

  # Do a fake sem() to obtain the correct settings to use in Step 2
  # The fake sem object determines which model we actually fit. All important setting must be done here!
  # just a single sample cov!
  fake <- lavaan::sem(
    model              = S2,
    sample.cov         = rep(cov_eta[1], nclus),
    sample.nobs        = rep(sum(N_gs), nclus),
    do.fit             = FALSE,
    baseline           = FALSE,
    h1                 = FALSE,
    check.post         = FALSE,
    loglik             = FALSE,
    sample.cov.rescale = FALSE,
    fixed.x            = TRUE
  )

  FakeprTbl           <- lavaan::parTable(fake)
  fake@Options$do.fit <- TRUE
  fake@Options$se     <- "none"
  fake@ParTable$start <- NULL
  fake@ParTable$est   <- NULL
  fake@ParTable$se    <- NULL
  fake@Options$start  <- "default"

  # Get the labels of the endogenous 1 and 2 factors
  lat_var <- lavNames(lavaanify(S2, auto = TRUE), "ov")
  endog1 <- lat_var[(lat_var %in% FakeprTbl$rhs[which(FakeprTbl$op == "~")]) &
                      (lat_var %in% FakeprTbl$lhs[which(FakeprTbl$op == "~")])]
  endog2 <- lat_var[!c(lat_var %in% FakeprTbl$rhs[which(FakeprTbl$op == "~")]) &
                      (lat_var %in% FakeprTbl$lhs[which(FakeprTbl$op == "~")])]
  endog <- c(endog1, endog2)
  exog  <- lat_var[!c(lat_var %in% endog)]


  # Do a fake model per endo LV (to avoid bias due to reconstruction of group-specific endo variances).
  # Please note that the index "lv" is used to identify each model per endo LV
  # Please also note that this is used AFTER the first iteration. In the first iteration we start with only one model
  fake_lv  <- vector(mode = "list", length = length(endog))
  prTbl_lv <- vector(mode = "list", length = length(endog))

  for (lv in 1:length(endog)) {
    # Create a parameter table per endogenous latent variables
    # Select the current latent variable
    this_lv <- endog[endog %in% endog[lv]]
    # browser()
    # Keep the (co)variances of the other latent variables and (if global) the measurement parameters
    var_not_this_lv <- which(FakeprTbl$lhs != this_lv & FakeprTbl$op == "~~") # keeps lv covariances and obs res variances
    # fac_load        <- which(FakeprTbl$lhs != this_lv & FakeprTbl$op == "=~") # Factor load of other factors

    # Get the new parameter table per endogenous latent variable
    prTbl_idx       <- c(which(FakeprTbl$lhs == this_lv), var_not_this_lv)
    prTbl_idx       <- sort(prTbl_idx)
    prTbl_lv[[lv]]  <- FakeprTbl[prTbl_idx, ]

    # Run the model per endo latent variable
    fake_lv[[lv]] <- lavaan::sem(
      model              = prTbl_lv[[lv]],
      sample.cov         = rep(cov_eta[1], nclus),
      sample.nobs        = rep(sum(N_gs), nclus),
      do.fit             = FALSE,
      baseline           = FALSE,
      h1                 = FALSE,
      check.post         = FALSE,
      loglik             = FALSE,
      sample.cov.rescale = FALSE,
      fixed.x            = TRUE
    )

    fake_lv[[lv]]@Options$do.fit <- TRUE
    fake_lv[[lv]]@Options$se     <- "none"
    fake_lv[[lv]]@ParTable$start <- NULL
    fake_lv[[lv]]@ParTable$est   <- NULL
    fake_lv[[lv]]@ParTable$se    <- NULL
    fake_lv[[lv]]@Options$start  <- "default"
  }

  # Re-order (order of columns and rows) cov_eta to make sure later computations are comparing correct matrices
  cov_eta <- lapply(1:(ngroups*ntimes), function(x) {
    reorder_long(cov_eta[[x]], exog = exog, endog = endog)
  })

  # ----------------------------------------------------------------------------
  # PARAMETER ESTIMATION
  # ----------------------------------------------------------------------------
  # Create an stacked posterior matrix to estimate the structural parameters ACROSS time points
  # This ensures that the cluster labels remain the same across time (i.e., avoid label switching)
  post_stacked <- do.call(rbind, lapply(seq_len(ntimes), function(t) {
    matrix(post[, , t], nrow = ngroups, ncol = nclus) # use matrix() to ensure correct dimensions even when k = 1
  }))
  # post_stacked <- c()
  # for(t in 1:ntimes){
  #   post_stacked <- rbind(post_stacked, post[, , t])
  # }
  # browser()
  # Obtain weighted group-cluster combination sample size
  N_gks <- post_stacked * N_gs
  N_gks <- c(N_gks)

  # Initialize weighted posteriors (weighted for each endo LV)
  # Done to avoid bias - NOT NECESSARY in the first iteration
  post_lv <- vector(mode = "list", length = length(endog))
  for (lv in 1:length(endog)) {
    post_lv[[lv]] <- matrix(data = NA, ncol = nclus, nrow = (ngroups*ntimes))
  }

  # Weight posteriors for each endogenous factor - Not necessary in the first iteration
  # To avoid bias when endo_group_specific == T. That is, when the endogenous variances are group-specific
  if (i > 1) {
    for (lv in 1:length(endog)) {
      for (k in 1:nclus) {
        for (g in 1:(ngroups*ntimes)) {
          # Correct bias by dividing by the correct endo LV
          post_lv[[lv]][g, k] <- post_stacked[g, k] / psi_gks[[g, k]][endog[lv], endog[lv]]
        }
      }
    }
  }

  # Trick to avoid slow multi-group estimation
  # Get a weighted averaged covariance matrix for each cluster
  if (i == 1) {
    # For the first iteration there is no weighted posteriors
    COV <- vector("list", length = nclus)

    for (k in 1:nclus) {
      # create 'averaged' sample cov for this cluster
      this_nobs <- post_stacked[, k] * N_gs
      this_w <- this_nobs / sum(this_nobs)
      tmp <- lapply(seq_along(cov_eta), function(g) {
        cov_eta[[g]] * this_w[g]
      })
      COV[[k]] <- Reduce("+", tmp)
    }
  } else if (i > 1) {
    # After the first iteration
    # Get one weighted cluster-specific COV per endo LV
    COV_lv <- vector("list", length = length(endog))
    for (lv in 1:length(endog)) {
      COV_lv[[lv]] <- vector("list", length = nclus)
    }

    for (lv in 1:length(endog)) {
      for (k in 1:nclus) {
        # create 'averaged' sample cov for this cluster
        this_nobs <- post_lv[[lv]][, k] * N_gs
        this_w <- this_nobs / sum(this_nobs)
        tmp <- lapply(seq_along(cov_eta), function(g) {
          cov_eta[[g]] * this_w[g]
        })
        COV_lv[[lv]][[k]] <- Reduce("+", tmp)
      }
    }
  }

  # Call lavaan to estimate the structural parameters
  # the 'groups' are the clusters
  # Note: this makes all resulting parameters to be cluster-specific (it is reconstructed later)
  if (i == 1) {
    # Do this when endo_group_specific is False OR when it is True and we are in the first iteration
    # For the first iteration, perform the full structural model estimation
    s2out <- lavaan::lavaan(
      slotOptions       = fake@Options,
      slotParTable      = fake@ParTable,
      sample.cov        = COV,
      sample.nobs       = rep(sum(N_gs), nclus)
      # slotModel       = slotModel,
      # slotData        = fake@Data,
      # slotSampleStats = fake@SampleStats
    )
  } else if (i > 1) {
    # After the first iteration
    # Run structural estimation once per endo LV
    s2out <- vector(mode = "list", length = length(endog))
    for (lv in 1:length(endog)) {
      s2out[[lv]] <- lavaan::lavaan(
        slotOptions       = fake_lv[[lv]]@Options,
        slotParTable      = fake_lv[[lv]]@ParTable,
        sample.cov        = COV_lv[[lv]],
        sample.nobs       = rep(sum(N_gs), nclus)
        # slotModel       = slotModel,
        # slotData        = fake@Data,
        # slotSampleStats = fake@SampleStats
      )
    }
  }

  # Parameter extraction
  if(i == 1){
    # Extract
    if (nclus == 1) {
      EST_s2  <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
      beta_ks <- EST_s2[["beta"]]
      psi_ks  <- EST_s2[["psi"]]
    } else if (nclus != 1) {
      EST_s2  <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
      beta_ks <- lapply(EST_s2, "[[", "beta") # Does not work with only one cluster
      psi_ks  <- lapply(EST_s2, "[[", "psi")
    }

    # Re order for correct comparisons
    if (nclus == 1) {
      beta_ks <- reorder_long(beta_ks, exog = exog, endog = endog)
      psi_ks <- reorder_long(psi_ks, exog = exog, endog = endog)
    } else if (nclus != 1) {
      beta_ks <- lapply(1:nclus, function(x) {
        reorder_long(beta_ks[[x]], exog = exog, endog = endog)
      }) # Does not work with only one cluster
      psi_ks <- lapply(1:nclus, function(x) {
        reorder_long(psi_ks[[x]], exog = exog, endog = endog)
      })
    }

  } else if(i > 1){
    # Extract the beta matrices per model (one per endo LV)
    # Initialize lists to store the parameters
    EST_s2_lv  <- vector(mode = "list", length = length(endog))
    beta_ks_lv <- vector(mode = "list", length = length(endog))
    psi_ks_lv  <- vector(mode = "list", length = length(endog))
    for (lv in 1:length(endog)) {
      if (nclus == 1) {
        EST_s2_lv[[lv]]  <- lavaan::lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
        beta_ks_lv[[lv]] <- EST_s2_lv[[lv]][["beta"]]
        psi_ks_lv[[lv]]  <- EST_s2_lv[[lv]][["psi"]]
      } else if (nclus != 1) {
        EST_s2_lv[[lv]]  <- lavaan::lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
        beta_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "beta") # Does not work with only one cluster
        psi_ks_lv[[lv]]  <- lapply(EST_s2_lv[[lv]], "[[", "psi")
      }
    }

    # Combine all the beta matrices into just one per cluster

    # Start with an empty beta
    # beta <- matrix(data = 0, nrow = length(lat_var), ncol = length(lat_var))
    # colnames(beta) <- rownames(beta) <- lat_var
    # beta_ks <- lapply(X = seq_along(beta_ks), FUN = function(k){beta_ks[[k]] * 0})

    # beta_ks will contain the regression parameters per cluster
    # beta_ks_lv contains regressions per cluster AND per model of each endo latent variables

    for (k in 1:nclus) {
      for (lv in 1:length(endog)) {
        # Select current endogenous latent variable
        this_lv <- endog[lv]
        col.idx <- if (nclus == 1) {
          colnames(beta_ks_lv[[lv]])
        } else {
          colnames(beta_ks_lv[[lv]][[k]])
        }

        # Extract the regression coefficients of each endogenous latent variables
        if (nclus == 1) {
          beta_ks[this_lv, col.idx] <- beta_ks_lv[[lv]][this_lv, col.idx]
        } else if (nclus != 1) {
          beta_ks[[k]][this_lv, col.idx] <- beta_ks_lv[[lv]][[k]][this_lv, col.idx] # Does not work with only 1 cluster
        }
      }
    }

    # Re-order the beta matrices to make sure we are comparing the correct matrices
    if (nclus == 1) {
      beta_ks <- reorder_long(beta_ks, exog = exog, endog = endog)
    } else if (nclus != 1) {
      beta_ks <- lapply(1:nclus, function(x) {
        reorder_long(beta_ks[[x]], exog = exog, endog = endog)
      }) # Does not work with only one cluster
    }
  }

  # Compute loglikelihood for all group/cluster combinations
  # Initialize matrices to store loglikelihoods
  loglik_gks  <- matrix(data = 0, nrow = (ngroups * ntimes), ncol = nclus)
  gk <- 0

  # Prepare Sigma
  # Initialize the object for estimating Sigma
  Sigma <- matrix(data = list(NA), nrow = (ngroups * ntimes), ncol = nclus)
  I <- diag(length(lat_var)) # Identity matrix based on number of latent variables. Used later

  # Create dummy psi_ks for future inputting
  if(i > 1){
    if(nclus > 1) {
      psi_ks <- psi_gks[, 1]
    } else {
      psi_ks <- psi_gks[, 1][[1]]
    }
  }

  # Previous matrices were only cluster-specific. We have to reconstruct the group-specific matrices (psi and sigma)
  for (k in 1:nclus) {
    ## Save the cluster-specific psi and beta
    ## ifelse() in case of only 1 cluster
    ifelse(test = (nclus == 1), yes = (psi <- psi_ks), no = (psi <- psi_ks[[k]]))
    ifelse(test = (nclus == 1), yes = (beta <- beta_ks), no = (beta <- beta_ks[[k]]))
    for (g in 1:(ngroups * ntimes)) {
      # Reconstruct psi and sigma so they are group- and cluster-specific again.
      # Replace the group-specific part of psi
      # Exogenous (co)variance is always group-specific
      psi[exog, exog] <- cov_eta[[g]][exog, exog] # Replace the group-specific part

      # If the user required group-specific endogenous covariances (endo_group_specific = T), do:
      # Take into account the effect of the cluster-specific regressions
      # cov_eta[[g]] = solve(I - beta) %*% psi %*% t(solve(I - beta))
      # If we solve for psi, then:
      solved_psi <- ((I - beta) %*% cov_eta[[g]] %*% t((I - beta))) # Extract group-specific endog cov

      # Replace endog 1
      g_endog1_cov <- solved_psi[endog1, endog1]
      if (length(endog1) > 1) { # Remove cov between endog 1 variables
        g_endog1_cov[row(g_endog1_cov) != col(g_endog1_cov)] <- 0
      }
      psi[endog1, endog1] <- g_endog1_cov

      # Replace endog 2
      g_endog2_cov <- solved_psi[endog2, endog2] # Extract group-specific endog cov
      psi[endog2, endog2] <- g_endog2_cov

      # If required by the user, set to 0 the covariance between endog 2 factors
      endogenous_cov <- F
      if (length(endog2) > 1 & isFALSE(endogenous_cov)) {
        psi[endog2, endog2][row(psi[endog2, endog2]) != col(psi[endog2, endog2])] <- 0
      }

      # Store for future weighting
      psi_gks[[g, k]] <- psi

      # Get log-likelihood by comparing factor covariance matrix of step 1 (cov_eta) and step 2 (Sigma)
      # Estimate Sigma (factor covariance matrix of step 2)
      Sigma[[g, k]] <- solve(I - beta) %*% psi %*% t(solve(I - beta))

      # Estimate the loglikelihood
      loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
        sample.mean = rep(0, nrow(cov_eta[[1]])),
        sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
        sample.cov  = cov_eta[[g]], # Factor covariance matrix from step 1
        Mu          = rep(0, nrow(cov_eta[[1]])),
        Sigma       = Sigma[[g, k]] # Factor covariance matrix from step 2
      )

      loglik_gks[g, k] <- loglik_gk

    }
  }

  # ----------------------------------------------------------------------------------
  # GET FINAL LOGLIKELIHOOD OBJECT PER TIME POINT
  # ----------------------------------------------------------------------------------
  # The resulting loglikelihood matrices MUST be separated per time points as it will be later used to compute the transitions
  # Initialize empty list
  log_liks <- vector(mode = "list", length = ntimes)
  # browser()
  for(t in 1:ntimes){
    start_idx <- (((t - 1) * ngroups) + 1)
    end_idx   <- t * ngroups
    log_liks[[t]] <- loglik_gks[start_idx:end_idx, , drop = FALSE]
  }

  # Return relevant objects
  return(list(log_liks = log_liks,
              s2out    = s2out,
              psi_gks  = psi_gks,
              beta_ks  = beta_ks))
}
