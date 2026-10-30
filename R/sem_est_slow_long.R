sem_est_slow_long <- function(S2, post, cov_eta, nstates, ntimes, ngroups, N_gs){
  # Swap back the dimensions of the posterior matrices for easier manipulation
  post <- aperm(a = post, perm = c(1,3,2))

  # Necessary objects
  nclus <- nstates
  gro_clu <- nclus*ngroups*ntimes

  # ----------------------------------------------------------------------------
  # CONSTRAINT DEFINITION
  # ----------------------------------------------------------------------------

  # Create a dummy Step 2 parameter table.
  fake_cov        <- rep(cov_eta, nclus) # Duplicate factor's covariances to match the number of clusters
  names(fake_cov) <- paste("group", seq_len(gro_clu))
  fake_model      <- lavaan::parTable(lavaan::sem(model = S2,
                                                  sample.cov = fake_cov,
                                                  sample.nobs = rep(N_gs, nclus),
                                                  do.fit = FALSE,
                                                  meanstructure = F,
                                                  h1 = FALSE,
                                                  check.post = FALSE,
                                                  loglik = FALSE,
                                                  sample.cov.rescale = FALSE,
                                                  fixed.x = TRUE
  ))

  # Remove unnecessary columns from fake_global
  fake_model$se      <- NULL
  fake_model$est     <- NULL
  fake_model$start   <- NULL

  # Get the labels of the endogenous 1 and 2 factors
  lat_var <- lavNames(lavaanify(S2, auto = TRUE), "ov")
  endog1  <- lat_var[(lat_var %in% fake_model$rhs[which(fake_model$op == "~")]) &
                       (lat_var %in% fake_model$lhs[which(fake_model$op == "~")])]
  endog2  <- lat_var[!c(lat_var %in% fake_model$rhs[which(fake_model$op == "~")]) &
                       (lat_var %in% fake_model$lhs[which(fake_model$op == "~")])]
  endog   <- c(endog1, endog2)
  exog    <- lat_var[!c(lat_var %in% endog)]

  # Add constraints, "normal" approach
  # Start the process to add cluster constraints
  # Create a constraint entry in the lavaan format
  constraints_row <- data.frame(
    id = "", lhs = "", op = "==", rhs = "",
    user = 2, block = 0, group = 0, free = 0,
    ustart = NA, exo = 0, label = "", plabel = "",
    cluster = NA
  )


  # constraints object refer to regression parameters constraints.
  # cons_exo_cov object refer to covariance parameters of exogenous variables (should be group-specific)

  # Identify regression parameters
  constraints <- fake_model$plabel[which(fake_model$op == "~")] # Get the regression parameters
  n_reg       <- length(fake_model$plabel[which(fake_model$op == "~" & fake_model$group == 1)]) # Number of reg PER GROUP


  # Identify independent latent variables in the model
  # Labels of the variance parameters of variables in exo
  cons_exo <- fake_model$plabel[which(fake_model$op == "~~" & fake_model$lhs %in% exog)]
  # Number of variance parameters that involve variables in exo PER GROUP
  n_exo    <- length(fake_model$plabel[which(fake_model$op == "~~" &
                                               fake_model$lhs %in% exog & fake_model$group == 1)])

  # Create matrices with the necessary constraints entries
  constraints_matrix <- constraints_row[rep(
    x = 1:nrow(constraints_row),
    times = (length(constraints))
  ), ]

  cons_exo_matrix <- constraints_row[rep(
    x = 1:nrow(constraints_row),
    times = (length(cons_exo))
  ), ]

  rownames(constraints_matrix) <- NULL
  rownames(cons_exo_matrix)    <- NULL

  # Get a cluster label for all groups (all combinations group*cluster)
  clus_label  <- rep(x = 1:nclus, each = (ngroups*ntimes))
  group_label <- rep(x = 1:(ngroups*ntimes), times = nclus)

  # Add cluster labels to the parameter table (not necessary, just for me)
  for (j in 1:length(clus_label)) {
    fake_model$cluster[fake_model$group == j] <- clus_label[j]
  }

  # Repeat each cluster label depending on the number of parameters per group.
  # i.e. Label each parameter per cluster
  reg_labels <- rep(clus_label, each = n_reg) # regressions
  exo_labels <- rep(group_label, each = n_exo) # variance of endog1

  # Add constraints per cluster (i.e. regression parameters are equal within cluster)
  for (k in 1:nclus) {
    # Regression constraints
    cluster_par <- constraints[reg_labels == k] # Identify regression parameters of cluster k
    constraints_matrix[(reg_labels == k), "lhs"] <- cluster_par[1:n_reg] # On the left hand side insert parameters of ONE group
    constraints_matrix[(reg_labels == k), "rhs"] <- cluster_par          # On the right hand side insert parameters of all groups
  }

  # Add constraints per group (i.e., exo covs are equal within a group - NOT a group-cluster parameter)
  for(g in 1:(ngroups*ntimes)){
    # Variances constraints (exo)
    group_par_exo <- cons_exo[exo_labels == g]
    cons_exo_matrix[(exo_labels == g), "lhs"] <- group_par_exo[1:n_exo]
    cons_exo_matrix[(exo_labels == g), "rhs"] <- group_par_exo
  }

  # Remove redundant constraints (e.g., p1 == p1)
  constraints_total <- rbind(constraints_matrix, cons_exo_matrix)
  redundant <- which(constraints_total$lhs == constraints_total$rhs) # Identify redundant
  constraints_total <- constraints_total[-redundant, ]
  rownames(constraints_total) <- NULL

  # Bind the model table with the constraints
  fake_model <- rbind(fake_model, constraints_total)
  fake_model$free <- 1:nrow(fake_model)

  # ----------------------------------------------------------------------------
  # PARAMETER ESTIMATION
  # ----------------------------------------------------------------------------
  # Create an stacked posterior matrix to estimate the structural parameters ACROSS time points
  # This ensures that the cluster labels remain the same across time (i.e., avoid label switching)
  post_stacked <- c()
  for(t in 1:ntimes){
    post_stacked <- rbind(post_stacked, post[, , t])
  }
  # browser()
  # Obtain weighted group-cluster combination sample size
  N_gks <- post_stacked * N_gs
  N_gks <- c(N_gks)

  # Duplicate the covariance matrices nclus times
  # fake_cov <- rep(cov_eta, nclus)
  # browser()
  # Call lavaan to estimate the structural parameters
  s2out <- sem(
    model = fake_model,    # 'fake' partable with duplicated parameters and so on
    sample.cov = fake_cov, # 'fake' duplicated factors' cov matrix (repeated by the number of clusters)
    sample.nobs = N_gks,   # Sample size per group-cluster combination weighted by the posteriors
    baseline = FALSE, se = "none",
    h1 = FALSE, check.post = FALSE,
    control = list(rel.tol = 1e-09),
    sample.cov.rescale = FALSE,
    fixed.x = FALSE
  ) # , control = list(max.iter = 50))

  # Compute loglikelihood for all group/cluster combinations
  # Initialize matrices to store loglikelihoods
  loglik_gks  <- matrix(data = 0, nrow = (ngroups * ntimes), ncol = nclus)
  gk <- 0

  for (k in 1:nclus) {
    for (g in 1:(ngroups * ntimes)) {
      gk <- gk + 1L
      loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
        sample.mean = s2out@SampleStats@mean[[gk]],
        sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
        sample.cov  = s2out@SampleStats@cov[[gk]],
        Mu          = s2out@SampleStats@mean[[gk]], # s2out@SampleStats@mean[[gk]],
        Sigma       = s2out@implied$cov[[gk]]
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
    log_liks[[t]] <- loglik_gks[start_idx:end_idx, ]
  }

  estimates <- lavaan::lavInspect(s2out, "est", add.class = TRUE,
                                  add.labels = TRUE)
  beta_ks <- vector("list", nclus)
  psi_gks <- matrix(vector("list", ngroups * ntimes * nclus),
                    nrow = ngroups * ntimes, ncol = nclus)
  for (k in seq_len(nclus)) {
    for (g in seq_len(ngroups * ntimes)) {
      estimate <- estimates[[(k - 1L) * ngroups * ntimes + g]]
      if (g == 1L) {
        beta_ks[[k]] <- estimate$beta
      }
      psi_gks[[g, k]] <- estimate$psi
    }
  }
  if (nclus == 1L) {
    beta_ks <- beta_ks[[1]]
  }

  # Return relevant objects
  return(list(log_liks = log_liks, s2out = s2out,
              beta_ks = beta_ks, psi_gks = psi_gks))
}
