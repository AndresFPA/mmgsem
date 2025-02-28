#' Standard Errors (SE) wrapper for MMGSEM
#'
#' @description
#' A wrapper to compute the SE of a fitted object by MMG-SEM. It will use the Hessian matrix to compute the SE of both the measurement and the structural parameters.
#'
#' @usage se(object)
#'
#' @param object A resulting fitted object from the MMGSEM function.
#'
#' @return SE: List with the SE for all parameters. The SE are inside matrices mimicking the matrices of the parameters themselves. It also contains the corrected SE.
#' @return HESS: Hessian matrix containing all the second derivatives of the parameters (equivalent to the SE, but in other presentation).
#'
#' @export
se <- function(object){
  # Prepare some objects
  # Remove unnecessary last column of the posteriors object
  object$posteriors <- object$posteriors[, 1:(ncol(object$posteriors) - 1)]
  ngroups <- nrow(object$posteriors)
  nclus   <- ncol(object$posteriors)

  # Calculate the Standard Errors of the structural parameters
  # To calculate the SE, we use the Hessian (second derivative) of the structural parameters.
  # For the derivative, we use a custom function which requires a parameter vector and an obj. function
  # (1) Get all necesseary objects for the objective function
  # Get the parameter vector
  ## Get the parameter vector ##
  # Get the measurement model parameters in a vector
  # Lambda
  # Deal with measurement blocks
  MM     <- object$MM
  partbl <- c()
  if (is.list(MM)){ # Extract the parameter table per block and bind it together for a complete table
    for(m in 1:length(MM)){
      tmp_partbl <- parTable(MM[[m]])
      partbl     <- rbind(partbl, tmp_partbl)
    }
  } else {
    partbl <- parTable(MM)
  }
  partbl$fullpar <- paste0(partbl$lhs, partbl$op, partbl$rhs)
  constrained    <- partbl$fullpar[which(partbl$op == "=~" & is.na(partbl$ustart) & partbl$group == 1 & partbl$label == partbl$plabel)]
  unconstrained  <- partbl$fullpar[which(partbl$op == "=~" & is.na(partbl$ustart) & partbl$group == 1 & partbl$label == "")]

  lambda_vec <- c()
  lambda_nam <- c()

  # Start with the lambdas
  for (g in 1:ngroups) {
    non_zer.idx <- which(unlist(object$param$lambda[[g]]) != 0 & unlist(object$param$lambda[[g]]) != 1)
    lambda_vec <- c(lambda_vec, unlist(object$param$lambda[[g]])[non_zer.idx])
    lambda_nam <- c(lambda_nam, as.vector(outer(X = rownames(object$param$lambda[[g]]),
                                                Y = colnames(object$param$lambda[[g]]),
                                                function(x, y) paste0(y, "=~", x, ".g", g)))[non_zer.idx])

  }

  lambda_vec <- setNames(object = lambda_vec, nm = lambda_nam)
  cons_vec <- lambda_vec[lambda_nam %in% paste0(constrained, ".g1")] # only get the first group (given they are constrained)
  if(length(unconstrained) == 0){unco_vec <- NULL} else {unco_vec <- lambda_vec[grepl(unconstrained, x = names(lambda_vec))]} # Get unconstrained vector only when there are actually unconstrained parameters
  lambda_vec <- c(cons_vec, unco_vec)

  # Theta
  theta_vec <- c()
  theta_nam <- c()

  # Start with the lambdas
  for (g in 1:ngroups) {
    non_zer.idx <- which(unlist(object$param$theta[[g]]) != 0)
    theta_vec <- c(theta_vec, unlist(object$param$theta[[g]])[non_zer.idx])
    theta_nam <- c(theta_nam, as.vector(outer(X = rownames(object$param$theta[[g]]),
                                              Y = colnames(object$param$theta[[g]]),
                                              function(x, y) paste0(y, "~~", x, ".g", g)))[non_zer.idx])

  }

  theta_vec <- setNames(object = theta_vec, nm = theta_nam)

  # Regressions (Cluster-specific!)
  beta_vec <- c() # Empty vector for the betas
  beta_nam <- c() # Empty vector for the names of such betas (not truly necessary, but makes everything easier to understand)

  # The loop identifies the non-zero parameters (i.e., the free parameters) in the beta matrices
  # It also saves the number of the parameter using the col and row names of the beta matrices (e.g., F3 ~ F1)
  for(k in 1:nclus){
    non_zer.idx <- which(unlist(object$param$beta_ks[[k]]) != 0)
    beta_vec <- c(beta_vec, unlist(object$param$beta_ks[[k]])[non_zer.idx])
    beta_nam <- c(beta_nam, as.vector(outer(X = rownames(object$param$beta_ks[[k]]),
                                            Y = colnames(object$param$beta_ks[[k]]),
                                            function(x, y) paste0(x, "~", y, ".k", k)))[non_zer.idx])
  }

  beta_vec <- setNames(object = beta_vec, nm = beta_nam) # Name the vector (the .1 and .2 represent the cluster)

  # Factor covariances (Group-cluster specific!)
  cov_vec <- c() # Empty vector for the covariances
  cov_nam <- c() # Empty vector for the names of such covs

  # The loop is similar the one of the betas.
  # Note there is an extra index (i.e., unique.idx). Given that the cov matrix is symmetric, some parameters are repeated.
  # However, even if repeated, it is only one parameter. The unique.idx removes the repeated parameter.
  for(k in 1:nclus){
    for(g in 1:ngroups){
      object$param$psi_gks[[g, k]] <- round(object$param$psi_gks[[g, k]], 7) # Avoid floting point difference
      non_zer.idx <- which(unlist(object$param$psi_gks[[g, k]]) != 0)
      unique.idx  <- !duplicated(unlist(object$param$psi_gks[[g, k]])[non_zer.idx])
      cov_vec <- c(cov_vec, unlist(object$param$psi_gks[[g, k]])[non_zer.idx][unique.idx])
      cov_nam <- c(cov_nam, as.vector(outer(X = rownames(object$param$psi_gks[[g, k]]),
                                            Y = rownames(object$param$psi_gks[[g, k]]),
                                            function(x, y) paste0(x, "~~", y, ".g", g, ".k", k)))[non_zer.idx][unique.idx])
    }
  }

  cov_vec <- setNames(object = cov_vec, nm = cov_nam) # Name the vector (the first number after the . is the group, the second number is the cluster)

  # Prepare final objects for the hessian function
  # How many lambdas and thetas?
  n_lambda <- length(cons_vec) + (length(unco_vec)/g)
  n_theta  <- length(theta_vec)/g

  # Indices of measurement parameters
  idx.cons <- which(object$param$lambda[[1]] %in% cons_vec)
  idx.unco <- which(object$param$lambda[[1]] %in% unco_vec)
  idx.theta  <- which(object$param$theta[[1]] %in% theta_vec[1:n_theta])

  # How many regressions and covariances?
  n_reg <- length(beta_vec)/nclus # nclus
  n_cov <- length(cov_vec)/(nclus*ngroups) # nclus*ngroups

  # Get indices needed for the objective function (used to input the value of the vector in corresponding place of the matrix)
  idx.beta <- which(object$param$beta_ks[[1]] %in% beta_vec[1:n_reg]) # indices of the free parameters in the beta matrix
  idx.psi  <- which(object$param$psi_gks[[1, 1]] %in% cov_vec) # indices of the free parameters in the psi matrix
  idx.psi_vec <- match(object$param$psi_gks[[1, 1]][idx.psi], cov_vec[1:n_cov]) # To repeat the free parameters that needs to be repeated (i.e, covariances)

  # I have group-specific covariances and group-cluster specific covariances.
  # Which ones are group-specific?
  g.covs  <- unique(cov_vec[duplicated(cov_vec)]) # Group-specific (why duplicated? I extract all possible group-cluster combination, but the group-specific ones will be repeated on the second cluster)
  g.covs  <- setNames(object = g.covs, nm = names(cov_vec[cov_vec %in% g.covs])[1:length(g.covs)])
  gk.covs <- cov_vec[!c(cov_vec %in% g.covs)] # Group-cluster specific

  # (2) Create the objective function of the step 2 parameters
  # The objective function simply takes the parameter vector and fills in the corresponding matrix to get the LL

  # (3) Create a function to compute the Hessian numerically, and use it.
  x <- c(cons_vec, unco_vec, theta_vec, beta_vec, g.covs, gk.covs)
  vec_ind <- list(which(x %in% cons_vec),
                  which(x %in% unco_vec),
                  which(x %in% theta_vec),
                  which(x %in% beta_vec),
                  which(x %in% g.covs),
                  which(x %in% gk.covs)
  )

  HESS <- compute_hessian(f         = obj.S2,
                          x         = x,
                          d         = 1e-02,
                          lambda_mat = object$param$lambda,
                          theta_mat  = object$param$theta,
                          beta_mat   = object$param$beta_ks,
                          psi_mat    = object$param$psi_gks,
                          pi_ks      = colMeans(object$posteriors),
                          S          = object$sample.stats$S,
                          nclus      = k, ngroups = g, n_cov_exo = object$sample.stats$n_cov_exo,
                          N_gs       = object$N_gs,
                          idx.beta   = idx.beta,
                          idx.psi    = idx.psi,
                          idx.psi_vec = idx.psi_vec,
                          idx.cons   = idx.cons,
                          idx.unco   = idx.unco,
                          idx.theta  = idx.theta,
                          vec_ind    = vec_ind)

  # (4) Organize the SE for each parameter
  vector_SE <- setNames(diag(sqrt(ginv(-HESS, tol = 1e-06))), colnames(HESS)) # The SE comes from the inverse of the negative hessian

  SE.S2 <- function(x, lambda_mat, theta_mat,
                    beta_mat, psi_mat,
                    pi_ks, S,
                    nclus, ngroups, N_gs,
                    idx.cons, idx.unco, idx.theta,
                    idx.beta, idx.psi, idx.psi_vec,
                    n_cov_exo, vec_ind){

    con_lam <- x[vec_ind[[1]]]
    unc_lam <- x[vec_ind[[2]]]
    thetas  <- x[vec_ind[[3]]]
    betas   <- x[vec_ind[[4]]]
    g.covs  <- x[vec_ind[[5]]]
    gk.covs <- x[vec_ind[[6]]]

    n_g.psi <- length(g.covs)
    n_gk.psi <- length(gk.covs)
    n_reg <- length(betas)

    # Relocating lambda
    g.lam.idx <- length(unc_lam)/ngroups

    for(g in 1:ngroups){
      unco_vec <- unc_lam[(((g - 1)*g.lam.idx) + 1):(g.lam.idx*g)]
      lambda_mat[[g]][idx.cons] <- con_lam
      lambda_mat[[g]][idx.unco] <- unco_vec
    }

    # Relocating theta
    g.theta.idx <- length(thetas)/ngroups
    # browser()
    for(g in 1:ngroups){
      theta_vec <- thetas[(((g - 1)*g.theta.idx) + 1):(g.theta.idx*g)]
      theta_mat[[g]][idx.theta] <- theta_vec
    }

    # Relocating beta
    clus_idx <- n_reg/nclus
    for(k in 1:nclus){
      beta_vec <- betas[(((k - 1)*clus_idx) + 1):(clus_idx*k)]
      beta_mat[[k]][idx.beta] <- beta_vec
    }

    # Relocating psi
    psi_gks <- psi_mat
    g.clus_idx <- n_g.psi/(ngroups)
    gk.clus_idx <- n_gk.psi/(nclus*ngroups)
    gk <- 0
    for(k in 1:nclus){
      for(g in 1:ngroups){
        gk <- gk + 1
        cov_vec <- c(g.covs[(((g - 1)*g.clus_idx) + 1):(g.clus_idx*g)], gk.covs[(((gk - 1)*gk.clus_idx) + 1):(gk.clus_idx*gk)])
        psi_mat[[g, k]][idx.psi] <- cov_vec[idx.psi_vec]
      }
    }

    return(list(
      lambda_SE = lambda_mat,
      theta_SE  = theta_mat,
      betas_SE  = beta_mat,
      psi_SE    = psi_mat,
      SE_vector = x
    ))
  }

  SE <- SE.S2(x          = vector_SE,
              lambda_mat = object$param$lambda,
              theta_mat  = object$param$theta,
              beta_mat   = object$param$beta_ks,
              psi_mat    = object$param$psi_gks,
              pi_ks      = colMeans(object$posteriors),
              S          = object$sample.stats$S,
              nclus      = k, ngroups = g, n_cov_exo = object$sample.stats$n_cov_exo,
              N_gs       = object$N_gs,
              idx.beta   = idx.beta,
              idx.psi    = idx.psi,
              idx.psi_vec = idx.psi_vec,
              idx.cons   = idx.cons,
              idx.unco   = idx.unco,
              idx.theta  = idx.theta,
              vec_ind    = vec_ind)

  # SE Correction ---
  Sigma_1 <- vcov(object$MM) # Extract Sigma step 1 directly from lavaan
  Sigma_1 <- Sigma_1[, !duplicated(colnames(Sigma_1))] # Lavaan returns the constrained parameters duplicated. Remove them
  colnames(Sigma_1)[grepl(pattern = ".p", x = colnames(Sigma_1))] <- names(cons_vec) # Rename the constrained parameter names to the actual parameter
  colnames(Sigma_1)[grep(pattern = ".g", x = colnames(Sigma_1), invert = T)] <- paste0(colnames(Sigma_1)[grep(pattern = ".g", x = colnames(Sigma_1), invert = T)], ".g1") # Add g1 to the first group

  rownames(Sigma_1)[grepl(pattern = ".p", x = rownames(Sigma_1))] <- names(cons_vec) # Rename the constrained parameter names to the actual parameter
  rownames(Sigma_1)[grep(pattern = ".g", x = rownames(Sigma_1), invert = T)] <- paste0(rownames(Sigma_1)[grep(pattern = ".g", x = rownames(Sigma_1), invert = T)], ".g1") # Add g1 to the first group

  Sigma_1 <- Sigma_1[c(names(cons_vec), names(unco_vec), names(theta_vec)),
                     c(names(cons_vec), names(unco_vec), names(theta_vec))] # Re-order in the same way as the one used in the joint hessian

  # Get indices
  step1.idx <- which(colnames(HESS) %in% c(names(lambda_vec), names(theta_vec)))
  step2.idx <- which(colnames(HESS) %in% c(names(beta_vec), names(cov_vec)))

  # Get matrix divided
  I_22 <- -HESS[step2.idx, step2.idx]
  I_21 <- -HESS[step2.idx, step1.idx]

  I_22.inv <- MASS::ginv(I_22, tol = 1e-06)
  Sigma_1  <- Sigma_1[step1.idx, step1.idx]
  # browser()

  # APPLY THE CORRECTION (note: for now, it is not changing too much the results)
  Sigma_2_corrected <- diag(sqrt(I_22.inv + I_22.inv %*% I_21 %*% Sigma_1 %*% t(I_21) %*% I_22.inv))

  # Add the corrected values to the final results
  SE$corrected <- setNames(object = Sigma_2_corrected, nm = colnames(HESS[step2.idx, step2.idx]))

  return(list(SE = SE, HESS = HESS, est_vec = x, se_vector = vector_SE))
}

# Objective function for the Hessian ---------------------------------------------------------------

obj.S2 <- function(x, lambda_mat, theta_mat,
                   beta_mat, psi_mat,
                   pi_ks, S,
                   nclus, ngroups, N_gs,
                   idx.cons, idx.unco, idx.theta,
                   idx.beta, idx.psi, idx.psi_vec,
                   n_cov_exo, vec_ind){

  con_lam <- x[vec_ind[[1]]]
  unc_lam <- x[vec_ind[[2]]]
  thetas  <- x[vec_ind[[3]]]
  betas   <- x[vec_ind[[4]]]
  g.covs  <- x[vec_ind[[5]]]
  gk.covs <- x[vec_ind[[6]]]

  n_g.psi <- length(g.covs)
  n_gk.psi <- length(gk.covs)
  n_reg <- length(betas)

  # Relocating lambda
  g.lam.idx <- length(unc_lam)/ngroups

  for(g in 1:ngroups){
    unco_vec <- unc_lam[(((g - 1)*g.lam.idx) + 1):(g.lam.idx*g)]
    lambda_mat[[g]][idx.cons] <- con_lam
    lambda_mat[[g]][idx.unco] <- unco_vec
  }

  # Relocating theta
  g.theta.idx <- length(thetas)/ngroups
  # browser()
  for(g in 1:ngroups){
    theta_vec <- thetas[(((g - 1)*g.theta.idx) + 1):(g.theta.idx*g)]
    theta_mat[[g]][idx.theta] <- theta_vec
  }

  # Relocating beta
  clus_idx <- n_reg/nclus
  for(k in 1:nclus){
    beta_vec <- betas[(((k - 1)*clus_idx) + 1):(clus_idx*k)]
    beta_mat[[k]][idx.beta] <- beta_vec
  }

  # Relocating psi
  psi_gks <- psi_mat
  g.clus_idx <- n_g.psi/(ngroups)
  gk.clus_idx <- n_gk.psi/(nclus*ngroups)
  gk <- 0
  for(k in 1:nclus){
    for(g in 1:ngroups){
      gk <- gk + 1
      cov_vec <- c(g.covs[(((g - 1)*g.clus_idx) + 1):(g.clus_idx*g)], gk.covs[(((gk - 1)*gk.clus_idx) + 1):(gk.clus_idx*gk)])
      psi_mat[[g, k]][idx.psi] <- cov_vec[idx.psi_vec]
    }
  }


  # compute loglikelihood for all group/cluster combinations
  # Initialize matrices to store loglikelihoods
  loglik_gks <- matrix(data = 0, nrow = ngroups, ncol = nclus)
  loglik_gksw <- matrix(data = 0, nrow = ngroups, ncol = nclus)

  Sigma <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
  I <- diag(nrow(psi_mat[[1, 1]]))

  for(k in 1:nclus){
    for(g in 1:ngroups){
      # Estimate Sigma (factor covariance matrix of step 2)
      Sigma[[g, k]] <- lambda_mat[[g]] %*% solve(I - beta_mat[[k]]) %*% psi_mat[[g, k]] %*% t(solve(I - beta_mat[[k]])) %*% t(lambda_mat[[g]]) + theta_mat[[g]]

      S_biased <- S[[g]] * (N_gs[[g]] - 1) / N_gs[[g]]

      # Estimate the loglikelihood
      loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
        sample.mean = rep(0, nrow(theta_mat[[1]])),
        sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
        sample.cov  = S_biased, # Factor covariance matrix from step 1
        Mu          = rep(0, nrow(theta_mat[[1]])),
        Sigma       = Sigma[[g, k]] # Factor covariance matrix from step 2
      )

      loglik_gks[g, k] <- loglik_gk
      loglik_gksw[g, k] <- log(pi_ks[k]) + loglik_gk # weighted loglik
    }
  }

  # Get total loglikelihood
  # First, deal with arithmetic underflow by subtracting the maximum value per group
  max_gs <- apply(loglik_gksw, 1, max) # Get max value per row
  minus_max <- sweep(x = loglik_gksw, MARGIN = 1, STATS = max_gs, FUN = "-") # Subtract the max per row
  exp_loglik <- exp(minus_max) # Exp before summing for total loglikelihood
  loglik_gsw <- log(apply(exp_loglik, 1, sum)) # Sum exp_loglik per row and then take the log again
  LL <- sum((loglik_gsw + max_gs)) # Add the maximum again and then sum them all for total loglikelihood

  return(LL)
}


# Hessian function to compute the second derivative ------------------------------------------------

compute_hessian <- function(f, x, d = 1e-5,
                            lambda_mat, theta_mat,
                            beta_mat, psi_mat,
                            pi_ks, S,
                            nclus, ngroups, N_gs,
                            idx.cons, idx.unco, idx.theta,
                            idx.beta, idx.psi, idx.psi_vec,
                            n_cov_exo, vec_ind) {
  # browser()
  n  <- length(x)
  H  <- matrix(0, n, n)
  colnames(H) <- rownames(H) <- names(x)
  f1 <- matrix(0, n)

  f0 <- f(x = x,
          lambda_mat = lambda_mat,
          theta_mat  = theta_mat,
          beta_mat = beta_mat,
          psi_mat = psi_mat,
          pi_ks = pi_ks,
          S = S,
          nclus = nclus,
          ngroups = ngroups,
          N_gs = N_gs,
          idx.beta = idx.beta,
          idx.cons = idx.cons,
          idx.unco = idx.unco,
          idx.theta = idx.theta,
          idx.psi = idx.psi,
          idx.psi_vec = idx.psi_vec,
          n_cov_exo = n_cov_exo,
          vec_ind = vec_ind)

  for (i in 1:n) {
    x_up    <- x
    x_up[i] <- x_up[i] + d

    f1[i]   <- f(x = x_up,
                 lambda_mat = lambda_mat,
                 theta_mat  = theta_mat,
                 beta_mat = beta_mat,
                 psi_mat = psi_mat,
                 pi_ks = pi_ks,
                 S = S,
                 nclus = nclus,
                 ngroups = ngroups,
                 N_gs = N_gs,
                 idx.beta = idx.beta,
                 idx.cons = idx.cons,
                 idx.unco = idx.unco,
                 idx.theta = idx.theta,
                 idx.psi = idx.psi,
                 idx.psi_vec = idx.psi_vec,
                 n_cov_exo = n_cov_exo,
                 vec_ind = vec_ind)

    for (j in 1:i) {
      x_up[j] <- x_up[j] + d

      f2      <- f(x = x_up,
                   lambda_mat = lambda_mat,
                   theta_mat  = theta_mat,
                   beta_mat = beta_mat,
                   psi_mat = psi_mat,
                   pi_ks = pi_ks,
                   S = S,
                   nclus = nclus,
                   ngroups = ngroups,
                   N_gs = N_gs,
                   idx.beta = idx.beta,
                   idx.cons = idx.cons,
                   idx.unco = idx.unco,
                   idx.theta = idx.theta,
                   idx.psi = idx.psi,
                   idx.psi_vec = idx.psi_vec,
                   n_cov_exo = n_cov_exo,
                   vec_ind = vec_ind)

      H[i, j] <- (f2 - f1[i] - f1[j] + f0) / (d^2)
      H[j, i] <- H[i, j]

      x_up[j] <- x_up[j] - d
    }
  }

  return(H)
}
