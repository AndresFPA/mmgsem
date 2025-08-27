#' Standard Errors (SE) wrapper for MMGSEM
#'
#' @description
#' A wrapper to compute the SE of a fitted object by mmgsem(). It obtains the SE by numerically computing the second derivative of the loglikelihood function.
#'
#' @usage se(object, d = 1e-03, naive = FALSE)
#'
#' @param object A fitted model returned by the mmgsem function.
#' @param d Step size for numerical approximation of the derivative. Default is 1e-03.
#' @param naive Logical. If TRUE, returns uncorrected SEs (faster but less accurate). If FALSE (default), applies the correction for step-wise estimation.
#'
#' @returns SE_vector: List. Includes several vectors that contain the standard errors per type of parameter (beta, psi, etc.)
#'
#' @export
compute_se <- function(object, d = 1e-03, naive = FALSE){

  # ngroups
  ngroups <- length(object$param$lambda)
  nclus   <- (ncol(object$posteriors) - 1)
  MM      <- object$MM

  # Prepare the necessary parameter matrices
  lambda <- object$param$lambda
  theta  <- object$param$theta
  beta   <- object$param$beta_ks

  # Psi requires extra peparation
  psi_gks <- object$param$psi_gks

  # We remove the 'nuisance' parameters. i.e., the parameters from the empty group-cluster combinations
  # Extract relevant psi_gks
  psi_idx <- which(x = mmgsem_fit$modal_post == 1, arr.ind = T) # array indices of the relevant group-cluster combinations

  # Reorder based on group number insted of column number
  correct_idx <- order(psi_idx[,1])
  psi_idx <- psi_idx[correct_idx,]

  psi <- vector("list", ngroups)
  for(g in 1:ngroups){
    psi[[g]] <- psi_gks[[psi_idx[g, 1], psi_idx[g, 2]]]
  }

  #################################
  ### Flatten all the matrices! ###
  #################################
  # First, get the necessary mapping matrices
  map_lambda <- mapping_params(mat_list = lambda)
  map_theta  <- mapping_params(mat_list = theta)
  map_beta   <- mapping_params(mat_list = beta)
  map_psi    <- mapping_params(mat_list = psi)

  # Put them together in a list to use when unflattening
  maps <- list(
    lambda = map_lambda,
    theta  = map_theta,
    beta   = map_beta,
    psi    = map_psi
  )

  # Now, flatten all the matrices
  flat_lambda <- flatten_params(param_list = lambda, idxs = map_lambda, type = "lambda")
  flat_theta  <- flatten_params(param_list = theta,  idxs = map_theta,  type = "theta")
  flat_beta   <- flatten_params(param_list = beta,   idxs = map_beta,   type = "beta")
  flat_psi    <- flatten_params(param_list = psi ,   idxs = map_psi,    type = "psi")

  # Put them together in a list to use when unflattening
  flats <- list(
    lambda = flat_lambda,
    theta  = flat_theta,
    beta   = flat_beta,
    psi    = flat_psi
  )

  ##############################
  ### Time to compute the SE ###
  ##############################

  ##########################
  # STEP 2 STANDARD ERRORS #
  ##########################
  # Get step 2 parameters into one single parameter vector (NAIVE standard errors)
  param_vec_S2 <- c(flat_beta$vec, flat_psi$vec)

  # We need the indices of each type of parameter before we continue (it is to reconstruct the matrices later)
  vec_ind_S2 <- vector("list", 3)
  vec_ind_S2[[1]] <- which(param_vec_S2 %in% flat_beta$vec)
  vec_ind_S2[[2]] <- which(param_vec_S2 %in% flat_psi$vec)
  vec_ind_S2[[3]] <- psi_idx # Psi matrix indices for psi_gks

  # First, let's compute the NAIVE standard errors
  HESS.S2 <- compute_hessian(f       = obj.S2,
                             x       = param_vec_S2,
                             d       = d,
                             psi_mat = psi_gks,
                             pi_ks   = colMeans(object$posteriors[, -1]), # We need to remove the first column as it is the group indicator
                             S       = object$param$cov_eta,
                             nclus   = nclus,
                             ngroups = ngroups,
                             N_gs    = object$N_gs,
                             vec_ind = vec_ind_S2,
                             maps    = maps,
                             flats   = flats)

  if(naive == TRUE){ # If the user only wants the naive standard error, stop here.
    # Vcov matrix of the step 2 parameters
    Sigma_2 <- MASS::ginv(-HESS.S2)
    colnames(Sigma_2) <- rownames(Sigma_2) <- colnames(HESS.S2)

    # vector form
    vector_SE <- sqrt(diag(Sigma_2))
    vector_SE <- setNames(object = vector_SE, nm = colnames(Sigma_2))

    # Organize results in vector form -----------------------------------
    SE_vector <- list(
      betas_se   = vector_SE[vec_ind_S2[[1]]],
      psi_se     = vector_SE[vec_ind_S2[[2]]]
    )

    estimates_vector <- list(
      betas_est   = param_vec_S2[vec_ind_S2[[1]]],
      psi_est     = param_vec_S2[vec_ind_S2[[2]]]
    )

    # Vcov matrix of the regression parameters
    vcov_betas <- Sigma_2[grepl("^[^~]*~[^~]*$", colnames(Sigma_2)), grepl("^[^~]*~[^~]*$", colnames(Sigma_2))]

    return(list(
      SE_vector  = SE_vector,
      HESS       = HESS.S2,
      vcov_betas = vcov_betas,
      estimates_vector = estimates_vector)
    )
  }

  ##########################
  # STEP 1 STANDARD ERRORS #
  ##########################
  # Get step 2 parameters into one single parameter vector
  param_vec_S1 <- c(flat_lambda$vec, flat_theta$vec)

  # We need the indices of each type of parameter before we continue (it is to reconstruct the matrices later)
  vec_ind_S1 <- vector("list", 3)
  vec_ind_S1[[1]] <- which(param_vec_S1 %in% flat_lambda$vec)
  vec_ind_S1[[2]] <- which(param_vec_S1 %in% flat_theta$vec)

  # Now, let's compute the standard errors of step 1
  HESS.S1 <- compute_hessian(f       = obj.S1,
                             x       = param_vec_S1,
                             d       = d,
                             psi_mat = object$param$cov_eta, # Factor's covariance matrix is necessary to compute the LL of Step 1
                             pi_ks   = colMeans(object$posteriors[, -1]), # We need to remove the first column as it is the group indicator
                             S       = object$sample.stats$S,
                             nclus   = nclus,
                             ngroups = ngroups,
                             N_gs    = object$N_gs,
                             vec_ind = vec_ind_S1,
                             maps    = maps,
                             flats   = flats)

  ##########################
  # JOINT STANDARD ERRORTS #
  ##########################
  # Get everything into one single parameter vector
  param_vec <- c(flat_lambda$vec, flat_theta$vec, flat_beta$vec, flat_psi$vec)

  # We need the indices of each type of parameter before we continue (it is to reconstruct the matrices later)
  vec_ind <- vector("list", 5)
  vec_ind[[1]] <- which(param_vec %in% flat_lambda$vec)
  vec_ind[[2]] <- which(param_vec %in% flat_theta$vec)
  vec_ind[[3]] <- which(param_vec %in% flat_beta$vec)
  vec_ind[[4]] <- which(param_vec %in% flat_psi$vec)
  vec_ind[[5]] <- psi_idx # Psi matrix indices for psi_gks

  # Now, compute the cross-derivatives
  HESS <- compute_hessian(f       = obj,
                          x       = param_vec,
                          d       = d,
                          psi_mat = psi_gks,
                          pi_ks   = colMeans(object$posteriors[, -1]), # We need to remove the first column as it is the group indicator
                          S       = object$sample.stats$S,
                          nclus   = nclus,
                          ngroups = ngroups,
                          N_gs    = object$N_gs,
                          vec_ind = vec_ind,
                          maps    = maps,
                          flats   = flats)

  # Organize the SE for each parameter
  # Compute the inverse of the -Hessian to obtain the SE
  HESS_inv  <- MASS::ginv(-HESS, tol = 1e-06)
  vector_SE <- sqrt(diag(HESS_inv))
  vector_SE <- setNames(vector_SE, colnames(HESS)) # This is done to get the cross-derivatives

  ##############################
  ### SE Correction ############
  ##############################
  # First, we need to get vcov of each step independently
  Sigma_1 <- MASS::ginv(-HESS.S1)
  Sigma_2 <- MASS::ginv(-HESS.S2)

  # Get indices for parameters of both steps from the JOINT Hessian
  step1.idx <- which(colnames(HESS) %in% c(names(flat_lambda$vec), names(flat_theta$vec)))
  step2.idx <- which(colnames(HESS) %in% c(names(flat_beta$vec), names(flat_psi$vec)))

  # Get matrix divided
  I_22 <- HESS[step2.idx, step2.idx]
  I_21 <- HESS[step2.idx, step1.idx]
  I_11 <- HESS[step1.idx, step1.idx]

  I_22.inv <- MASS::ginv(I_22)

  # browser()
  # if(is.list(MM)){
  #   Sigma_1 <- I_11.inv <- MASS::ginv(-I_11, tol = 1e-06) # Equivalent(?) to Sigma_11
  # } else { # If we do not have measurement blocks, lavaan already gives me the parameter sampling variance
  #   Sigma_1 <- lavaan::vcov(object$MM) # Extract Sigma step 1 directly from lavaan
  #   Sigma_1 <- Sigma_1[, !duplicated(colnames(Sigma_1))] # Lavaan returns the constrained parameters duplicated. Remove them
  #
  #   #colnames(Sigma_1)[grepl(pattern = ".p", x = colnames(Sigma_1))] <- names(cons_vec) # Rename the constrained parameter names to the actual parameter
  #   #colnames(Sigma_1)[grep(pattern = ".g", x = colnames(Sigma_1), invert = T)] <- paste0(colnames(Sigma_1)[grep(pattern = ".g", x = colnames(Sigma_1), invert = T)], ".g1") # Add g1 to the first group
  #   #rownames(Sigma_1)[grepl(pattern = ".p", x = rownames(Sigma_1))] <- names(cons_vec) # Rename the constrained parameter names to the actual parameter
  #   #rownames(Sigma_1)[grep(pattern = ".g", x = rownames(Sigma_1), invert = T)] <- paste0(rownames(Sigma_1)[grep(pattern = ".g", x = rownames(Sigma_1), invert = T)], ".g1") # Add g1 to the first group
  #
  #   Sigma_1 <- Sigma_1[c(names(cons_vec), names(unco_vec), names(theta_vec)),
  #                      c(names(cons_vec), names(unco_vec), names(theta_vec))] # Re-order in the same way as the one used in the joint hessian
  #
  #   Sigma_1  <- Sigma_1[step1.idx, step1.idx]
  # }

  # APPLY THE CORRECTION (note: for now, it is not changing too much the results)
  Sigma_2_corrected           <- Sigma_2 + I_22.inv %*% I_21 %*% Sigma_1 %*% t(I_21) %*% I_22.inv
  colnames(Sigma_2_corrected) <- rownames(Sigma_2_corrected) <- colnames(HESS[step2.idx, step2.idx])
  vcov_betas                  <- Sigma_2_corrected[grepl("^[^~]*~[^~]*$", colnames(Sigma_2_corrected)), grepl("^[^~]*~[^~]*$", colnames(Sigma_2_corrected))]

  # vector form
  vector_corrected <- diag(sqrt(Sigma_2_corrected))
  vector_corrected <- setNames(object = vector_corrected, nm = colnames(HESS[step2.idx, step2.idx]))
  betas_se_corrected <- vector_corrected[grepl("^[^~]*~[^~]*$", names(vector_corrected))]

  # Organize results in vector form -----------------------------------
  SE_vector <- list(
    lambda_se  = vector_SE[vec_ind[[1]]],
    thetas_se  = vector_SE[vec_ind[[2]]],
    betas_se   = vector_SE[vec_ind[[3]]],
    psi_se     = vector_SE[vec_ind[[4]]],
    betas_se_corrected = betas_se_corrected
  )

  estimates_vector <- list(
    lambda_est  = param_vec[vec_ind[[1]]],
    thetas_est  = param_vec[vec_ind[[2]]],
    betas_est   = param_vec[vec_ind[[3]]],
    psi_est     = param_vec[vec_ind[[4]]]
  )

  return(list(
    SE_vector  = SE_vector,
    HESS       = HESS,
    vcov_betas = vcov_betas,
    estimates_vector = estimates_vector)
  )

}

##### SECONDARY FUNCTIONS ----------------------------------------------------------------------------------------------
#### HESSIAN -----------------------------------------------------------------------------------------------------------
# Hessian function to compute the numerical second derivative of our likelihood function
compute_hessian <- function(f,
                            x,
                            d,
                            psi_mat,
                            pi_ks,
                            S,
                            nclus,
                            ngroups,
                            N_gs,
                            vec_ind,
                            maps,
                            flats) {

  # Prepare the matrices necessary to store the results
  n  <- length(x)
  H  <- matrix(0, n, n)
  colnames(H) <- rownames(H) <- names(x)
  f1 <- matrix(0, n)

  # Loglikelihood result without adding the small number
  f0 <- f(x       = x,
          psi_mat = psi_mat,
          pi_ks   = pi_ks,
          S       = S,
          nclus   = nclus,
          ngroups = ngroups,
          N_gs    = N_gs,
          vec_ind = vec_ind,
          maps    = maps,
          flats   = flats)

  for (i in 1:n) {
    x_up    <- x
    x_up[i] <- x_up[i] + d

    f1[i]   <- f(x       = x_up,
                 psi_mat = psi_mat,
                 pi_ks   = pi_ks,
                 S       = S,
                 nclus   = nclus,
                 ngroups = ngroups,
                 N_gs    = N_gs,
                 vec_ind = vec_ind,
                 maps    = maps,
                 flats   = flats)

    for (j in 1:i) {
      x_up[j] <- x_up[j] + d

      f2      <- f(x       = x_up,
                   psi_mat = psi_mat,
                   pi_ks   = pi_ks,
                   S       = S,
                   nclus   = nclus,
                   ngroups = ngroups,
                   N_gs    = N_gs,
                   vec_ind = vec_ind,
                   maps    = maps,
                   flats   = flats)

      H[i, j] <- (f2 - f1[i] - f1[j] + f0) / (d^2)
      H[j, i] <- H[i, j]

      x_up[j] <- x_up[j] - d
    }
  }

  return(H)
}
#### OBJECTIVE FUNCTION ---------------------------------------------------------------------------------------------------
# The function that computes the general loglikelihood
obj <- function(x,
                psi_mat,
                pi_ks,
                S,
                nclus,
                ngroups,
                N_gs,
                vec_ind,
                maps,
                flats){

  lambda_vec <- x[vec_ind[[1]]]
  theta_vec  <- x[vec_ind[[2]]]
  beta_vec   <- x[vec_ind[[3]]]
  psi_vec    <- x[vec_ind[[4]]]
  psi_idx    <- vec_ind[[5]]

  lambda_mat <- unflatten_params(flat = lambda_vec, meta = flats$lambda$meta, idxs = maps$lambda, type = "lambda")
  theta_mat  <- unflatten_params(flat = theta_vec, meta = flats$theta$meta, idxs = maps$theta, type = "theta")
  beta_mat   <- unflatten_params(flat = beta_vec, meta = flats$beta$meta, idxs = maps$beta, type = "beta")
  psi_gs     <- unflatten_params(flat = psi_vec, meta = flats$psi$meta,idxs = maps$psi, type = "psi") # Psi is an special case

  # Re-input the (possibly) changed psi matrices into the psi_gks matrices.
  # psi_gks is necessary to compute the loglikelihood. However, some parameters (from the empty group-cluster combinations)
  # are treated as nuisance parameters and are fixed.
  for(g in 1:ngroups){
    psi_mat[[psi_idx[g, 1], psi_idx[g, 2]]] <- psi_gs[[g]]
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

      # S_biased <- S[[g]] * (N_gs[[g]] - 1) / N_gs[[g]]
      S_biased <- S[[g]]

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

#### OBJECTIVE FUNCTION (ONLY STEP 2) ----------------------------------------------------------------------------------
obj.S2 <- function(x,
                   psi_mat,
                   pi_ks,
                   S,
                   nclus,
                   ngroups,
                   N_gs,
                   vec_ind,
                   maps,
                   flats){

  beta_vec   <- x[vec_ind[[1]]]
  psi_vec    <- x[vec_ind[[2]]]
  psi_idx    <- vec_ind[[3]]

  beta_mat   <- unflatten_params(flat = beta_vec, meta = flats$beta$meta, idxs = maps$beta, type = "beta")
  psi_gs     <- unflatten_params(flat = psi_vec, meta = flats$psi$meta,idxs = maps$psi, type = "psi") # Psi is an special case

  # Re-input the (possibly) changed psi matrices into the psi_gks matrices.
  # psi_gks is necessary to compute the loglikelihood. However, some parameters (from the empty group-cluster combinations)
  # are treated as nuisance parameters and are fixed.
  for(g in 1:ngroups){
    psi_mat[[psi_idx[g, 1], psi_idx[g, 2]]] <- psi_gs[[g]]
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
      Sigma[[g, k]] <- solve(I - beta_mat[[k]]) %*% psi_mat[[g, k]] %*% t(solve(I - beta_mat[[k]]))

      # Estimate the loglikelihood
      loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
        sample.mean = rep(0, nrow(beta_mat[[1]])),
        sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
        sample.cov  = S[[g]], # Factor covariance matrix from step 1
        Mu          = rep(0, nrow(beta_mat[[1]])),
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
#### OBJECTIVE FUNCTION (ONLY STEP 1) ----------------------------------------------------------------------------------
obj.S1 <- function(x,
                   psi_mat,
                   pi_ks,
                   S,
                   nclus,
                   ngroups,
                   N_gs,
                   vec_ind,
                   maps,
                   flats){

  lambda_vec <- x[vec_ind[[1]]]
  theta_vec  <- x[vec_ind[[2]]]

  lambda_mat <- unflatten_params(flat = lambda_vec, meta = flats$lambda$meta, idxs = maps$lambda, type = "lambda")
  theta_mat  <- unflatten_params(flat = theta_vec, meta = flats$theta$meta, idxs = maps$theta, type = "theta")

  # compute loglikelihood for all group/cluster combinations
  # Initialize matrices to store loglikelihoods
  loglik_gs <- numeric(ngroups)
  Sigma <- vector("list", ngroups)

  for(g in 1:ngroups){
    # Estimate Sigma (factor covariance matrix of step 2)
    Sigma[[g]] <- lambda_mat[[g]] %*% psi_mat[[g]] %*% t(lambda_mat[[g]]) + theta_mat[[g]]

    # Estimate the loglikelihood
    loglik_g <- lavaan:::lav_mvnorm_loglik_samplestats(
      sample.mean = rep(0, nrow(theta_mat[[1]])),
      sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
      sample.cov  = S[[g]], # Factor covariance matrix from step 1
      Mu          = rep(0, nrow(theta_mat[[1]])),
      Sigma       = Sigma[[g]] # Factor covariance matrix from step 2
    )

    loglik_gs[g] <- loglik_g
  }

  # Get total loglikelihood
  LL <- sum(loglik_gs) # Add the maximum again and then sum them all for total loglikelihood

  return(LL)
}

#### MAPPING -----------------------------------------------------------------------------------------------------------
# Map the the matrices to get correct indices for all parameters
mapping_params <- function(mat_list){
  # First round the estimations to avoid differences in floating points
  rounded <- lapply(mat_list, round, digits = 4)

  # Identify constrained parameters
  # Which parameters are different from 1 and 0 AND are the same for the first two groups?
  idx_constrained   <- rounded[[1]] != 1 & rounded[[1]] != 0 & rounded[[1]] == rounded[[2]]
  idx_unconstrained <- rounded[[1]] != 1 & rounded[[1]] != 0 & rounded[[1]] != rounded[[2]]

  # Where are the ones for lambda?
  idx_marker <- rounded[[1]] == 1

  idxs <- list(
    idx_constrained   = idx_constrained,
    idx_unconstrained = idx_unconstrained,
    idx_marker        = idx_marker
  )

  return(idxs)
}

#### FLATTENING --------------------------------------------------------------------------------------------------------
# Flatten the matrices to a single vector
flatten_params <- function(param_list, idxs, type) {
  # First, remove the constrained parameters from all groups (except group 1)
  param_list_copy <- param_list[-1]
  param_list_copy <- lapply(param_list_copy, function(mat){
    mat[idxs$idx_constrained] <- 0 # Set the constrained parameters to 0. They will be later removed
    return(mat)
  })

  # Go back to the complete list
  param_list[2:length(param_list)] <- param_list_copy

  # Get the "meta-data" to later unflatten the matrices
  meta <- lapply(param_list, function(mat) {
    list(dim = dim(mat), rn = rownames(mat), cn = colnames(mat))
  })

  # Flatten the matrices
  # Define operator for named vectors
  operator <- ifelse(type == "lambda", "=~", ifelse(type == "beta", "~", "~~"))
  vec <- lapply(seq_along(param_list), function(g) {
    mat <- param_list[[g]]
    rn <- rownames(mat)
    cn <- colnames(mat)
    if(type != "beta"){
      nm <- as.vector(outer(rn, cn, function(r, c) paste0(c, operator, r, ".g", g)))
    } else {
      nm <- as.vector(outer(rn, cn, function(r, c) paste0(r, operator, c, ".k", g)))
    }
    v <- as.vector(mat)
    names(v) <- nm
    return(v)
  })

  # Remove the non-relevant values from the vector
  vec <- unlist(vec)
  vec <- vec[vec != 0 & vec != 1]

  return(list(vec = vec, meta = meta))
}

#### UNFLATTENING ------------------------------------------------------------------------------------------------------
# Going back to matrix form
# Reconstructs list of matrices from vector and metadata
unflatten_params <- function(flat, meta, idxs, type) {
  # Extract necessary info
  dims <- meta[[1]]$dim
  rn   <- meta[[1]]$rn
  cn   <- meta[[1]]$cn

  # Create empty matrices with correct dimensions (and names)
  mats <- vector("list", length(meta)) # Length is the number of groups
  mats <- lapply(mats, function(mat){
    mat <- matrix(data = 0,
                   nrow = dims[1],
                   ncol = dims[2])
    colnames(mat) <- cn
    rownames(mat) <- rn
    return(mat)
  })

  # Add marker variable for lambdas
  if(type == "lambda"){
    mats <- lapply(mats, function(mat){
      mat[idxs$idx_marker] <- 1
        return(mat)
    })
  }

  # Identify the constrained parameters in the parameter vector
  # Get which parameters only occur once (the constrained ones)
  clean_names <- sub("\\..*", "", names(flat))
  constrained_names <- names(table(clean_names)[which(table(clean_names) == 1)])
  constrained_idxs  <- clean_names %in% constrained_names

  # Separate the constrained from the unconstrained parameters in the parameter vector
  constrained_vec   <- flat[constrained_idxs]
  unconstrained_vec <- flat[!constrained_idxs]

  # The unconstrained parameters are easier to work with if they are separated per group
  ngroups <- length(meta)
  n_unconstrained      <- length(unconstrained_vec)/ngroups # How many unconstrained parameters per group?
  unconstrained_matrix <- matrix(data  = unconstrained_vec,
                                 nrow  = n_unconstrained)
  unconstrained_list   <- split(unconstrained_matrix, col(unconstrained_matrix))

  # Input the constrained values in all the matrices
  mats <- lapply(mats, function(mat){
    mat[idxs$idx_constrained] <- constrained_vec
    return(mat)
  })

  if (isTRUE(any(idxs$idx_unconstrained))){ # Only do this when we actually have non-invariant parameters
    # Input the unconstrained values in their respective matrix
    mats <- lapply(seq_along(mats), function(i) {
      mat <- mats[[i]]
      vals <- unconstrained_list[[i]]
      mat[idxs$idx_unconstrained] <- vals
      return(mat)
    })
  }

  return(mats)
}

