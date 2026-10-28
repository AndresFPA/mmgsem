#' Standard Errors (SE) wrapper for MMGSEM
#'
#' @description
#' A wrapper to compute the SE of a fitted object by mmgsem(). It obtains the SE by numerically computing the second derivative of the loglikelihood function.
#'
#' @usage compute_se(object, d = 1e-03, naive = FALSE, nuisance = FALSE,
#'   hessian_tol = 1e-05)
#'
#' @param object A fitted model returned by the mmgsem function.
#' @param d Step size for numerical approximation of the derivative. Default is 1e-03.
#' @param naive Logical. If TRUE, returns uncorrected SEs (faster but less accurate). If FALSE (default), applies the correction for step-wise estimation.
#' @param nuisance Logical. If TRUE, includes nuisance parameters (parameters from empty group-cluster combinations). Default is FALSE.
#' @param hessian_tol Relative tolerance for the generalized inverse of the Hessians. Singular values at or below this fraction of the largest are excluded. Must be between 0 and 1; the default is 1e-05.
#'
#' @returns SE_vector: List. Includes several vectors that contain the standard errors per type of parameter (beta, psi, etc.)
#'
#' @export
compute_se <- function(object, d = 1e-03, naive = FALSE, nuisance = FALSE,
                       hessian_tol = 1e-05){

  if (!is.numeric(hessian_tol) || length(hessian_tol) != 1L ||
      !is.finite(hessian_tol) || hessian_tol <= 0 || hessian_tol >= 1) {
    stop("'hessian_tol' must be a finite number strictly between 0 and 1.",
         call. = FALSE)
  }

  # ngroups
  ngroups <- length(object$param$lambda)
  nclus   <- (ncol(object$posteriors) - 1)
  gro_clu   <- ngroups * nclus
  pi_ks   <- colMeans(as.matrix(object$posteriors[, -1, drop = FALSE]))
  MM      <- object$MM

  # Prepare the necessary parameter matrices
  lambda <- object$param$lambda
  theta  <- object$param$theta
  beta   <- object$param$beta_ks
  if (is.matrix(beta)) {
    beta <- list(beta)
  }

  # Psi requires extra peparation
  psi_gks <- object$param$psi_gks

  # We remove the 'nuisance' parameters. i.e., the parameters from the empty group-cluster combinations
  # modal_post includes the group label as its first column; only cluster columns index psi_gks.
  modal_post <- as.matrix(object$modal_post[, -1, drop = FALSE])
  psi_idx <- which(x = modal_post == 1, arr.ind = TRUE) # array indices of the relevant group-cluster combinations

  # Reorder based on group number insted of column number
  correct_idx <- order(psi_idx[,1])
  psi_idx <- psi_idx[correct_idx,]

  psi <- vector("list", ngroups)
  for(g in 1:ngroups){
    psi[[g]] <- psi_gks[[psi_idx[g, 1], psi_idx[g, 2]]]
  }

  # If the user requires the nuisance parameters (theoretically, more correct), include also those parameters.
  # It may be necessary with complicated models that involve many parameters
  # Which ones are NOT the correct cluster?
  psi_empty_idx <- which(x = modal_post != 1, arr.ind = TRUE)

  if(isTRUE(nuisance)){
    # How long must the object be? ngroups*nclus - ngroups (only empty combinations)
    gk_minus  <- gro_clu - ngroups
    psi_empty <- vector("list", gk_minus)

    ## Extract the "empty" group-cluster combinations
    # We also remove the exogenous factors (co)variances.
    # Those parameters are already included in the regular "psi" object (the non-nuisance one)
    # Use beta to identify the exogenous variables
    tmp_beta <- beta[[1]]
    exog  <- rownames(tmp_beta)[rowSums(tmp_beta) == 0]
    endog <- !c(rownames(tmp_beta) %in% exog)

    # Proceed to extract the correct object and remove the undesired variances
    for (gk in seq_len(gk_minus)) {
      psi_empty[[gk]] <- psi_gks[[psi_empty_idx[gk, 1], psi_empty_idx[gk, 2]]] # Correct object
      psi_empty[[gk]][exog, exog] <- 0 # Exog variances become 0
    }

  }

  #################################
  ### Flatten all the matrices! ###
  #################################
  # First, get the necessary mapping matrices
  map_lambda <- mapping_params(mat_list = lambda)
  map_theta  <- mapping_params(mat_list = theta)
  map_beta   <- mapping_params(mat_list = beta)
  map_psi    <- mapping_params(mat_list = psi)

  if(isTRUE(nuisance) && length(psi_empty) > 0){
    map_psi_empty <- mapping_params(mat_list = psi_empty)
  } else {
    map_psi_empty <- NULL
  }

  # Put them together in a list to use when unflattening
  maps <- list(
    lambda    = map_lambda,
    theta     = map_theta,
    beta      = map_beta,
    psi       = map_psi,
    psi_empty = map_psi_empty
  )

  # Add the nuisance argument within the psi map
  maps$psi$nuisance <- nuisance

  # Remove the null entry of psi_empty (when nuisance is FALSE)
  maps <- maps[!sapply(maps, is.null)]

  # Now, flatten all the matrices
  flat_lambda <- flatten_params(param_list = lambda, idxs = map_lambda, type = "lambda")
  flat_theta  <- flatten_params(param_list = theta,  idxs = map_theta,  type = "theta")
  flat_beta   <- flatten_params(param_list = beta,   idxs = map_beta,   type = "beta")
  flat_psi    <- flatten_params(param_list = psi,    idxs = map_psi,    type = "psi")

  if(isTRUE(nuisance) && length(psi_empty) > 0){
    flat_psi_empty <- flatten_params(param_list = psi_empty, idxs = map_psi_empty, type = "psi")
  } else {
    flat_psi_empty <- NULL
  }

  # Put them together in a list to use when unflattening
  flats <- list(
    lambda    = flat_lambda,
    theta     = flat_theta,
    beta      = flat_beta,
    psi       = flat_psi,
    psi_empty = flat_psi_empty
  )

  # Remove the null entry of psi_empty (when nuisance is FALSE)
  flats <- flats[!sapply(flats, is.null)]

  ##############################
  ### Time to compute the SE ###
  ##############################

  ##########################
  # STEP 2 STANDARD ERRORS #
  ##########################
  # Get step 2 parameters into one single parameter vector (NAIVE standard errors)
  ifelse(test = isTRUE(nuisance),
         yes  = param_vec_S2 <- c(flat_beta$vec, flat_psi$vec, flat_psi_empty$vec),
         no   = param_vec_S2 <- c(flat_beta$vec, flat_psi$vec))

  # These are positions in param_vec_S2, whose sections are concatenated beta, psi, then
  # optional nuisance psi. Keep the offsets aligned with that order; matching values can
  # select the wrong section when estimates happen to be equal.
  vec_ind_S2 <- vector("list", 5)
  vec_ind_S2[[1]] <- seq_along(flat_beta$vec)
  vec_ind_S2[[2]] <- length(flat_beta$vec) + seq_along(flat_psi$vec)
  vec_ind_S2[[3]] <- psi_idx # Psi matrix indices for psi_gks
  vec_ind_S2[[4]] <- length(flat_beta$vec) + length(flat_psi$vec) +
    seq_along(flat_psi_empty$vec)
  vec_ind_S2[[5]] <- psi_empty_idx # Psi matrix indices for psi_gks

  # First, let's compute the NAIVE standard errors
  message("Naive Standard error computation (Step 2 parameters):\n")
  HESS.S2 <- compute_hessian(f       = obj.S2,
                             x       = param_vec_S2,
                             d       = d,
                             psi_mat = psi_gks,
                             pi_ks   = pi_ks,
                             S       = object$param$cov_eta,
                             nclus   = nclus,
                             ngroups = ngroups,
                             N_gs    = object$N_gs,
                             vec_ind = vec_ind_S2,
                             maps    = maps,
                             flats   = flats)

  if(naive == TRUE){ # If the user only wants the naive standard error, stop here.
    # Vcov matrix of the step 2 parameters
    Sigma_2 <- invert_hessian(HESS.S2, "Step 2", hessian_tol)
    colnames(Sigma_2) <- rownames(Sigma_2) <- colnames(HESS.S2)

    # vector form
    vector_SE <- covariance_se(Sigma_2, "Step 2")
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
  # Get step 1 parameters into one single parameter vector
  param_vec_S1 <- c(flat_lambda$vec, flat_theta$vec)

  # We need the indices of each type of parameter before we continue (it is to reconstruct the matrices later)
  vec_ind_S1 <- vector("list", 3)
  vec_ind_S1[[1]] <- seq_along(flat_lambda$vec)
  vec_ind_S1[[2]] <- length(flat_lambda$vec) + seq_along(flat_theta$vec)

  # Now, let's compute the standard errors of step 1
  message("Used for the standard error correction (Step 1 parameters):\n")
  HESS.S1 <- compute_hessian(f       = obj.S1,
                             x       = param_vec_S1,
                             d       = d,
                             psi_mat = object$param$cov_eta, # Factor's covariance matrix is necessary to compute the LL of Step 1
                             pi_ks   = pi_ks,
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
  ifelse(test = isTRUE(nuisance),
         yes  = param_vec <- c(flat_lambda$vec, flat_theta$vec, flat_beta$vec, flat_psi$vec, flat_psi_empty$vec),
         no   = param_vec <- c(flat_lambda$vec, flat_theta$vec, flat_beta$vec, flat_psi$vec))

  # These are positions in param_vec, whose sections are concatenated lambda, theta,
  # beta, psi, then optional nuisance psi. The offsets must follow that construction
  # order so each objective function receives the correct parameter slice.
  vec_ind <- vector("list", 7)
  vec_ind[[1]] <- seq_along(flat_lambda$vec)
  vec_ind[[2]] <- length(flat_lambda$vec) + seq_along(flat_theta$vec)
  vec_ind[[3]] <- length(flat_lambda$vec) + length(flat_theta$vec) +
    seq_along(flat_beta$vec)
  vec_ind[[4]] <- length(flat_lambda$vec) + length(flat_theta$vec) +
    length(flat_beta$vec) + seq_along(flat_psi$vec)
  vec_ind[[5]] <- psi_idx # Psi matrix indices for psi_gks
  vec_ind[[6]] <- length(flat_lambda$vec) + length(flat_theta$vec) +
    length(flat_beta$vec) + length(flat_psi$vec) + seq_along(flat_psi_empty$vec)
  vec_ind[[7]] <- psi_empty_idx # Psi matrix indices for psi_gks

  # Now, compute the cross-derivatives
  message("Used for the standard error correction (Cross-derivatives):\n")
  HESS <- compute_hessian(f       = obj,
                          x       = param_vec,
                          d       = d,
                          psi_mat = psi_gks,
                          pi_ks   = pi_ks,
                          S       = object$sample.stats$S,
                          nclus   = nclus,
                          ngroups = ngroups,
                          N_gs    = object$N_gs,
                          vec_ind = vec_ind,
                          maps    = maps,
                          flats   = flats)

  ##############################
  ### SE Correction ############
  ##############################
  # First, we need to get vcov of each step independently
  # The joint likelihood is not maximized by the step-wise estimates. Use each
  # step's own Hessian for its marginal SEs; the joint Hessian is only used below
  # for the cross-derivative correction.
  Sigma_1 <- invert_hessian(HESS.S1, "Step 1", hessian_tol)
  Sigma_2 <- invert_hessian(HESS.S2, "Step 2", hessian_tol)
  step1_SE <- covariance_se(Sigma_1, "Step 1")
  step2_SE <- covariance_se(Sigma_2, "Step 2")
  vector_SE <- setNames(
    c(step1_SE, step2_SE),
    c(colnames(HESS.S1), colnames(HESS.S2))
  )

  # Get indices for parameters of both steps from the JOINT Hessian
  step1.idx <- which(colnames(HESS) %in% c(names(flat_lambda$vec), names(flat_theta$vec)))
  ifelse(test = isTRUE(nuisance),
         yes  = step2.idx <- which(colnames(HESS) %in% c(names(flat_beta$vec), names(flat_psi$vec), names(flat_psi_empty$vec))),
         no   = step2.idx <- which(colnames(HESS) %in% c(names(flat_beta$vec), names(flat_psi$vec))))

  # Get matrix divided
  I_22 <- HESS[step2.idx, step2.idx]
  I_21 <- HESS[step2.idx, step1.idx]
  I_11 <- HESS[step1.idx, step1.idx]

  I_22.inv <- MASS::ginv(I_22, tol = hessian_tol)

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
  vector_corrected <- covariance_se(Sigma_2_corrected, "Corrected Step 2")
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

  # Prepare object for the progress bar
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  i_deriv <- 0

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

    # Update the progress bar
    i_deriv <- i_deriv + 1
    utils::setTxtProgressBar(pb, i_deriv)
  }

  return(H)
}

covariance_se <- function(vcov, label) {
  vcov <- (vcov + t(vcov)) / 2
  variances <- diag(vcov)
  scale <- max(1, abs(variances), na.rm = TRUE)
  roundoff <- sqrt(.Machine$double.eps) * scale
  materially_negative <- !is.na(variances) & variances < -roundoff

  if (any(materially_negative, na.rm = TRUE)) {
    warning(
      sprintf(
        "%s covariance has %d negative variance estimate(s); corresponding standard errors are NA.",
        label,
        sum(materially_negative, na.rm = TRUE)
      ),
      call. = FALSE
    )
  }

  variances[is.finite(variances) & variances < 0 & !materially_negative] <- 0
  standard_errors <- sqrt(variances)
  standard_errors[materially_negative] <- NA_real_
  standard_errors[!is.finite(variances)] <- NA_real_
  standard_errors
}

invert_hessian <- function(hessian, label, tol) {
  hessian <- (hessian + t(hessian)) / 2
  eigenvalues <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
  cutoff <- tol * max(abs(eigenvalues))
  weak_directions <- sum(abs(eigenvalues) <= cutoff)
  nonconcave_directions <- sum(eigenvalues > cutoff)

  if (weak_directions > 0) {
    warning(
      sprintf(
        "%s Hessian has %d near-zero direction(s), excluded from the generalized inverse; associated SEs may be unreliable.",
        label,
        weak_directions
      ),
      call. = FALSE
    )
  }
  if (nonconcave_directions > 0) {
    warning(
      sprintf(
        "%s Hessian has %d positive-curvature direction(s); standard errors may be unreliable.",
        label,
        nonconcave_directions
      ),
      call. = FALSE
    )
  }

  MASS::ginv(-hessian, tol = tol)
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

  lambda_vec    <- x[vec_ind[[1]]]
  theta_vec     <- x[vec_ind[[2]]]
  beta_vec      <- x[vec_ind[[3]]]
  psi_vec       <- x[vec_ind[[4]]]
  psi_idx       <- vec_ind[[5]]
  psi_empty_vec <- x[vec_ind[[6]]]
  psi_empty_idx <- vec_ind[[7]]

  lambda_mat <- unflatten_params(flat = lambda_vec, meta = flats$lambda$meta, idxs = maps$lambda, type = "lambda")
  theta_mat  <- unflatten_params(flat = theta_vec, meta = flats$theta$meta, idxs = maps$theta, type = "theta")
  beta_mat   <- unflatten_params(flat = beta_vec, meta = flats$beta$meta, idxs = maps$beta, type = "beta")
  psi_gs     <- unflatten_params(flat = psi_vec, meta = flats$psi$meta,idxs = maps$psi, type = "psi") # Psi is an special case

  nuisance <- maps$psi$nuisance
  if(isTRUE(nuisance) && length(psi_empty_idx) > 0){
    psi_empty <- unflatten_params(flat = psi_empty_vec, meta = flats$psi_empty$meta, idxs = maps$psi_empty, type = "psi") # Psi is an special case
  }

  # Re-input the (possibly) changed psi matrices into the psi_gks matrices.
  # psi_gks is necessary to compute the loglikelihood. However, some parameters (from the empty group-cluster combinations)
  # are treated as nuisance parameters and are fixed.
  for(g in 1:ngroups){
    if(isTRUE(nuisance)){
      # Input the already modified psi_gs for every row in psi_mat (ONLY when nuisance is TRUE)
      psi_mat[g, ] <- psi_gs[g]
    } else {
      # If nuisance is false all the "empty" elements are treated as fixed parameters
      psi_mat[[psi_idx[g, 1], psi_idx[g, 2]]] <- psi_gs[[g]]
    }
  }

  if(isTRUE(nuisance) && length(psi_empty_idx) > 0){
    for(empty in seq_along(psi_empty)){
      psi_mat[[psi_empty_idx[empty, 1], psi_empty_idx[empty, 2]]][maps$psi_empty$idx_unconstrained] <- psi_empty[[empty]][maps$psi_empty$idx_unconstrained]
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

  beta_vec      <- x[vec_ind[[1]]]
  psi_vec       <- x[vec_ind[[2]]]
  psi_idx       <- vec_ind[[3]]
  psi_empty_vec <- x[vec_ind[[4]]]
  psi_empty_idx <- vec_ind[[5]]

  beta_mat   <- unflatten_params(flat = beta_vec, meta = flats$beta$meta, idxs = maps$beta, type = "beta")
  psi_gs     <- unflatten_params(flat = psi_vec, meta = flats$psi$meta, idxs = maps$psi, type = "psi") # Psi is an special case

  nuisance <- maps$psi$nuisance
  if(isTRUE(nuisance) && length(psi_empty_idx) > 0){
    psi_empty <- unflatten_params(flat = psi_empty_vec, meta = flats$psi_empty$meta, idxs = maps$psi_empty, type = "psi") # Psi is an special case
  }

  # Re-input the (possibly) changed psi matrices into the psi_gks matrices.
  # psi_gks is necessary to compute the loglikelihood. However, some parameters (from the empty group-cluster combinations)
  # are treated as nuisance parameters and are fixed.
  for(g in 1:ngroups){
    if(isTRUE(nuisance)){
      # Input the already modified psi_gs for every row in psi_mat (ONLY when nuisance is TRUE)
      psi_mat[g, ] <- psi_gs[g]
    } else {
      # If nuisance is false all the "empty" elements are treated as fixed parameters
      psi_mat[[psi_idx[g, 1], psi_idx[g, 2]]] <- psi_gs[[g]]
    }
  }

  if(isTRUE(nuisance) && length(psi_empty_idx) > 0){
    for(empty in seq_along(psi_empty)){
      psi_mat[[psi_empty_idx[empty, 1], psi_empty_idx[empty, 2]]][maps$psi_empty$idx_unconstrained] <- psi_empty[[empty]][maps$psi_empty$idx_unconstrained]
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

  # Identify free parameters and whether they are shared across all matrices.
  idx_free <- Reduce(`|`, lapply(rounded, function(mat) {
    !is.na(mat) & mat != 0 & mat != 1
  }))
  idx_shared <- matrix(TRUE, nrow = nrow(rounded[[1]]), ncol = ncol(rounded[[1]]))
  if (length(rounded) > 1) {
    for (mat in rounded[-1]) {
      idx_shared <- idx_shared & !is.na(mat) & !is.na(rounded[[1]]) &
        mat == rounded[[1]]
    }
  }
  idx_constrained   <- idx_free & idx_shared
  idx_unconstrained <- idx_free & !idx_shared

  idxs <- list(
    idx_constrained   = idx_constrained,
    idx_unconstrained = idx_unconstrained
  )

  return(idxs)
}

#### FLATTENING --------------------------------------------------------------------------------------------------------
# Flatten the matrices to a single vector
flatten_params <- function(param_list, idxs, type) {
  # Get the "meta-data" to later unflatten the matrices
  meta <- lapply(param_list, function(mat) {
    list(dim = dim(mat), rn = rownames(mat), cn = colnames(mat), values = mat)
  })

  # Flatten the matrices
  # Define operator for named vectors
  operator <- ifelse(type == "lambda", "=~", ifelse(type == "beta", "~", "~~"))
  vec <- lapply(seq_along(param_list), function(g) {
    mat <- param_list[[g]]
    rn <- rownames(mat)
    cn <- colnames(mat)
    idx_parameter <- idxs$idx_unconstrained
    if (g == 1) {
      idx_parameter <- idx_parameter | idxs$idx_constrained
    }
    if (type %in% c("theta", "psi")) {
      idx_parameter <- idx_parameter & lower.tri(mat, diag = TRUE)
    }
    if(type != "beta"){
      nm <- as.vector(outer(rn, cn, function(r, c) paste0(c, operator, r, ".g", g)))
    } else {
      nm <- as.vector(outer(rn, cn, function(r, c) paste0(r, operator, c, ".k", g)))
    }
    v <- as.vector(mat)
    names(v) <- nm
    return(v[as.vector(idx_parameter)])
  })

  vec <- unlist(vec, use.names = TRUE)

  return(list(vec = vec, meta = meta))
}

#### UNFLATTENING ------------------------------------------------------------------------------------------------------
# Going back to matrix form
# Reconstructs list of matrices from vector and metadata
unflatten_params <- function(flat, meta, idxs, type) {
  # Start from the fitted matrices so fixed entries (including unit variances) survive
  # reconstruction; the flattened vector contains only free parameters.
  mats <- lapply(meta, `[[`, "values")

  operator <- ifelse(type == "lambda", "=~", ifelse(type == "beta", "~", "~~"))
  symmetric <- type %in% c("theta", "psi")
  idx_constrained <- idxs$idx_constrained
  idx_unconstrained <- idxs$idx_unconstrained
  if (symmetric) {
    lower_triangle <- lower.tri(mats[[1]], diag = TRUE)
    idx_constrained <- idx_constrained & lower_triangle
    idx_unconstrained <- idx_unconstrained & lower_triangle
  }
  parameter_names <- function(mat, i) {
    rn <- rownames(mat)
    cn <- colnames(mat)
    if (type != "beta") {
      as.vector(outer(rn, cn, function(r, c) paste0(c, operator, r, ".g", i)))
    } else {
      as.vector(outer(rn, cn, function(r, c) paste0(r, operator, c, ".k", i)))
    }
  }

  # Separate the constrained from the unconstrained parameters in the parameter vector
  constrained_values <- flat[match(
    parameter_names(mats[[1]], 1)[as.vector(idx_constrained)],
    names(flat)
  )]
  if (anyNA(constrained_values)) {
    stop("Could not match constrained parameters while reconstructing matrices.")
  }
  unconstrained_values <- lapply(seq_along(mats), function(i) {
    flat[match(
      parameter_names(mats[[i]], i)[as.vector(idx_unconstrained)],
      names(flat)
    )]
  })
  if (any(vapply(unconstrained_values, anyNA, logical(1)))) {
    stop("Could not match unconstrained parameters while reconstructing matrices.")
  }

  for (i in seq_along(mats)) {
    if (isTRUE(any(idx_constrained))) {
      mats[[i]][idx_constrained] <- constrained_values
    }
    if (isTRUE(any(idx_unconstrained))) {
      mats[[i]][idx_unconstrained] <- unconstrained_values[[i]]
    }
    if (symmetric) {
      mats[[i]][upper.tri(mats[[i]])] <- t(mats[[i]])[upper.tri(mats[[i]])]
    }
  }

  return(mats)
}
