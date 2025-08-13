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
compute_se <- function(object, d = 1e-03){
  # ngroups
  ngroups <- length(object$param$lambda)

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

  ### Flatten all the matrices! ###
  # First, get the necessary mapping matrices
  map_lambda <- mapping_params(mat_list = lambda)
  map_theta  <- mapping_params(mat_list = theta)
  map_beta   <- mapping_params(mat_list = beta)
  map_psi    <- mapping_params(mat_list = psi)

  # Now, flatten all the matrices
  flat_lambda <- flatten_params(param_list = lambda, idxs = map_lambda, type = "lambda")
  flat_theta  <- flatten_params(param_list = theta,  idxs = map_theta,  type = "theta")
  flat_beta   <- flatten_params(param_list = beta,   idxs = map_beta,   type = "beta")
  flat_psi    <- flatten_params(param_list = psi ,   idxs = map_psi,    type = "psi")

  # Get everything into one single parameter vector
  param_vec <- c(flat_lambda$vec, flat_theta$vec, flat_beta$vec, flat_psi$vec)

  # We need the indices of each type of parameter before we continue (it is to reconstruct the matrices later)
  vec_ind <- vector("list", 4)
  vec_ind[[1]] <- which(param_vec %in% flat_lambda$vec)
  vec_ind[[2]] <- which(param_vec %in% flat_theta$vec)
  vec_ind[[3]] <- which(param_vec %in% flat_beta$vec)
  vec_ind[[4]] <- which(param_vec %in% flat_psi$vec)

  ### Time to compute the SE ###
  HESS <- compute_hessian(f = obj,
                          x,
                          d,
                          psi_mat,
                          pi_ks,
                          S,
                          nclus,
                          ngroups,
                          N_gs,
                          vec_ind)

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
                            vec_ind) {

  # Prepare the matrices necessary to store the results
  n  <- length(x)
  H  <- matrix(0, n, n)
  colnames(H) <- rownames(H) <- names(x)
  f1 <- matrix(0, n)

  # Loglikelihood result without adding the small number
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
#### OBJECTIVE FUNCTION ---------------------------------------------------------------------------------------------------
# The function that computes the loglikelihood
obj <- function(x,
                psi_mat,
                pi_ks,
                S,
                nclus,
                ngroups,
                N_gs,
                vec_ind){

  lambda_vec <- x[vec_ind[[1]]]
  theta_vec  <- x[vec_ind[[2]]]
  beta_vec   <- x[vec_ind[[3]]]
  psi_vec    <- x[vec_ind[[4]]]

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
unflatten_params <- function(flat, idxs, type) {
  # Extract necessary info
  dims <- flat$meta[[1]]$dim
  rn   <- flat$meta[[1]]$rn
  cn   <- flat$meta[[1]]$cn

  # Create empty matrices with correct dimensions (and names)
  mats <- vector("list", length(flat$meta)) # Length is the number of groups
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
  clean_names <- sub("\\..*", "", names(flat$vec))
  constrained_names <- names(table(clean_names)[which(table(clean_names) == 1)])
  constrained_idxs  <- clean_names %in% constrained_names

  # Separate the constrained from the unconstrained parameters in the parameter vector
  constrained_vec   <- flat$vec[constrained_idxs]
  unconstrained_vec <- flat$vec[!constrained_idxs]

  # The unconstrained parameters are easier to work with if they are separated per group
  ngroups <- length(flat$meta)
  n_unconstrained      <- length(unconstrained_vec)/ngroups # How many unconstrained parameters per group?
  unconstrained_matrix <- matrix(data  = unconstrained_vec,
                                 nrow  = n_unconstrained)
  unconstrained_list   <- split(unconstrained_matrix, col(unconstrained_matrix))

  # Input the constrained values in all the matrices
  mats <- lapply(mats, function(mat){
    mat[idxs$idx_constrained] <- constrained_vec
    return(mat)
  })

  # Input the unconstrained values in their respective matrix
  mats <- lapply(seq_along(mats), function(i) {
    mat <- mats[[i]]
    vals <- unconstrained_list[[i]]
    mat[idxs$idx_unconstrained] <- vals
    return(mat)
  })

  return(mats)
}

