#' Mixture Multi-Group Structural Equation Modelling (MMGSEM)
#'
#' Performs a mixture clustering based on the structural parameters (i.e., regressions) of a SEM model. 
#' The estimation is done in a step-wise fashion and uses an expectation-maximization (EM) algorithm in the second step.
#'
#' INPUT: Arguments required by the function
#' @param dat Observed data of interest for the MMGSEM model.
#' @param S1 Measurement model (MM). Used in step 1. Must be a string (like in lavaan). 
#'                   Can be a list of strings determing the number of measurement blocks (e.g., one string for the MM of 
#'                   factor 1, and a second string for the MM of factor 2)  
#' @param S2 Structural model. Used in step 2. Must be a string (like in lavaan).
#' @param group Name of the group variable. Must be a string.
#' @param nclus Pre-specified number of clusters.
#' @param seed Pre-defined seed for the random start in case a replication is needed in the future.
#' @param userStart Pre-defined start provided by the user. Must be a matrix of dimensions ngroups*nclus.
#'                  A 1 represents that the row group is in the column cluster. There must be only one 1
#'                  in each row. The function will return results only for this start.
#'                  Example for 6 groups and 2 clusters:
#'
#'                                   [,1] [,2]
#'                             [1,]    1    0
#'                             [2,]    1    0
#'                             [3,]    1    0
#'                             [4,]    0    1
#'                             [5,]    0    1
#'                             [6,]    0    1
#'
#' @param s1out Resulting lavaan object from a cfa() analysis. Can be used to directly input the results from the step 1
#'              if the user only wants to use MMGSEM() to estimate the structural model (Step 2). If not NULL, MMGSEM()
#'              will skip the estimation of Step 1 and use the s1out object as the input for Step 2.
#' @param max_it Maximum number of iterations for each start.
#' @param nstarts Number of random starts that will be performed (in an attempt to avoid local maxima).
#' @param partition Type of random partition used to start the EM algorithm. Can be "hard" and "soft".
#' @param constraints String containing the constrains of the model. Receives the same input as group.equal in lavaan.
#' @param NonInv String containing the non-invariances of the model. Similar to group.partial in lavaan.
#' @param Endo2Cov TRUE or FALSE argument to determine whether to allow or not covariance between endogenous 2 variables.
#'                 If TRUE (the default), the covariance between endogenous factors is allowed.
#' @param allG TRUE or FALSE. Determines whether the endogenous covariances are group (T) or cluster (F) specific.
#'             By default, is TRUE.
#' @param fit String argument. Either "factors" or "observed". Determines which loglikelihood will be used to fit the model.
#' @param est_method either "local" or "global. Follows local and global approaches from the SAM method. GLOBAL NOT FUNCTIONAL YET.
#'
#' OUTPUT: The function will return a list with the following results:
#' @return z_gks: Posterior classification probabilities or cluster memberships (ngroups*nclus matrix).
#' @return final_fit: Lavaan fit of the best and final model (not useful, only for testing purposes).
#' @return psi_gks: Resulting psi matrices (i.e., residual factor covariances) for each group-cluster combination.
#' @return lambda: Resulting invariant lambda matrix (i.e., factor loadings) for all groups.
#' @return theta_gs: Resulting theta matrices (i.e., items' residual/unique variances) for each group.
#' @return beta_ks: Resulting beta matrices (i.e., regression coefficients) for each cluster.
#' @return loglikelihood: Total loglikelihood of the best model.
#' @return logliks_starts: Total loglikelihood of the models of each run (allows to check sensitivity to local maxima).
#' @return Obs_loglik: Reconstruction of the observed loglikelihood after the model has converged. Only valid when
#'                     fit is "factors".
#'
#' PLEASE NOTE: This function requires 'lavaan' package to work. It also requires the function 'EStep.R'.
#' 
#' @export
MMGSEM <- function(dat, S1 = NULL, S2 = NULL,
                   group, nclus, seed = NULL, userStart = NULL, s1_fit = NULL,
                   max_it = 10000L, nstarts = 20L, printing = FALSE,
                   partition = "hard", Endo2Cov = TRUE, allG = TRUE, fit = "factors",
                   est_method = "local", meanstr = FALSE, ...) {
  #browser()
  
  # Get arguments in ...
  dots_args <- list(...)
  constraints <- dots_args$group.equal
  
  # Add a warning in case there is a pre-defined start and the user also requires a multi-start
  if (!(is.null(userStart)) && nstarts > 1) {
    warning("If a start is defined by the user, no multi-start is performed. The results correspond to the one start used an input")
    nstarts <- 1
  }

  
  # Get several values relevant for future steps
  g_name  <- as.character(unique(dat[, group]))
  vars    <- lavaan::lavNames(lavaanify(S1, auto = TRUE))
  lat_var <- lavaan::lavNames(lavaanify(S1, auto = TRUE), "lv")
  n_var   <- length(vars)

  # # Center the data per group (so that the mean for all variables in each group is 0)
  centered <- dat
  
  # if the mean structure is not required, then remove the mean structure of the data (i.e., center the data)
  # but, if the intercepts are required, then meanstr changes to TRUE
  if(isFALSE(meanstr) & "intercepts" %in% constraints){
    warning("If the intercepts are included in the constraints, then meanstr automatically changes to TRUE to include the mean structure.")
    meanstr <- T
  }
  
  if(isFALSE(meanstr)){ 
    group.idx <- match(dat[, group], g_name)
    group.sizes <- tabulate(group.idx)
    group.means <- rowsum.default(as.matrix(dat[, vars]),
                                  group = group.idx, reorder = FALSE,
                                  na.rm = FALSE
    ) / group.sizes
    centered[, vars] <- dat[, vars] - group.means[group.idx, , drop = FALSE]
  } 
  
  
  # Get sample covariance matrix per group (used later)
  S_unbiased <- lapply(X = unique(centered[, group]), FUN = function(x) {
    cov(centered[centered[, group] == x, vars])
  })
  
  ## STEP 1 - MMG-SEM ----------------------------------------------------------------------------------------
  # Save the measurement model results
  # Call function to run Step 1 of MMG-SEM (estimates CFA)
  Step1_args <- list(S1 = S1, s1_fit = s1_fit, centered = centered, 
                     group = group, S_unbiased = S_unbiased)
  Step1_args <- c(dots_args, Step1_args)
  MM <- do.call(what = Step1, args = Step1_args)
  
  # Extract necessary objects
  ngroups   <- MM$ngroups
  S1output  <- MM$S1output
  lambda_gs <- MM$lambda_gs
  theta_gs  <- MM$theta_gs
  cov_eta   <- MM$cov_eta
  N_gs      <- MM$N_gs
  S_biased  <- MM$S_biased
  
  gro_clu   <- ngroups * nclus
  
  # STEP 2 (EM algorithm for model estimation) -----------------------------------------------------
  # We perform a MULTI-START procedure to avoid local maxima.
  # Initialize objects to store results per random start.

  results_nstarts <- vector(mode = "list", length = nstarts)
  z_gks_nstarts <- vector(mode = "list", length = nstarts) # z_gks refer to posteriors
  loglik_nstarts <- numeric(nstarts)

  # Start using a pre-defined seed for the random partitions
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Create a function to reorder matrices (used later). Done due to:
  # (1) It allows us to organize the matrices in an easier to understand order. Exogenous variables first, and endogeonus later.
  # (2) Used to make sure we are multiplying matrices in the correct way
  
  # Re-order for factors
  reorder <- function(x) {
    x <- x[c(exog, endog), c(exog, endog)]
    return(x)
  }

  # Re-order for measurement model matrices (observed variables must be in the same order as the factors)
  reorder_obs <- function(x, matrix) {
    # Re-write model
    # Split into lines
    lines_model <- unlist(strsplit(unlist(S1), "\n"))
    rewritten <- c()
    # browser()
    # Extract the observed variable names per factor
    # Exogenous factors
    for (i in 1:length(exog)) {
      rewritten <- c(rewritten, lines_model[grepl(exog[i], lines_model)])
    }

    # Endogenous 1 factors
    if(length(endog1 > 0)){
      for (i in 1:length(endog1)) {
        rewritten <- c(rewritten, lines_model[grepl(endog1[i], lines_model)])
      }
    }

    # Endogenous 2 factors
    for (i in 1:length(endog2)) {
      rewritten <- c(rewritten, lines_model[grepl(endog2[i], lines_model)])
    }

    # Run fake measur_model
    fake_measur  <- lavaan::cfa(model = rewritten, data = dat, do.fit = F)
    correct_vars <- lavNames(fake_measur)

    if (matrix == "lambda") {
      # Reorder lambda
      x <- x[correct_vars, c(exog, endog)]
    } else if (matrix == "theta") {
      # Reorder theta
      x <- x[correct_vars, correct_vars]
    }

    return(x)
  }

  # Do a fake sem() to obtain the correct settings to use in Step 2
  # just a single sample cov!
  fake <- lavaan::sem(
    model = S2, sample.cov = rep(cov_eta[1], nclus),
    sample.nobs = rep(nrow(dat), nclus), do.fit = FALSE,
    baseline = FALSE,
    h1 = FALSE, check.post = FALSE,
    loglik = FALSE,
    sample.cov.rescale = FALSE,
    fixed.x = TRUE,
    information = "observed"
  )
  FakeprTbl <- parTable(fake)
  fake@Options$do.fit <- TRUE
  fake@Options$se     <- "none"
  fake@ParTable$start <- NULL
  fake@ParTable$est   <- NULL
  fake@ParTable$se    <- NULL
  fake@Options$start  <- "default"

  # Get the labels of the endogenous 1 and 2 factors
  endog1 <- lat_var[(lat_var %in% FakeprTbl$rhs[which(FakeprTbl$op == "~")]) &
                      (lat_var %in% FakeprTbl$lhs[which(FakeprTbl$op == "~")])]
  endog2 <- lat_var[!c(lat_var %in% FakeprTbl$rhs[which(FakeprTbl$op == "~")]) &
                      (lat_var %in% FakeprTbl$lhs[which(FakeprTbl$op == "~")])]
  endog <- c(endog1, endog2)
  exog <- lat_var[!c(lat_var %in% endog)]
  
  # endog2 <- lavNames(fake, "ov.y")
  # endog1 <- lavNames(fake, "eqs.y")
  # both.idx <- which(endog1 %in% endog2)
  # if (length(both.idx) > 0L) {
  #   endog1 <- endog1[-both.idx]
  # }
  # endog <- c(endog1, endog2)
  # exog <- lavNames(fake, "lv")[!c(lavNames(fake, "lv") %in% endog)]

  # Do a fake model per endo LV (to avoid bias due to reconstruction of group-specific endo variances). For more info, please see the Appendix of the paper related to this code.
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
    fac_load        <- which(FakeprTbl$lhs != this_lv & FakeprTbl$op == "=~") # Factor load of other factors

    # Get the new parameter table per endogenous latent variable
    prTbl_idx       <- c(which(FakeprTbl$lhs == this_lv), fac_load, var_not_this_lv)
    prTbl_idx       <- sort(prTbl_idx)
    prTbl_lv[[lv]]  <- FakeprTbl[prTbl_idx, ]

    # Run the model per endo latent variable
    fake_lv[[lv]] <- lavaan::sem(
      model = prTbl_lv[[lv]], sample.cov = rep(cov_eta[1], nclus),
      sample.nobs = rep(nrow(dat), nclus), do.fit = FALSE,
      baseline = FALSE,
      h1 = FALSE, check.post = FALSE,
      loglik = FALSE,
      sample.cov.rescale = FALSE,
      fixed.x = TRUE,
      information = "observed"
    )
    
    fake_lv[[lv]]@Options$do.fit <- TRUE
    fake_lv[[lv]]@Options$se     <- "none"
    fake_lv[[lv]]@ParTable$start <- NULL
    fake_lv[[lv]]@ParTable$est   <- NULL
    fake_lv[[lv]]@ParTable$se    <- NULL
    fake_lv[[lv]]@Options$start  <- "default"
  }
 
  # Re-order (order of columns and rows) cov_eta to make sure later computations are comparing correct matrices
  cov_eta <- lapply(1:ngroups, function(x) {
    reorder(cov_eta[[x]])
  })
  lambda_gs <- lapply(1:ngroups, function(x) {
    reorder_obs(lambda_gs[[x]], matrix = "lambda")
  })
  theta_gs <- lapply(1:ngroups, function(x) {
    reorder_obs(theta_gs[[x]], matrix = "theta")
  })
  S_unbiased <- lapply(1:ngroups, function(x) {
    reorder_obs(S_unbiased[[x]], matrix = "theta")
  })

  # Multi-start
  for (s in 1:nstarts) {
    if (printing == T) {
      print(paste("Start", s, "-----------------"))
    }

    # Random Start
    if (!is.null(userStart)) {
      # In case the user inputs a pre-defined start, use it for z_gks
      z_gks <- userStart
    } else if (partition == "hard") {
      # Create initial random partition. Hard partition. (z_gks)
      cl <- 0
      while (cl < 1) { # "while loop" to make sure all clusters get at least one group
        z_gks <- t(replicate(ngroups, sample(x = c(rep(0, (nclus - 1)), 1))))
        cl <- min(colSums(z_gks))
      }
    } else if (partition == "soft") {
      z_gks <- matrix(data = runif(n = c(nclus * ngroups)), ncol = nclus, nrow = ngroups)
      z_gks <- z_gks / rowSums(z_gks)
    }

    # Initialize psi_gks
    psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)

    # Initialize weighted z_gks (weighted for each endo LV)
    # Done to avoid bias - NOT NECESSARY in the first iteration
    z_gks_lv <- vector(mode = "list", length = length(endog))
    for (lv in 1:length(endog)) {
      z_gks_lv[[lv]] <- matrix(data = NA, ncol = nclus, nrow = ngroups)
    }

    # Prepare objects for the while loop
    i <- 0 # iteration initialization
    prev_LL <- 0 # previous loglikelihood initialization
    diff_LL <- 1 # Set a diff of 1 just to start the while loop
    log_test <- T # TEMPORARY - To check if there is decreasing loglikelihood

    # Run full-convergence multi-start
    while (diff_LL > 1e-6 & i < max_it & isTRUE(log_test)) {
      i <- i + 1
      pi_ks <- colMeans(z_gks) # Prior probabilities
      N_gks <- z_gks * N_gs # Sample size per group-cluster combination
      # N_gks <- as.vector(N_gks)

      # Weight posteriors for each endogenous factor - Not necessary in the first iteration
      # To avoid bias when allG == T. That is, when the endogenous variances are group-specific
      if (isTRUE(allG) & i > 1) {
        for (lv in 1:length(endog)) {
          for (k in 1:nclus) {
            for (g in 1:ngroups) {
              # Correct bias by dividing by the correct endo LV
              z_gks_lv[[lv]][g, k] <- z_gks[g, k] / psi_gks[[g, k]][endog[lv], endog[lv]]
            }
          }
        }
      }

      # M-Step ---------

      # Trick to avoid slow multi-group estimation
      # Get a weighted averaged covariance matrix for each cluster
      if (isFALSE(allG) | i == 1) {
        
        # Do this when allG is False OR when it is True and we are in the first iteration
        # For the first iteration there is no weighted z_gks
        COV <- vector("list", length = nclus)

        for (k in 1:nclus) {
          # create 'averaged' sample cov for this cluster
          this_nobs <- z_gks[, k] * N_gs
          this_w <- this_nobs / sum(this_nobs)
          tmp <- lapply(seq_along(cov_eta), function(g) {
            cov_eta[[g]] * this_w[g]
          })
          COV[[k]] <- Reduce("+", tmp)
          # COV[[k]] <- 0.5 * (COV[[k]] + t(COV[[k]])) # Force the matrix to be symmetric
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
            this_nobs <- z_gks_lv[[lv]][, k] * N_gs
            this_w <- this_nobs / sum(this_nobs)
            tmp <- lapply(seq_along(cov_eta), function(g) {
              cov_eta[[g]] * this_w[g]
            })
            COV_lv[[lv]][[k]] <- Reduce("+", tmp)
            # COV_lv[[lv]][[k]] <- 0.5 * (COV_lv[[lv]][[k]] + t(COV_lv[[lv]][[k]])) # Force symmetry
          }
        }
      }

      # PARAMETER ESTIMATION
      # Call lavaan to estimate the structural parameters
      # the 'groups' are the clusters
      # Note: this makes all resulting parameters to be cluster-specific (it is reconstructed later)

      if (isFALSE(allG) | i == 1) {
        # Do this when allG is False OR when it is True and we are in the first iteration
        # For the first iteration, perform the full structural model estimation
        s2out <- lavaan::lavaan(
          slotOptions       = fake@Options,
          slotParTable      = fake@ParTable,
          sample.cov        = COV,
          sample.nobs       = rep(nrow(dat), nclus),
          # slotModel       = slotModel,
          # slotData        = fake@Data,
          # slotSampleStats = fake@SampleStats
        )
      } else if (i > 1) {
        # After the first iteration
        # Run structural estimation once per endo LV
        s2out <- vector(mode = "list", length = length(endog))
        for (lv in 1:length(endog)) {
          s2out[[lv]] <- lavaan(
            slotOptions       = fake_lv[[lv]]@Options,
            slotParTable      = fake_lv[[lv]]@ParTable,
            sample.cov        = COV_lv[[lv]],
            sample.nobs       = rep(nrow(dat), nclus)
            # slotModel       = slotModel,
            # slotData        = fake@Data,
            # slotSampleStats = fake@SampleStats
          )
        }
      }


      # start <- partable(s2out)$est

      # compute loglikelihood for all group/cluster combinations
      # Initialize matrices to store loglikelihoods
      loglik_gks  <- matrix(data = 0, nrow = ngroups, ncol = nclus)
      loglik_gksw <- matrix(data = 0, nrow = ngroups, ncol = nclus)
      gk <- 0

      # Prepare Sigma
      # Initialize the object for estimating Sigma
      Sigma <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
      I <- diag(length(lat_var)) # Identity matrix based on number of latent variables. Used later

      # Extract cluster-specific parameters from step 2
      if (isFALSE(allG) | i == 1) {
        # Do this when allG is False OR when it is True and we are in the first iteration
        # In iteration one, there is only model (s2out) from which we extract the parameters.
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
          beta_ks <- reorder(beta_ks)
          psi_ks <- reorder(psi_ks)
        } else if (nclus != 1) {
          beta_ks <- lapply(1:nclus, function(x) {
            reorder(beta_ks[[x]])
          }) # Does not work with only one cluster
          psi_ks <- lapply(1:nclus, function(x) {
            reorder(psi_ks[[x]])
          })
        }
      } else if (i > 1) {
        # After iteration 1 we have several models (s2out[[lv]]) from which we extract the parameters
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
            col.idx <- colnames(beta_ks_lv[[lv]][[k]])

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
          beta_ks <- reorder(beta_ks)
        } else if (nclus != 1) {
          beta_ks <- lapply(1:nclus, function(x) {
            reorder(beta_ks[[x]])
          }) # Does not work with only one cluster
        }
      }

      for (k in 1:nclus) {
        # Previous matrices were only cluster-specific. We have to reconstruct the group-specific matrices (psi and sigma)

        ## Save the cluster-specific psi and beta
        ## ifelse() in case of only 1 cluster
        ifelse(test = (nclus == 1), yes = (psi <- psi_ks), no = (psi <- psi_ks[[k]]))
        ifelse(test = (nclus == 1), yes = (beta <- beta_ks), no = (beta <- beta_ks[[k]]))
        for (g in 1:ngroups) {
          # Reconstruct psi and sigma so they are group- and cluster-specific again.
          # Replace the group-specific part of psi
          # Exogenous (co)variance is always group-specific
          psi[exog, exog] <- cov_eta[[g]][exog, exog] # Replace the group-specific part

          # If the user required group-specific endogenous covariances (allG = T), do:
          if (allG == T) {
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
          }

          # If required by the user, set to 0 the covariance between endog 2 factors
          if (isFALSE(Endo2Cov)) {
            psi[endog2, endog2][row(psi[endog2, endog2]) != col(psi[endog2, endog2])] <- 0
          }

          # Store for future check
          psi_gks[[g, k]] <- psi

          # Get log-likelihood by comparing factor covariance matrix of step 1 (cov_eta) and step 2 (Sigma)
          if (fit == "factors") {
            # Estimate Sigma (factor covariance matrix of step 2)
            Sigma[[g, k]] <- solve(I - beta) %*% psi %*% t(solve(I - beta))
            Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]])) # Force to be symmetric
            
            # Estimate the loglikelihood
            loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
              sample.mean = rep(0, nrow(cov_eta[[1]])),
              sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
              # sample.nobs = N_gks[g, k],
              sample.cov  = cov_eta[[g]], # Factor covariance matrix from step 1
              Mu          = rep(0, nrow(cov_eta[[1]])),
              Sigma       = Sigma[[g, k]] # Factor covariance matrix from step 2
            )
          } else if (fit == "observed") {
            # browser()
            S_biased <- S_unbiased[[g]] * (N_gs[[g]] - 1) / N_gs[[g]]
            Sigma[[g, k]] <- lambda_gs[[g]] %*% solve(I - beta) %*% psi %*% t(solve(I - beta)) %*% t(lambda_gs[[g]]) + theta_gs[[g]]
            Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]]))
            # Sigma[[g, k]][lower.tri(Sigma[[g, k]])] <- t(Sigma[[g, k]])[lower.tri(Sigma[[g, k]])]
            loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
              sample.mean = rep(0, length(vars)),
              sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
              sample.cov  = S_biased, # Item (observed) covariance matrix from step 1
              Mu          = rep(0, length(vars)),
              Sigma       = Sigma[[g, k]] # Item (observed) covariance matrix from step 2
            )
          }
          loglik_gks[g, k] <- loglik_gk
          loglik_gksw[g, k] <- log(pi_ks[k]) + loglik_gk # weighted loglik
        } # ngroups
      } # cluster

      # Get total loglikelihood
      # First, deal with arithmetic underflow by subtracting the maximum value per group
      max_gs <- apply(loglik_gksw, 1, max) # Get max value per row
      minus_max <- sweep(x = loglik_gksw, MARGIN = 1, STATS = max_gs, FUN = "-") # Subtract the max per row
      exp_loglik <- exp(minus_max) # Exp before summing for total loglikelihood
      loglik_gsw <- log(apply(exp_loglik, 1, sum)) # Sum exp_loglik per row and then take the log again
      LL <- sum((loglik_gsw + max_gs)) # Add the maximum again and then sum them all for total loglikelihood

      # Now, do E-step
      E_out <- EStep(
        pi_ks = pi_ks, ngroup = ngroups,
        nclus = nclus, loglik = loglik_gks
      )

      z_gks <- E_out
      diff_LL <- abs(LL - prev_LL)
      log_test <- prev_LL < LL | isTRUE(all.equal(prev_LL, LL))
      if (i == 1) {
        log_test <- T
      }
      if (log_test == F) {
        print(paste("Start", s, "; Iteration", i, "-------"))
        print(paste("Difference", LL - prev_LL))
      } # ; browser()}
      log_test <- T
      prev_LL <- LL
      if (printing == T) {
        print(i)
        print(LL)
      }
    }

    results_nstarts[[s]] <- s2out
    z_gks_nstarts[[s]] <- z_gks
    loglik_nstarts[s] <- LL
  } # multistart

  # Get best fit and z_gks based on the loglikelihood
  best_idx <- which.max(loglik_nstarts)
  s2out    <- results_nstarts[[best_idx]]
  LL       <- loglik_nstarts[best_idx]
  z_gks    <- z_gks_nstarts[[best_idx]]
  colnames(z_gks) <- paste("Cluster", seq_len(nclus))

  # Extract matrices from final step 2 output
  if (isFALSE(allG)) {
    if (nclus == 1) {
      EST_s2  <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE) # Estimated matrices step 2
      beta_ks <- EST_s2[["beta"]]
      psi_ks  <- EST_s2[["psi"]]
    } else if (nclus != 1) {
      EST_s2  <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE) # Estimated matrices step 2
      beta_ks <- lapply(EST_s2, "[[", "beta")
      psi_ks  <- lapply(EST_s2, "[[", "psi")
    }
  } else if (isTRUE(allG)) {
    EST_s2_lv <- vector(mode = "list", length = length(endog))
    beta_ks_lv <- vector(mode = "list", length = length(endog))
    psi_ks_lv <- vector(mode = "list", length = length(endog))
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

    # Re-construct beta_ks
    for (k in 1:nclus) {
      for (lv in 1:length(endog)) {
        this_lv <- endog[lv]
        col.idx <- colnames(beta_ks_lv[[lv]][[k]])
        if (nclus == 1) {
          beta_ks[this_lv, col.idx] <- beta_ks_lv[[lv]][this_lv, col.idx]
        } else if (nclus != 1) {
          beta_ks[[k]][this_lv, col.idx] <- beta_ks_lv[[lv]][[k]][this_lv, col.idx]
        }
      }
    }
  }

  # Re-order betas
  if (nclus == 1) {
    beta_ks <- reorder(beta_ks)
  } else if (nclus != 1) {
    beta_ks <- lapply(1:nclus, function(x) {
      reorder(beta_ks[[x]])
    }) # Does not work with only one cluster
  }

  # Get the group- and cluster-specific psi_gks matrices
  psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)

  for (k in 1:nclus) {
    ifelse(test = (nclus == 1), yes = (psi_k <- psi_ks), no = (psi_k <- psi_ks[[k]]))
    ifelse(test = (nclus == 1), yes = (beta <- beta_ks), no = (beta <- beta_ks[[k]]))
    for (g in 1:ngroups) {
      psi_gks[[g, k]] <- psi_k
      psi_gks[[g, k]][exog, exog] <- cov_eta[[g]][exog, exog]

      # If the user required group-specific endogenous covariances (allG = T), do:
      if (allG == T) {
        # Take into account the effect of the cluster-specific regressions
        # cov_eta[[g]] = solve(I - beta) %*% psi %*% t(solve(I - beta))
        # If we solve for psi, then:
        g_endog1_cov <- ((I - beta) %*% cov_eta[[g]] %*% t((I - beta)))[endog1, endog1] # Extract group-specific endog cov
        if (length(endog1) > 1) {
          g_endog1_cov[row(g_endog1_cov) != col(g_endog1_cov)] <- 0
        }
        # browser()
        psi_gks[[g, k]][endog1, endog1] <- g_endog1_cov

        g_endog2_cov <- ((I - beta) %*% cov_eta[[g]] %*% t((I - beta)))[endog2, endog2] # Extract group-specific endog cov
        psi_gks[[g, k]][endog2, endog2] <- g_endog2_cov
      }

      # If required by the user, endog2 covariances are set to 0
      if (isFALSE(Endo2Cov)) {
        offdiag <- row(psi_gks[[g, k]][endog2, endog2]) != col(psi_gks[[g, k]][endog2, endog2])
        psi_gks[[g, k]][endog2, endog2][offdiag] <- 0
      }
    } # groups
  } # cluster

  ############################
  ##### GLOBAL ESTIMATION ####
  ############################
  # If the user requires global estimation of SAM:
  #   - Use the final local estimates to start global SAM.
  #   - Prepare the parameter table
  if(est_method == "global"){
    # Create a dummy complete parameter table.
    fake_S <- rep(S_unbiased, nclus) # Duplicate observed covariances to match the number of clusters
    names(fake_S) <- paste("group", seq_len(gro_clu))
    fake_global   <- lavaan::parTable(lavaan::sem(model = c(S1, S2), 
                                                  sample.cov = fake_S,
                                                  sample.nobs = rep(N_gs, nclus), 
                                                  do.fit = FALSE,
                                                  meanstructure = F,
                                                  ceq.simple = T))
    
    # Update parameter table with only needed free parameters
    fake_global$free    <- 0 # Set all parameters to be not freely estimated 
    # fake_global$free[fake_global$rhs %in% lat_var] <- 1:length(which(fake_global$rhs %in% lat_var))
    fake_global$par     <- paste0(fake_global$lhs, fake_global$op, fake_global$rhs, ".g", fake_global$group)
    fake_global$cluster <- rep(1:nclus, each = length(fake_global$id[fake_global$group %in% 1:ngroups])) # Add cluster column
    
    
    # Fill the measurement model parameters in the fake_global table
    # Get parameter table of the measurement model
    if(is.list(S1)){ # If we have measurement blocks...
      s1table <- c()
      for(m in 1:M){
        s1table <- rbind(s1table, lavaan::parTable(S1output[[m]]))
      }
    } else {
      s1table <- lavaan::partable(S1output)
    }
    
    # Remove the structural parameters of the s1 parameter table
    idx.str <- which(s1table$rhs %in% lat_var)
    s1table <- s1table[-idx.str, ]
    
    # Extend the measurement par table to have enough group-cluster combinations
    s1table   <- s1table[-which(s1table$op == "=="), ] # Remove constrain rows (they don't matter anymore, mm is fixed)
    nr.par.mm <- nrow(s1table)/ngroups
    s1table   <- s1table[rep(seq_len(nrow(s1table)), nclus), ]
    s1table$group <- rep(x = 1:gro_clu, each = nr.par.mm)
    s1table$par <- paste0(s1table$lhs, s1table$op, s1table$rhs, ".g", s1table$group)
    
    # FIX THE MM PARAMETERS IN FAKE_GLOBAL #
    mm.idx <- match(s1table$par, fake_global$par)
    fake_global$ustart[mm.idx] <- s1table$est # Measurement parameters are now fixed
    
    # Fill in, the structural parameters in the table (not fixed, only as starting values)
    # Get structural model parameter table
    fake_global$parK <- paste0(fake_global$lhs, fake_global$op, fake_global$rhs, ".k", fake_global$cluster)
    
    # Add the estimates of the structural parameters as start of the fake_global table
    # Extract the regression parameters in a vector
    # Regressions (Cluster-specific!)
    beta_vec <- c() # Empty vector for the betas
    beta_nam <- c() # Empty vector for the names of such betas (not truly necessary, but makes everything easier to understand)
    
    # The loop identifies the non-zero parameters (i.e., the free parameters) in the beta matrices
    # It also saves the number of the parameter using the col and row names of the beta matrices (e.g., F3 ~ F1)
    for(k in 1:nclus){
      ifelse(test = (nclus == 1), yes = (beta <- beta_ks), no = (beta <- beta_ks[[k]]))
      non_zer.idx <- which(unlist(beta) != 0)
      beta_vec <- c(beta_vec, unlist(beta)[non_zer.idx])
      beta_nam <- c(beta_nam, as.vector(outer(X = rownames(beta), 
                                              Y = colnames(beta), 
                                              function(x, y) paste0(x, "~", y, ".k", k)))[non_zer.idx])
    }
    
    beta_vec <- setNames(object = beta_vec, nm = beta_nam) # Name the vector (the .1 and .2 represent the cluster)
    beta.idx <- match(fake_global$parK, names(beta_vec)) # Indices of beta_vec in fake_global
    fake_global$ustart <- ifelse(test = !is.na(beta.idx), yes = beta_vec[beta.idx], no = fake_global$ustart)
    
    # Factor covariances (Group-cluster specific!)
    cov_vec <- c() # Empty vector for the covariances
    cov_nam <- c() # Empty vector for the names of such covs
    
    # The loop is similar the one of the betas.
    # Note there is an extra index (i.e., unique.idx). Given that the cov matrix is symmetric, some parameters are repeated.
    # However, even if repeated, it is only one parameter. The unique.idx removes the repeated parameter.
    gk <- 0
    for(k in 1:nclus){
      for(g in 1:ngroups){
        gk <- gk + 1
        non_zer.idx <- which(unlist(psi_gks[[g, k]]) != 0)
        unique.idx  <- !duplicated(unlist(psi_gks[[g, k]])[non_zer.idx])
        cov_vec <- c(cov_vec, unlist(psi_gks[[g, k]])[non_zer.idx][unique.idx])
        cov_nam <- c(cov_nam, as.vector(outer(X = rownames(psi_gks[[g, k]]),
                                              Y = colnames(psi_gks[[g, k]]),
                                              function(x, y) paste0(y, "~~", x, ".g", gk)))[non_zer.idx][unique.idx])
      }
    }
    
    cov_vec <- setNames(object = cov_vec, nm = cov_nam) # Name the vector (the number after the . is the group)
    cov.idx <- match(fake_global$par, names(cov_vec)) # Indices of cov_vec in fake_global
    fake_global$ustart <- ifelse(test = !is.na(cov.idx), yes = cov_vec[cov.idx], no = fake_global$ustart)
    
    # Remove unnecessary columns from fake_global
    fake_global$se      <- NULL
    fake_global$cluster <- NULL
    fake_global$par     <- NULL
    fake_global$parK    <- NULL
    fake_global$est     <- NULL
    fake_global$start   <- NULL
    
    # Add constrain to the regressions per cluster (ceq.simple approach) - NOT WORKING!
    # n_reg <- length(fake_global$id[which(fake_global$group == 1 & fake_global$op == "~")])
    # for(k in 1:nclus){
    #   beta.idx <- which(fake_global$op == "~" & fake_global$cluster == k) # Get index of regression parameters
    #   fake_global$free[beta.idx] <- (n_reg*(k-1)+1):(n_reg*k) # Add constrain to free column
    #   fake_global$label[beta.idx] <- fake_global$plabel[which(fake_global$op == "~" & fake_global$cluster == k & fake_global$group == 1 + (ngroups*(k-1)))]
    # }
    # 
    # # Add constrain to the exogenous covariances per group
    # max_free <- max(fake_global$free)
    # n_cov_exo <- ((length(exog) * (length(exog) + 1)) / 2)
    # for(g in 1:ngroups){
    #   exo.idx <- which(fake_global$op == "~~" &
    #                    fake_global$group %in% c(g, (ngroups)*(1:nclus) + g) &
    #                    fake_global$lhs %in% exog) # Get index of exo cov parameters
    # 
    #   fake_global$free[exo.idx] <- (max_free+(n_cov_exo*(g-1))+1):(max_free+(n_cov_exo*g)) # Add constrain to free column
    # 
    #   fake_global$label[exo.idx] <- fake_global$plabel[which(fake_global$op == "~~" &
    #                                                          fake_global$group == g &
    #                                                          fake_global$lhs %in% exog)]
    # }
    # 
    # # Add free (unconstrained) estimation for the group-cluster-specific endog covariances
    # max_free <- max(fake_global$free)
    # free_idx <- which(fake_global$free == 0 & fake_global$op == "~~" & fake_global$lhs %in% endog)
    # fake_global$free[free_idx] <- (max_free+1):(max_free+length(fake_global$free[free_idx]))
    
    # Add free parameters
    free_idx <- which(fake_global$rhs %in% lat_var)
    fake_global$free[free_idx] <- 1:length(free_idx)
    
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
    constraints <- fake_global$plabel[which(fake_global$op == "~")] # Get the regression parameters
    n_reg <- length(fake_global$plabel[which(fake_global$op == "~" & fake_global$group == 1)]) # Number of reg PER GROUP
    
    
    # Identify latent variables that are both independent and dependent variables in the model
    # Labels of the variance parameters of variables in exo
    cons_exo <- fake_global$plabel[which(fake_global$op == "~~" & fake_global$lhs %in% exog)]
    # Number of variance parameters that involve variables in exo PER GROUP
    n_exo <- length(fake_global$plabel[which(fake_global$op == "~~" &
                                               fake_global$lhs %in% exog & fake_global$group == 1)])
    
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
    clus_label <- rep(x = 1:nclus, each = ngroups)
    group_label <- rep(x = 1:ngroups, times = nclus)
    
    # Add cluster labels to the parameter table (not necessary, just for me)
    for (j in 1:length(clus_label)) {
      fake_global$cluster[fake_global$group == j] <- clus_label[j]
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
      constraints_matrix[(reg_labels == k), "rhs"] <- cluster_par # On the right hand side insert parameters of all groups
    }
    
    for(g in 1:ngroups){
      # Variances constraints (exo)
      group_par_exo <- cons_exo[exo_labels == g]
      cons_exo_matrix[(exo_labels == g), "lhs"] <- group_par_exo[1:n_exo]
      cons_exo_matrix[(exo_labels == g), "rhs"] <- group_par_exo
    }
    
    # Remove redundant constraints (i.e. p1 == p1)
    constraints_total <- rbind(constraints_matrix, cons_exo_matrix)
    redundant <- which(constraints_total$lhs == constraints_total$rhs) # Identify redundant
    constraints_total <- constraints_total[-redundant, ]
    rownames(constraints_total) <- NULL
    
    # Bind the model table with the constraints
    # browser()
    fake_global <- rbind(fake_global, constraints_total)
    
    # Update fake_global with correct order of parameter numbers in the free column
    
    # Now, actually start the global estimation
    # Step 2: Iterative EM Algorithm ----
    # Create initial random probabilities (z_gks)
    
    # Start using a pre-defined seed for the partition
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    i <- 0 # iteration initialization
    prev_LL <- 0 # previous loglikelihood initialization
    diff_LL <- 1 # Set a diff of 1 just to start the while loop
    
    while (diff_LL > 1e-6 & i < max_it) {
      i <- i + 1
      pi_ks <- colMeans(z_gks)
      N_gks <- z_gks * N_gs # Sample size per group-cluster combination
      N_gks <- c(N_gks)
      # browser()
      # M-Step
      if (i == 1) {
        s2out <- lavaan::sem(
            model = fake_global, sample.cov = rep(S_biased, nclus), sample.nobs = N_gks,
            baseline = FALSE, se = "none",
            h1 = FALSE, check.post = FALSE,
            control = list(rel.tol = 1e-06),
            sample.cov.rescale = FALSE,
            fixed.x = FALSE, ceq.simple = T
        )
      } else {
        # fake_global$ustart <- NA
        s2out <- lavaan::sem(
          model = fake_global, sample.cov = rep(S_biased, nclus), sample.nobs = N_gks,
          baseline = FALSE,
          h1 = FALSE, check.post = FALSE,
          control = list(rel.tol = 1e-06), start = start,
          sample.cov.rescale = FALSE,
          fixed.x = FALSE, ceq.simple = T
        ) # , control = list(max.iter = 50))
      }
      
      start <- coef(s2out)
      
      # E-Step
      # Prepare the log likelihood used in the E-step
      # Get the log likelihood from s2out results
      # 1. Estimate log-likelihood using Lavaan
      
      gk <- 0
      
      global_loglik_gks <- matrix(data = 0, nrow = ngroups, ncol = nclus)
      global_loglik_gksw <- matrix(data = 0, nrow = ngroups, ncol = nclus)
      global_LL <- 0
      
      for (k in 1:nclus) {
        for (g in 1:ngroups) {
          # browser()
          gk <- gk + 1L
          global_loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
            sample.mean = rep(0, length(vars)), # s2out@SampleStats@mean[[gk]],
            sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
            sample.cov  = S_biased[[g]], # s2out@SampleStats@cov[[gk]],
            Mu          = rep(0, length(vars)), # s2out@SampleStats@mean[[gk]],
            Sigma       = s2out@implied$cov[[gk]]
          )
          global_loglik_gks[g, k] <- global_loglik_gk
          global_loglik_gksw[g, k] <- log(pi_ks[k]) + global_loglik_gk # weighted loglik
        }
      }
      
      # Get total loglikelihood
      # First, deal with arithmetic underflow by subtracting the maximum value per group
      max_gs <- apply(global_loglik_gksw, 1, max) # Get max value per row
      minus_max <- sweep(x = global_loglik_gksw, MARGIN = 1, STATS = max_gs, FUN = "-") # Subtract the max per row
      exp_loglik <- exp(minus_max) # Exp before summing for total loglikelihood
      loglik_gsw <- log(apply(exp_loglik, 1, sum)) # Sum exp_loglik per row and then take the log again
      global_LL <- sum((loglik_gsw + max_gs)) # Add the maximum again and then sum them all for total loglikelihood
      
      # Now, do E-step
      E_out <- EStep(
        pi_ks = pi_ks, ngroup = ngroups,
        nclus = nclus, loglik = global_loglik_gks
      )
      # browser()
      z_gks <- E_out
      diff_LL <- abs(global_LL - prev_LL)
      prev_LL <- global_LL
      # print(LL)
      # print(diff_LL)
    }
    
    # Extract matrices from final step 2 output
    EST_s2      <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE) # Estimated matrices step 2
    beta_gks    <- lapply(EST_s2, "[[", "beta")
    psi_gks_tmp <- lapply(EST_s2, "[[", "psi")
    
    # Select useful betas
    k.idx <- (seq_len(nclus) - 1) * ngroups + 1L
    beta_ks <- beta_gks[k.idx]
    
    # Re-order betas
    if (nclus == 1) {
      beta_ks <- beta_ks[[1]]
      beta_ks <- reorder(beta_ks)
    } else if (nclus != 1) {
      beta_ks <- lapply(1:nclus, function(x) {
        reorder(beta_ks[[x]])
      }) # Does not work with only one cluster
    }
    
    # Re-order psi
    # Put them in a matrix of matrices
    psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
    for(gk in 1:(ngroups*nclus)){
      psi_gks[[gk]] <- psi_gks_tmp[[gk]]
    }
    
  }
  
  # MODEL SELECTION
  # Get observed data log-likelihood for model selection purposes)
  Sigma_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
  Obs.loglik_gks <- matrix(data = 0, nrow = ngroups, ncol = nclus)
  Obs.loglik_gksw <- matrix(data = 0, nrow = ngroups, ncol = nclus)
  pi_ks <- colMeans(z_gks)

  for (k in 1:nclus) {
    ifelse(test = (nclus == 1), yes = (beta <- beta_ks), no = (beta <- beta_ks[[k]]))
    for (g in 1:ngroups) {
      S_biased <- S_unbiased[[g]] * (N_gs[[g]] - 1) / N_gs[[g]]
      var_eta <- solve(I - beta) %*% psi_gks[[g, k]] %*% t(solve(I - beta))
      Sigma_gks[[g, k]] <- lambda_gs[[g]] %*% var_eta %*% t(lambda_gs[[g]]) + theta_gs[[g]]
      Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]]))
      # Sigma[[g, k]][lower.tri(Sigma[[g, k]])] <- t(Sigma[[g, k]])[lower.tri(Sigma[[g, k]])]
      Obs.loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
        sample.mean = rep(0, length(vars)),
        sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
        # sample.nobs = N_gks[g, k],
        sample.cov  = S_biased, # Item (observed) covariance matrix from step 1
        Mu          = rep(0, length(vars)),
        Sigma       = Sigma_gks[[g, k]] # Item (observed) covariance matrix from step 2
      )

      Obs.loglik_gks[g, k] <- Obs.loglik_gk
      Obs.loglik_gksw[g, k] <- log(pi_ks[k]) + Obs.loglik_gk
    }
  }

  # Get total observed loglikelihood
  # First, deal with arithmetic underflow by subtracting the maximum value per group
  Obs.max_gs <- apply(Obs.loglik_gksw, 1, max) # Get max value per row
  Obs.minus_max <- sweep(x = Obs.loglik_gksw, MARGIN = 1, STATS = Obs.max_gs, FUN = "-") # Subtract the max per row
  Obs.exp_loglik <- exp(Obs.minus_max) # Exp before summing for total loglikelihood
  Obs.loglik_gsw <- log(apply(Obs.exp_loglik, 1, sum)) # Sum exp_loglik per row and then take the log again
  Obs.LL <- sum((Obs.loglik_gsw + Obs.max_gs)) # Add the maximum again and then sum them all for total loglikelihood

  # Calculate BIC for model selection
  # Four types of BIC:
  # (1) BIC(N) based on the log-likelihood from the factors
  # (2) BIC(G) based on the log-likelihood from the factors
  # (3) BIC(N) based on the log-likelihood from the observed data
  # (4) BIC(G) based on the log-likelihood from the observed data

  # Get important values
  Q <- length(lat_var)
  J <- length(vars)

  # Structural parameters
  ifelse(test = (nclus == 1), yes = (n_reg <- sum(beta_ks != 0)), no = (n_reg <- sum(beta_ks[[1]] != 0)))
  Q_exo <- length(exog)
  n_cov_exo <- ((Q_exo * (Q_exo + 1)) / 2)
  Q_endo1 <- length(endog1)
  Q_endo2 <- length(endog2)
  n_cov_endo2 <- ((Q_endo2 * (Q_endo2 + 1)) / 2)

  # Measurement parameters
  n_res <- sum(theta_gs[[1]] != 0)
  n_load <- sum(lambda_gs[[1]] != 0)

  # How many free loadings?
  # Identify the free loadings using the parameter table from lavaan
  n_free <- 0
  if (is.list(S1)) {
    for (m in 1:m) {
      partbl    <- lavaan::parTable(S1output[[m]])
      free_load <- which(partbl$op == "=~" & is.na(partbl$ustart) & partbl$group == 1 & partbl$label != partbl$plabel)
      n_free    <- n_free + length(free_load)
    }
  } else if (!is.list(S1)) {
    partbl    <- lavaan::parTable(S1output)
    free_load <- which(partbl$op == "=~" & is.na(partbl$ustart) & partbl$group == 1 & partbl$label != partbl$plabel)
    n_free    <- length(free_load)
  }

  # Get the correct number of free parameters depending on the possible combinations
  # browser()
  if (allG == F) { # Is endogenous covariance group-specific?
    nr_par_factors <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus)
    nr_pars <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus) + (n_res * ngroups) + (n_load - Q - n_free) + (n_free * ngroups)
  } else if (allG == T) {
    nr_par_factors <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * ngroups) + (n_cov_endo2 * ngroups)
    nr_pars <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * ngroups) + (n_cov_endo2 * ngroups) + (n_res * ngroups) + (n_load - Q - n_free) + (n_free * ngroups)
  }

  # Calculate BIC
  # Observed
  Obs.BIC_N <- (-2 * Obs.LL) + (nr_pars * log(sum(N_gs)))
  Obs.BIC_G <- (-2 * Obs.LL) + (nr_pars * log(ngroups))

  # Factors
  BIC_N <- (-2 * LL) + (nr_par_factors * log(sum(N_gs)))
  BIC_G <- (-2 * LL) + (nr_par_factors * log(ngroups))
  
  # Calculate AIC (and AIC3). 
  # Observed
  Obs.AIC  <- (-2 * Obs.LL) + (nr_pars * 2)
  Obs.AIC3 <- (-2 * Obs.LL) + (nr_pars * 3)
  
  # Factors
  AIC  <- (-2 * LL) + (nr_par_factors * 2)
  AIC3 <- (-2 * LL) + (nr_par_factors * 3)
  
  # Calculate entropy and ICL
  # Entropy
  # Code from github user daob (Oberski, 2019): https://gist.github.com/daob/c2b6d83815ddd57cde3cebfdc2c267b3
  # p is the prior or posterior probabilities
  entropy <- function(p) {
      p <- p[p > sqrt(.Machine$double.eps)] # since Lim_{p->0} p log(p) = 0
      sum(-p * log(p))
  }
  # browser()
  sum_entropy <- sum(apply(z_gks, 1, entropy)) # Total entropy
  
  # Entropy R2
  entropy.R2 <- function(pi, post) { 
      error_prior <- entropy(pi) # Class proportions
      error_post <- mean(apply(post, 1, entropy))
      R2_entropy <- (error_prior - error_post) / error_prior
      R2_entropy
  }
  
  R2_entropy <- entropy.R2(pi = pi_ks, post = z_gks)
  # browser()
  # ICL
  ICL     <- BIC_G + (sum_entropy * 2)
  Obs.ICL <- Obs.BIC_G + (sum_entropy * 2)
  
  # Re order matrices so that we get them in the following order:
  # (1) Exogenous latent variables
  # (2) Endogenous latent variables: independent and dependent variables at the same time
  # (3) Endogenous latent variables: only dependent variables

  # Reoder psi_ks and beta_ks by using the reorder function in the lapply function
  psi_gks <- array(lapply(1:gro_clu, function(x) {
    reorder(psi_gks[[x]])
  }), dim = c(ngroups, nclus))
  if (nclus == 1) {
    beta_ks <- reorder(beta_ks)
  } else if (nclus != 1) {
    beta_ks <- lapply(1:nclus, function(x) {
      reorder(beta_ks[[x]])
    }) # Does not work with only one cluster
  }

  names(beta_ks) <- paste("Cluster", seq_len(nclus))

  return(list(
    posteriors    = z_gks,
    final_fit     = s2out, # Final fit of step 2 (contains all group-cluster combinations)
    MM            = S1output, # Output of step 1 (measurement model)
    param         = list(psi_gks = psi_gks, lambda = lambda_gs, # Lambda is invariant across all groups
                         theta = theta_gs, beta_ks = beta_ks, cov_eta = cov_eta), # Factor covariance matrix from first step
    logLik        = list(loglik        = LL, # Final logLik of the model (its meaning depends on arguments "fit" and "est_method")
                         global_loglik = ifelse(test = est_method == "global", global_LL, NA), # Only valid if est_method = "global"
                         loglik_gksw   = loglik_gksw, # Weighted logLik per group-cluster combinations
                         runs_loglik   = loglik_nstarts, # loglik for each start
                         obs_loglik    = Obs.LL), # Only useful if fit = "local"
    model_sel     = list(BIC        = list(observed = list(BIC_N = Obs.BIC_N, BIC_G = Obs.BIC_G),
                                           Factors = list(BIC_N = BIC_N, BIC_G = BIC_G)),
                         AIC        = list(observed = Obs.AIC, Factors = AIC),
                         AIC3       = list(observed = Obs.AIC3, Factors = AIC3),
                         R2_entropy = R2_entropy,
                         ICL        = list(observed = Obs.ICL, Factors = ICL)),
    sample.stats  = list(S = S_unbiased, n_cov_exo = n_cov_exo),
    NrPar         = list(Obs.nrpar = nr_pars, Fac.nrpar = nr_par_factors),
    N_gs          = N_gs
  ))
}


# Step 1 function ----------------------------------------------------------------------------------

Step1 <- function(S1 = S1, s1_fit = s1_fit, centered = centered,
                  group = group, S_unbiased = S_unbiased, ...){
  # Step 1: Get group-specific factor covariances
  # Perform Step 1 according to the number of measurement blocks
  
  if (is.list(S1)) { # Do we have measurement blocks?
    M <- length(S1) # How many measurement blocks?
    
    if (!is.null(s1_fit)) {
      # If the user inputs their own step 1 results, use it
      S1output <- s1_fit
    } else if (is.null(s1_fit)) {
      # If not, estimate step 1 using cfa()
      S1output <- vector(mode = "list", length = length(S1))
      for (m in 1:M) {
        # Estimate one cfa per measurement block
        S1output[[m]] <- lavaan::cfa(
          model = S1[[m]], data = centered, group = group,
          estimator = "ML", se = "none", test = "none", 
          baseline = FALSE, h1 = FALSE,
          implied = FALSE, loglik = FALSE,
          meanstructure = FALSE, ...
        )
      }
    }
    
    # How many groups?
    ngroups <- lavInspect(S1output[[1]], "ngroups")
    
    # Extract measurement parameters per measurement block
    # Extract Lambda & Theta for each group in all blocks
    # Initialize lists to store lambdas and thetas per block
    lambda_block <- vector(mode = "list", length = M)
    theta_block  <- vector(mode = "list", length = M)
    for (m in 1:M) {
      EST_block         <- lavaan::lavInspect(S1output[[m]], "est")
      lambda_block[[m]] <- lapply(X = EST_block, "[[", "lambda")
      theta_block[[m]]  <- lapply(X = EST_block, "[[", "theta")
    }
    
    # Put together lambda & theta for all groups
    # We should end with one lambda and theta matrix per group
    lambda_group <- vector(mode = "list", length = ngroups)
    theta_group  <- vector(mode = "list", length = ngroups)
    for (g in 1:ngroups) {
      for (m in 1:M) { # Put matrices of the same group in the same list
        lambda_group[[g]][[m]] <- lambda_block[[m]][[g]]
        theta_group[[g]][[m]]  <- theta_block[[m]][[g]]
      }
      
      # Put together the matrices per group
      # Lambda
      lambda_group[[g]] <- lavaan::lav_matrix_bdiag(lambda_group[[g]])
      
      # Theta
      theta_group[[g]]  <- lavaan::lav_matrix_bdiag(theta_group[[g]])
      
      # Label correctly the rows and columns of the resulting matrices
      # Lambda
      rownames(lambda_group[[g]]) <- vars
      colnames(lambda_group[[g]]) <- lat_var
      
      # Theta
      rownames(theta_group[[g]]) <- colnames(theta_group[[g]]) <- vars
    }
    
    # Change names and get matrices/values relevant for future steps
    lambda_gs <- lambda_group
    theta_gs  <- theta_group
    N_gs      <- lavaan::lavInspect(S1output[[1]], "nobs") # nobs per group
    
    # Estimate cov_eta (Covariance between the factors)
    M_mat   <- vector(mode = "list", length = ngroups) # M matrices from SAM
    cov_eta <- vector(mode = "list", length = ngroups)
    
    for (g in 1:ngroups) {
      # Compute the M (mapping) matrix in case we have different blocks
      lambda_g <- lambda_gs[[g]]
      theta_g  <- theta_gs[[g]]
      M_mat[[g]] <- solve(t(lambda_g) %*% solve(theta_g) %*% lambda_g) %*% t(lambda_g) %*% solve(theta_g)
      
      # Get the covariance of the factors (cov_eta)
      # First, get biased sample covariance matrix per group (S)
      S <- S_unbiased[[g]] * (N_gs[[g]] - 1) / N_gs[[g]]
      cov_eta[[g]] <- M_mat[[g]] %*% (S - theta_g) %*% t(M_mat[[g]])
    }
  } else if (!is.list(S1)) {
    # If not a list, then we only have one measurement block (all latent variables at the same time)
    if (!is.null(s1_fit)) {
      # If the user input their own step 1 results, use it
      S1output <- s1_fit
    } else if (is.null(s1_fit)) {
      S1output <- lavaan::cfa(
        model = S1, data = centered, group = group,
        estimator = "ML", se = "none", test = "none", 
        baseline = FALSE, h1 = FALSE,
        implied = FALSE, loglik = FALSE,
        meanstructure = FALSE, ...
      )
    }
    
    # Define some important objects
    # How many groups?
    ngroups <- lavaan::lavInspect(S1output, "ngroups")
    N_gs    <- lavaan::lavInspect(S1output, "nobs") # nobs per group
    
    # all estimated model matrices, per group
    EST       <- lavaan::lavInspect(S1output, "est", add.class = FALSE, add.labels = TRUE)
    theta_gs  <- lapply(EST, "[[", "theta")
    lambda_gs <- lapply(EST, "[[", "lambda")
    cov_eta   <- lapply(EST, "[[", "psi") # cov_eta name refers to Variance of eta (eta being the latent variables)
  }
  
  # Biased cov matrix
  S_biased <- vector(mode = "list", length = ngroups)
  for(g in 1:ngroups){S_biased[[g]] <- S_unbiased[[g]] * (N_gs[[g]] - 1) / N_gs[[g]]}
  
  # Return all the useful objects
  return(list(
    S1output  = S1output,
    lambda_gs = lambda_gs,
    theta_gs  = theta_gs,
    cov_eta   = cov_eta,
    ngroups   = ngroups,
    N_gs      = N_gs,
    S_biased  = S_biased
  ))
}

# E-step --------------------------------------------------------
EStep <- function(pi_ks, ngroup, nclus, loglik){
  
  max_g <-rep(0,ngroup)
  z_gks <- matrix(NA,nrow = ngroup,ncol = nclus)
  
  for(g in 1:ngroup){
    for(k in 1:nclus){
      z_gks[g,k] <- log(pi_ks[k])+loglik[g,k]
    }
    max_g[g] <- max(z_gks[g,]) # prevent arithmetic underflow 
    z_gks[g,] <- exp(z_gks[g,]-rep(max_g[g],nclus))
  }
  
  # divide by the rowwise sum of the above calculated part 
  z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
  # z_gks <- round(z_gks, digits = 16)
  # z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
  
  return(z_gks)
}

