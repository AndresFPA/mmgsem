#' Mixture Multi-Group Structural Equation Modelling (MMGSEM)
#'
#' Performs a mixture clustering based on the structural parameters (i.e., regressions) of a SEM model. 
#' The estimation is done in a step-wise fashion and uses an expectation-maximization (EM) algorithm in the second step.
#'
#' INPUT: Arguments required by the function
#' @param dat Observed data of interest for the MMGSEM model.
#' @param step1model Measurement model (MM). Used in step 1. Must be a string (like in lavaan). 
#'                   Can be a list of strings determing the number of measurement blocks (e.g., one string for the MM of 
#'                   factor 1, and a second string for the MM of factor 2)  
#' @param step2model Structural model. Used in step 2. Must be a string (like in lavaan).
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
#' ################################# ALL GLOBAL SAM CODE IS NOT FUNCTIONAL YET #########################################

MMGSEM <- function(dat, step1model = NULL, step2model = NULL,
                   group, nclus, seed = NULL, userStart = NULL, s1out = NULL,
                   max_it = 10000L, nstarts = 20L, printing = FALSE,
                   partition = "hard", NonInv = NULL, constraints = "loadings",
                   Endo2Cov = TRUE, allG = TRUE, fit = "factors", 
                   se = "standard", est_method = "local", meanstr = FALSE,
                   do.se = F) {
  
  # Add a warning in case there is a pre-defined start and the user also requires a multi-start
  if (!(is.null(userStart)) && nstarts > 1) {
    warning("If a start is defined by the user, no multi-start is performed. The results correspond to the one start used an input")
    nstarts <- 1
  }

  # Get several values relevant for future steps
  g_name <- as.character(unique(dat[, group]))
  vars <- lavNames(lavaanify(step1model, auto = TRUE))
  lat_var <- lavNames(lavaanify(step1model, auto = TRUE), "lv")
  n_var <- length(vars)

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
    centered[, vars] <- dat[, vars] - group.means[group.idx, drop = FALSE]
  } 
  
  # Get sample covariance matrix per group (used later)
  S_unbiased <- lapply(X = unique(centered[, group]), FUN = function(x) {
    cov(centered[centered[, group] == x, vars])
  })

  ## STEP 1 - MMG-SEM ----------------------------------------------------------------------------------------

  # Step 1: Get group-specific factor covariances
  # Perform Step 1 according to the number of measurement blocks
  
  if (is.list(step1model)) { # Do we have measurement blocks?
    M <- length(step1model) # How many measurement blocks?

    if (!is.null(s1out)) {
      # If the user inputs their own step 1 results, use it
      S1output <- s1out
    } else if (is.null(s1out)) {
      # If not, estimate step 1 using cfa()
      S1output <- vector(mode = "list", length = length(step1model))
      for (m in 1:M) {
        # Estimate one cfa per measurement block
        S1output[[m]] <- lavaan::cfa(
          model = step1model[[m]], data = centered, group = group,
          estimator = "ML", group.equal = constraints,
          se = "none", test = "none", baseline = FALSE, h1 = FALSE,
          implied = FALSE, loglik = FALSE,
          meanstructure = FALSE, group.partial = NonInv
        )
      }
    }

    # How many groups?
    ngroups <- lavInspect(S1output[[1]], "ngroups")

    # Extract measurement parameters per measurement block
    # Extract Lambda & Theta for each group in all blocks
    # Initialize lists to store lambdas and thetas per block
    lambda_block <- vector(mode = "list", length = M)
    theta_block <- vector(mode = "list", length = M)
    for (m in 1:M) {
      EST_block <- lavInspect(S1output[[m]], "est")
      lambda_block[[m]] <- lapply(X = EST_block, "[[", "lambda")
      theta_block[[m]] <- lapply(X = EST_block, "[[", "theta")
    }

    # Put together lambda & theta for all groups
    # We should end with one lambda and theta matrix per group
    lambda_group <- vector(mode = "list", length = ngroups)
    theta_group <- vector(mode = "list", length = ngroups)
    for (g in 1:ngroups) {
      for (m in 1:M) { # Put matrices of the same group in the same list
        lambda_group[[g]][[m]] <- lambda_block[[m]][[g]]
        theta_group[[g]][[m]] <- theta_block[[m]][[g]]
      }
      # Put together the matrices per group
      # Lambda
      lambda_group[[g]] <- lav_matrix_bdiag(lambda_group[[g]])

      # Theta
      theta_group[[g]] <- lav_matrix_bdiag(theta_group[[g]])

      # Label correctly the rows and columns of the resulting matrices
      # Lambda
      rownames(lambda_group[[g]]) <- vars
      colnames(lambda_group[[g]]) <- lat_var

      # Theta
      rownames(theta_group[[g]]) <- colnames(theta_group[[g]]) <- vars
    }

    # Change names and get matrices/values relevant for future steps
    lambda_gs <- lambda_group
    theta_gs <- theta_group
    N_gs <- lavInspect(S1output[[1]], "nobs") # nobs per group

    # Estimate cov_eta (Covariance between the factors)
    M_mat <- vector(mode = "list", length = ngroups) # M matrices from SAM
    cov_eta <- vector(mode = "list", length = ngroups)

    for (g in 1:ngroups) {
      # Compute the M (mapping) matrix in case we have different blocks
      lambda_g <- lambda_gs[[g]]
      theta_g <- theta_gs[[g]]
      M_mat[[g]] <- solve(t(lambda_g) %*% solve(theta_g) %*% lambda_g) %*% t(lambda_g) %*% solve(theta_g)

      # Get the covariance of the factors (cov_eta)
      # First, get biased sample covariance matrix per group (S)
      S <- S_unbiased[[g]] * (N_gs[[g]] - 1) / N_gs[[g]]
      cov_eta[[g]] <- M_mat[[g]] %*% (S - theta_g) %*% t(M_mat[[g]])
    }
  } else if (!is.list(step1model)) {
    # If not a list, then we only have one measurement block (all latent variables at the same time)
    if (!is.null(s1out)) {
      # If the user input their own step 1 results, use it
      S1output <- s1out
    } else if (is.null(s1out)) {
      S1output <- lavaan::cfa(
        model = step1model, data = centered, group = group,
        estimator = "ML", group.equal = constraints,
        se = "none", test = "none", baseline = FALSE, h1 = FALSE,
        implied = FALSE, loglik = FALSE,
        meanstructure = FALSE, group.partial = NonInv
      )
    }

    # Define some important objects
    # How many groups?
    ngroups <- lavInspect(S1output, "ngroups")
    N_gs <- lavInspect(S1output, "nobs") # nobs per group

    # all estimated model matrices, per group
    EST <- lavInspect(S1output, "est", add.class = FALSE, add.labels = TRUE)
    theta_gs <- lapply(EST, "[[", "theta")
    lambda_gs <- lapply(EST, "[[", "lambda")
    cov_eta <- lapply(EST, "[[", "psi") # cov_eta name refers to Variance of eta (eta being the latent variables)
  }

  # Which are our latent variables?
  lat_var <- lavNames(lavaanify(step1model, auto = TRUE), "lv")

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
    lines_model <- unlist(strsplit(unlist(step1model), "\n"))
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
    fake_measur <- cfa(model = rewritten, data = dat, do.fit = F)
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
  if (est_method == "local"){
    fake <- sem(
      model = step2model, sample.cov = rep(cov_eta[1], nclus),
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
    fake@Options$se     <- se
    fake@ParTable$start <- NULL
    fake@ParTable$est   <- NULL
    fake@ParTable$se    <- NULL
    fake@Options$start  <- "default"
    
  } else if (est_method == "global"){
    model.comb <- paste(step1model, step2model)
    fake <- sem(
      model = model.comb, data = dat, group = "group",
      do.fit = FALSE,
      baseline = FALSE,
      h1 = FALSE, check.post = FALSE,
      loglik = FALSE,
      sample.cov.rescale = FALSE,
      fixed.x = TRUE,
      meanstructure = meanstr
    )
    
    # Prepare a proper parameter table for later. We have to fix the MM parameters in the parameter table
    FakeprTbl <- parTable(fake)
    PT.MM     <- partable(S1output) # Parameter table MM
    idx.str   <- which(PT.MM$rhs %in% lat_var)
    PT.MM     <- PT.MM[-idx.str, ] # Remove structural parameters from PT.MM
    
    # Fill out MM parameters in the fake table
    # First make a full parameter column including group
    PT.MM$par     <- paste0(PT.MM$lhs, PT.MM$op, PT.MM$rhs, ".g", PT.MM$group)
    FakeprTbl$par <- paste0(FakeprTbl$lhs, FakeprTbl$op, FakeprTbl$rhs, ".g", FakeprTbl$group)
    
    # Introduce the correct MM parameters
    idx.par <- match(PT.MM$par, FakeprTbl$par) # indices of parameters from PT.MM on FakeprTbl
    idx.par <- idx.par[!is.na(idx.par)]
    FakeprTbl$ustart[idx.par] <- PT.MM$est[1:length(idx.par)] # The index makes sure we do not include est from constraints in the partable
    
    # Leave only structural parameters as free parameters
    FakeprTbl$free  <- 0
    FakeprTbl$free[is.na(FakeprTbl$ustart)] <- 1:sum(is.na(FakeprTbl$ustart)) 
    FakeprTbl$par   <- NULL
    FakeprTbl$est   <- NULL
    FakeprTbl$se    <- NULL
    FakeprTbl$start    <- NULL
    
    fake@Options$do.fit <- TRUE
    # fake@ParTable$start <- NULL
    # fake@ParTable$est   <- NULL
    # fake@ParTable$se    <- NULL
    fake@Options$start  <- "default"
  }
  
 
  # fake@Options$verbose <- TRUE

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
    if(est_method == "local"){
      fake_lv[[lv]] <- sem(
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
      fake_lv[[lv]]@Options$se     <- se
      fake_lv[[lv]]@ParTable$start <- NULL
      fake_lv[[lv]]@ParTable$est   <- NULL
      fake_lv[[lv]]@ParTable$se    <- NULL
      fake_lv[[lv]]@Options$start  <- "default"
      
    } else if (est_method == "global"){
      fake_lv[[lv]] <- sem(
        model = prTbl_lv[[lv]], data = centered, group = "group",
        do.fit = FALSE,
        baseline = FALSE,
        h1 = FALSE, check.post = FALSE,
        loglik = FALSE,
        sample.cov.rescale = FALSE,
        fixed.x = TRUE,
        meanstructure = meanstr
      )
      # Leave only structural parameters as free parameters
      prTbl_lv[[lv]]$free  <- 0
      prTbl_lv[[lv]]$free[is.na(prTbl_lv[[lv]]$ustart)] <- 1:sum(is.na(prTbl_lv[[lv]]$ustart)) 
      prTbl_lv[[lv]]$par   <- NULL
      prTbl_lv[[lv]]$est   <- NULL
      prTbl_lv[[lv]]$se    <- NULL
      prTbl_lv[[lv]]$start    <- NULL
      
      fake_lv[[lv]]@Options$do.fit <- TRUE
      # fake_lv[[lv]]@ParTable$start <- NULL
      # fake_lv[[lv]]@ParTable$est <- NULL
      # fake_lv[[lv]]@ParTable$se <- NULL
      fake_lv[[lv]]@Options$start <- "default"
    }
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
        if (est_method == "local"){
          s2out <- lavaan(
            slotOptions       = fake@Options,
            slotParTable      = fake@ParTable,
            sample.cov        = COV,
            sample.nobs       = rep(nrow(dat), nclus),
            # slotModel       = slotModel,
            # slotData        = fake@Data,
            # slotSampleStats = fake@SampleStats
          )
        } else if (est_method == "global"){
          s2out <- lavaan(
            slotOptions       = fake@Options,
            model             = FakeprTbl,
            sample.cov        = COV,
            sample.nobs       = rep(nrow(dat), nclus),
            meanstructure     = meanstr,
            se                = se,
            information       = "observed"
            # slotModel       = slotModel,
            # slotData        = fake@Data,
            # slotSampleStats = fake@SampleStats
          )
        }
        
      } else if (i > 1) {
        # After the first iteration
        # Run structural estimation once per endo LV
        s2out <- vector(mode = "list", length = length(endog))
        
        if (est_method == "local"){
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
        } else if (est_method == "global"){
          # browser()
          for (lv in 1:length(endog)) {
            s2out[[lv]] <- lavaan(
              slotOptions       = fake_lv[[lv]]@Options,
              model             = prTbl_lv[[lv]],
              sample.cov        = COV_lv[[lv]],
              sample.nobs       = rep(nrow(dat), nclus),
              meanstructure     = meanstr,
              se                = se,
              information       = "observed"
              # slotModel       = slotModel,
              # slotData        = fake@Data,
              # slotSampleStats = fake@SampleStats
            )
          }
        }
        
      }


      # start <- partable(s2out)$est

      # compute loglikelihood for all group/cluster combinations
      # Initialize matrices to store loglikelihoods
      loglik_gks <- matrix(data = 0, nrow = ngroups, ncol = nclus)
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
          EST_s2 <- lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
          beta_ks <- EST_s2[["beta"]]
          psi_ks <- EST_s2[["psi"]]
        } else if (nclus != 1) {
          EST_s2 <- lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
          beta_ks <- lapply(EST_s2, "[[", "beta") # Does not work with only one cluster
          psi_ks <- lapply(EST_s2, "[[", "psi")
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
        EST_s2_lv <- vector(mode = "list", length = length(endog))
        beta_ks_lv <- vector(mode = "list", length = length(endog))
        psi_ks_lv <- vector(mode = "list", length = length(endog))
        for (lv in 1:length(endog)) {
          if (nclus == 1) {
            EST_s2_lv[[lv]] <- lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
            beta_ks_lv[[lv]] <- EST_s2_lv[[lv]][["beta"]]
            psi_ks_lv[[lv]] <- EST_s2_lv[[lv]][["psi"]]
          } else if (nclus != 1) {
            EST_s2_lv[[lv]] <- lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
            beta_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "beta") # Does not work with only one cluster
            psi_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "psi")
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
          
          if (est_method == "local"){
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
          } else if (est_method == "global"){
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
  s2out <- results_nstarts[[best_idx]]
  LL <- loglik_nstarts[best_idx]
  z_gks <- z_gks_nstarts[[best_idx]]
  colnames(z_gks) <- paste("Cluster", seq_len(nclus))

  # Extract matrices from final step 2 output
  if (isFALSE(allG)) {
    if (nclus == 1) {
      EST_s2 <- lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE) # Estimated matrices step 2
      beta_ks <- EST_s2[["beta"]]
      psi_ks <- EST_s2[["psi"]]
    } else if (nclus != 1) {
      EST_s2 <- lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE) # Estimated matrices step 2
      beta_ks <- lapply(EST_s2, "[[", "beta")
      psi_ks <- lapply(EST_s2, "[[", "psi")
    }
  } else if (isTRUE(allG)) {
    EST_s2_lv <- vector(mode = "list", length = length(endog))
    beta_ks_lv <- vector(mode = "list", length = length(endog))
    psi_ks_lv <- vector(mode = "list", length = length(endog))
    for (lv in 1:length(endog)) {
      if (nclus == 1) {
        EST_s2_lv[[lv]] <- lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
        beta_ks_lv[[lv]] <- EST_s2_lv[[lv]][["beta"]]
        psi_ks_lv[[lv]] <- EST_s2_lv[[lv]][["psi"]]
      } else if (nclus != 1) {
        EST_s2_lv[[lv]] <- lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
        beta_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "beta") # Does not work with only one cluster
        psi_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "psi")
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

  # MODEL SELECTION
  # Get observed data log-likelihood using Kim's code (for model selection purposes)
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
  if (is.list(step1model)) {
    for (m in 1:m) {
      partbl <- parTable(S1output[[m]])
      free_load <- which(partbl$op == "=~" & is.na(partbl$ustart) & partbl$group == 1 & partbl$label != partbl$plabel)
      n_free <- n_free + length(free_load)
    }
  } else if (!is.list(step1model)) {
    partbl <- parTable(S1output)
    free_load <- which(partbl$op == "=~" & is.na(partbl$ustart) & partbl$group == 1 & partbl$label != partbl$plabel)
    n_free <- length(free_load)
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
  
  # Calculate AIC. 
  # Observed
  Obs.AIC <- (-2 * LL) + (nr_pars * 2)
  
  # Factors
  AIC <- (-2 * LL) + (nr_par_factors * 2)
  
  # Calculate the Standard Errors of the structural parameters
  if (isTRUE(do.se)){
    # To calculate the SE, we use the Hessian (second derivative) of the structural parameters.
    # For the derivative, we use the numDeriv package which requires a parameter vector and an obj. function
    # (1) Get all necesseary objects for the objective function
    # Get the parameter vector
    # Regressions (Cluster-specific!)
    beta_vec <- c() # Empty vector for the betas
    beta_nam <- c() # Empty vector for the names of such betas (not truly necessary, but makes everything easier to understand)
    
    # The loop identifies the non-zero parameters (i.e., the free parameters) in the beta matrices
    # It also saves the number of the parameter using the col and row names of the beta matrices (e.g., F3 ~ F1)
    for(k in 1:nclus){
      non_zer.idx <- which(unlist(beta_ks[[k]]) != 0)
      beta_vec <- c(beta_vec, unlist(beta_ks[[k]])[non_zer.idx])
      beta_nam <- c(beta_nam, as.vector(outer(X = rownames(beta_ks[[k]]), 
                                              Y = colnames(beta_ks[[k]]), 
                                              function(x, y) paste0(x, "~", y, ".", k)))[non_zer.idx])
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
        non_zer.idx <- which(unlist(psi_gks[[g, k]]) != 0)
        unique.idx  <- !duplicated(unlist(psi_gks[[g, k]])[non_zer.idx])
        cov_vec <- c(cov_vec, unlist(psi_gks[[g, k]])[non_zer.idx][unique.idx])
        cov_nam <- c(cov_nam, as.vector(outer(X = rownames(psi_gks[[g, k]]),
                                              Y = rownames(psi_gks[[g, k]]),
                                              function(x, y) paste0(x, "~~", y, ".", g, ".", k)))[non_zer.idx][unique.idx])
      }
    }
    
    cov_vec <- setNames(object = cov_vec, nm = cov_nam) # Name the vector (the first number after the . is the group, the second number is the cluster)
    
    # How many regressions and covariances?
    n_reg <- length(beta_vec)/nclus # nclus
    n_cov <- length(cov_vec)/(nclus*ngroups) # nclus*ngroups
    
    # Get indices needed for the objective function (used to input the value of the vector in corresponding place of the matrix)
    idx.beta <- which(beta_ks[[1]] %in% beta_vec[1:n_reg]) # indices of the free parameters in the beta matrix 
    idx.psi  <- which(psi_gks[[1, 1]] %in% cov_vec) # indices of the free parameters in the psi matrix 
    idx.psi_vec <- match(psi_gks[[1, 1]][idx.psi], cov_vec[1:n_cov]) # To repeat the free parameters that needs to be repeated (i.e, covariances)
    
    # I have group-specific covariances and group-cluster specific covariances.
    # Which ones are group-specific? 
    g.covs  <- unique(cov_vec[duplicated(cov_vec)]) # Group-specific (why duplicated? I extract all possible group-cluster combination, but the group-specific ones will be repeated on the second cluster) 
    g.covs  <- setNames(object = g.covs, nm = names(cov_vec[cov_vec %in% g.covs])[1:length(g.covs)])
    gk.covs <- cov_vec[!c(cov_vec %in% g.covs)] # Group-cluster specific
    
    # (2) Create the objective function of the step 2 parameters
    # The objective function simply takes the parameter vector and fills in the corresponding matrix to get the LL
    obj.S2 <- function(x, beta_mat, psi_mat, pi_ks, cov_eta, 
                       nclus, ngroups, N_gs, idx.beta, idx.psi, idx.psi_vec,
                       n_cov_exo){
      # browser()
      n_reg <- length(unlist(beta_mat)[unlist(beta_mat) != 0]) # How many free regression do we have?
      
      betas <- x[1:n_reg] # The first parameters are always the regressions. Extract them.
      covs  <- x[(n_reg + 1):length(x)] # The rest of the parameters are the covariances
      g.covs  <- covs[1:(n_cov_exo*ngroups)] # Which ones are the group-specific? Use the number of (co)variances * ngroups 
      gk.covs <- covs[!c(covs %in% g.covs)] # The rest are the group-cluster specific covariances
      
      # How many?
      n_g.psi <- length(g.covs) 
      n_gk.psi <- length(gk.covs)
      
      # Input the correct beta parameter in the matrix
      clus_idx <- n_reg/nclus # How many regression per cluster?
      for(k in 1:nclus){
        beta_vec <- betas[(((k - 1)*clus_idx) + 1):(clus_idx*k)]
        beta_mat[[k]][idx.beta] <- beta_vec
      }
      
      # Input the correct psi parameter in the matrix
      psi_gks <- psi_mat
      g.clus_idx <- n_g.psi/(ngroups) # How many g.psi per group?
      gk.clus_idx <- n_gk.psi/(nclus*ngroups) # How many gk.psi per group-cluster?
      gk <- 0 # group-cluster counter
      for(k in 1:nclus){
        for(g in 1:ngroups){
          gk <- gk + 1
          cov_vec <- c(g.covs[(((g - 1)*g.clus_idx) + 1):(g.clus_idx*g)], gk.covs[(((gk - 1)*gk.clus_idx) + 1):(gk.clus_idx*gk)])
          psi_mat[[g, k]][idx.psi] <- cov_vec[idx.psi_vec] # A part is group-specific and the other part is group-cluster-specific
        }
      }
      
      # Compute loglikelihood. Unrelated to the assignment of the parameters. This simply calculate the logL
      # compute loglikelihood for all group/cluster combinations
      # Initialize matrices to store loglikelihoods
      loglik_gks <- matrix(data = 0, nrow = ngroups, ncol = nclus)
      loglik_gksw <- matrix(data = 0, nrow = ngroups, ncol = nclus)
      
      Sigma <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
      I <- diag(nrow(cov_eta[[1]]))
      
      for(k in 1:nclus){
        for(g in 1:ngroups){
          # Estimate Sigma (factor covariance matrix of step 2)
          Sigma[[g, k]] <- solve(I - beta_mat[[k]]) %*% psi_mat[[g, k]] %*% t(solve(I - beta_mat[[k]]))
          # Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]])) # Force to be symmetric
          
          # Estimate the loglikelihood
          loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
            sample.mean = rep(0, nrow(cov_eta[[1]])),
            sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
            # sample.nobs = N_gks[g, k],
            sample.cov  = cov_eta[[g]], # Factor covariance matrix from step 1
            Mu          = rep(0, nrow(cov_eta[[1]])),
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
     
    # (3) Create a function to compute the Hessian numerically, and use it.
    compute_hessian <- function(f, x, d = 1e-5,
                                beta_mat, psi_mat, pi_ks, cov_eta, 
                                nclus, ngroups, N_gs, idx.beta, idx.psi, idx.psi_vec,
                                n_cov_exo) {
      # browser()
      n  <- length(x)
      H  <- matrix(0, n, n)
      colnames(H) <- rownames(H) <- names(x)
      f1 <- matrix(0, n)
      
      f0 <- f(x = x, 
              beta_mat = beta_mat, 
              psi_mat = psi_mat, 
              pi_ks = pi_ks, 
              cov_eta = cov_eta, 
              nclus = nclus, 
              ngroups = ngroups, 
              N_gs = N_gs, 
              idx.beta = idx.beta, 
              idx.psi = idx.psi, 
              idx.psi_vec = idx.psi_vec,
              n_cov_exo = n_cov_exo)
      
      for (i in 1:n) {
        x_up    <- x
        x_up[i] <- x_up[i] + d
        
        f1[i]   <- f(x = x_up, 
                     beta_mat = beta_mat, 
                     psi_mat = psi_mat, 
                     pi_ks = pi_ks, 
                     cov_eta = cov_eta, 
                     nclus = nclus, 
                     ngroups = ngroups, 
                     N_gs = N_gs, 
                     idx.beta = idx.beta, 
                     idx.psi = idx.psi, 
                     idx.psi_vec = idx.psi_vec,
                     n_cov_exo = n_cov_exo)
        
        for (j in 1:i) {
          x_up[j] <- x_up[j] + d
          
          f2      <- f(x = x_up, 
                       beta_mat = beta_mat, 
                       psi_mat = psi_mat, 
                       pi_ks = pi_ks, 
                       cov_eta = cov_eta, 
                       nclus = nclus, 
                       ngroups = ngroups, 
                       N_gs = N_gs, 
                       idx.beta = idx.beta, 
                       idx.psi = idx.psi, 
                       idx.psi_vec = idx.psi_vec,
                       n_cov_exo = n_cov_exo)
          
          H[i, j] <- (f2 - f1[i] - f1[j] + f0) / (d^2)
          H[j, i] <- H[i, j]
          
          x_up[j] <- x_up[j] - d
        }
      }
      
      return(H)
    }
    
    x <- c(beta_vec, g.covs, gk.covs)
    
    HESS <- compute_hessian(f         = obj.S2,
                            x         = x, 
                            d         = 1e-05, 
                            beta_mat  = beta_ks, 
                            psi_mat   = psi_gks, 
                            pi_ks     = colMeans(z_gks), 
                            cov_eta   = cov_eta,
                            nclus     = nclus,
                            ngroups   = ngroups,
                            n_cov_exo = n_cov_exo,
                            N_gs      = N_gs, 
                            idx.beta  = idx.beta, 
                            idx.psi   = idx.psi, 
                            idx.psi_vec = idx.psi_vec)
    
    # (4) Organize the SE for each parameter
    vector_SE <- setNames(diag(sqrt(ginv(-HESS, tol = 1e-05))), colnames(HESS)) # The SE comes from the inverse of the negative hessian
    
    SE.S2 <- function(x, beta_mat, psi_mat, cov_eta, 
                      nclus, ngroups, N_gs, 
                      idx.beta, idx.psi, idx.psi_vec, n_cov_exo){
      
      n_reg <- length(unlist(beta_mat)[unlist(beta_mat) != 0]) # How many free regression do we have?
      
      betas <- x[1:n_reg] # The first parameters are always the regressions. Extract them.
      covs  <- x[(n_reg + 1):length(x)] # The rest of the parameters are the covariances
      g.covs  <- covs[1:(n_cov_exo*ngroups)] # Which ones are the group-specific? Use the number of (co)variances * ngroups 
      gk.covs <- covs[!c(covs %in% g.covs)] # The rest are the group-cluster specific covariances
      
      # How many?
      n_g.psi <- length(g.covs) 
      n_gk.psi <- length(gk.covs)
      
      # Input the correct beta parameter in the matrix
      clus_idx <- n_reg/nclus # How many regression per cluster?
      for(k in 1:nclus){
        beta_vec <- betas[(((k - 1)*clus_idx) + 1):(clus_idx*k)]
        beta_mat[[k]][idx.beta] <- beta_vec
      }
      
      # Input the correct psi parameter in the matrix
      psi_gks <- psi_mat
      g.clus_idx <- n_g.psi/(ngroups) # How many g.psi per group?
      gk.clus_idx <- n_gk.psi/(nclus*ngroups) # How many gk.psi per group-cluster?
      gk <- 0 # group-cluster counter
      for(k in 1:nclus){
        for(g in 1:ngroups){
          gk <- gk + 1
          cov_vec <- c(g.covs[(((g - 1)*g.clus_idx) + 1):(g.clus_idx*g)], gk.covs[(((gk - 1)*gk.clus_idx) + 1):(gk.clus_idx*gk)])
          psi_mat[[g, k]][idx.psi] <- cov_vec[idx.psi_vec] # A part is group-specific and the other part is group-cluster-specific
        }
      }
      
      return(list(
        betas_SE  = beta_mat,
        psi_SE    = psi_mat,
        SE_vector = x
      ))
    }
    
    SE <- SE.S2(x         = vector_SE, 
                beta_mat  = beta_ks, 
                psi_mat   = psi_gks,
                cov_eta   = cov_eta,
                nclus     = nclus,
                ngroups   = ngroups,
                n_cov_exo = n_cov_exo,
                N_gs      = N_gs, 
                idx.beta  = idx.beta, 
                idx.psi   = idx.psi, 
                idx.psi_vec = idx.psi_vec)
  } else {
    SE <- NULL
  }
  
  # Re order matrices so that we get them in the following order:
  # (1) Exogenous latent variables
  # (2) Endogenous latent variables: independent and dependent variables at the same time
  # (3) Endogenous latent variables: only dependent variables

  # Reoder psi_ks and beta_ks by using the reorder function in the lapply function
  gro_clu <- ngroups * nclus
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
    z_gks         = z_gks,
    final_fit     = s2out,
    psi_gks       = psi_gks,
    lambda        = lambda_gs, # Lambda is invariant across all groups
    MM            = S1output, # Output of step 1 (measurement model)
    cov_eta       = cov_eta, # Factor covariance matrix from first step
    beta_ks       = beta_ks,
    loglikelihood = LL,
    loglik_gkw    = loglik_gksw,
    runs_loglik   = loglik_nstarts,
    obs_loglik    = Obs.LL,
    BIC = list(
      observed = list(BIC_N = Obs.BIC_N, BIC_G = Obs.BIC_G),
      Factors = list(BIC_N = BIC_N, BIC_G = BIC_G)
    ),
    AIC = list(
      observed = Obs.AIC,
      Factors = AIC
    ),
    NrPar         = list(Obs.nrpar = nr_pars, Fac.nrpar = nr_par_factors),
    N_gs          = N_gs,
    SE            = SE,
    est.vec       = list(beta = beta_vec, psi = cov_vec)
  ))
}
