#' Mixture Multi-Group Structural Equation Modelling (MMGSEM)
#'
#' Performs a mixture clustering based on the structural parameters (i.e., regressions) of a SEM model.
#' The estimation is done in a step-wise fashion and uses an expectation-maximization (EM) algorithm in the second step.
#'
#' @param dat Data frame containing observed variables and a grouping variable..
#' @param S1 Step 1 (measurement model) specification using lavaan syntax. Can be a list of strings determing the number of measurement blocks (e.g., one string for the MM of factor 1, and a second string for the MM of factor 2).
#'           If s1_type = "mplus", S1 must be written in mplus syntax.
#' @param S2 Step 2 (structural model) specification using lavaan syntax.
#' @param group Name of the grouping variable (as a character string).
#' @param nclus Integer. Number of clusters to be estimated in Step 2.
#' @param seed Optional. Random seed for replicable results.
#' @param userStart Optional. A user-defined cluster membership matrix (dimensions: number of groups × number of clusters) with binary values (1 or 0; 1 indicates cluster membership). To be used when the user has prior knowledge about the data. There must be only one 1 for each row. Skips random starts.
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
#' @param s1_type String. Determines which model is used when estimating step 1. Can be "lavaan" (CFA), "blavaan" (BCFA), or "mplus" (ML-CFA).
#' @param s1_fit Optional. A fitted model object for Step 1. Can be from lavaan (CFA), blavaan (BCFA), or MplusAutomation (ML-CFA).
#' @param max_it Maximum number of iterations for Step 2 (default = 10000).
#' @param nstarts Number of random starts for Step 2 (default = 20)..
#' @param partition Initialisation partition for random starts: "hard" (default) or "soft".
#' @param endogenous_cov Logical. If TRUE (default), residual covariances among purely endogenous latent variables are estimated. If FALSE, they are fixed to 0 and only the residual variances are estimated.
#' @param endo_group_specific Logical. If TRUE (default), residual covariances are group- and cluster-specific (i.e., they are estimated for every group-cluster combination). If FALSE, they are cluster-specific (i.e., they are fixed to be equal across groups within a cluster). Note that, if FALSE, the residual covariances will also influence the clustering.
#' @param sam_method either "local" or "global. Follows local and global approaches from the SAM method. GLOBAL NOT FUNCTIONAL YET.
#' @param rescaling Only used when data is ordered. By default, MMGSEM uses the marker variable scaling approach. But identification
#'                  constraints with ordinal data (by default) are handled by standardizing the factors' variance in the first step.
#'                  The rescaling argument (either T or F) rescales the factor variances and loadings to the marker variable scaling
#'                  before running step 2. It is set to T by default (rescaling happens). If set to F, the factor variances are kept fixed to 1.
#' @param meanstr Logical. If TRUE, includes the mean structure in the model (e.g., for scalar invariance).
#' @param ordinal Logical. If TRUE, observed variables are treated as ordinal (default = FALSE)
#' @param ... MMGSEM relies on lavaan for the estimation of the first step (i.e., CFA). If needed, the users can pass any lavaan argument to MMGSEM
#'            and it will be considered when estimating the CFA. For instance, std.lv if users want standardized latent variables,
#'            group.equal for constraints, group.partial for non-invariances, etc.
#'
#' OUTPUT:
#' @return The function will return a list with the following results:
#' @return posteriors: A groups × clusters matrix of posterior membership probabilities.
#' @return modal_post: A hard classification matrix indicating the most likely cluster for each group.
#' @return final_fit: Lavaan fit of the best and final model (not useful, only for testing purposes).
#' @return MM: Fitted measurement model (Step 1; i.e., cfa, bcfa, or ml-cfa). Returns the user-supplied model if provided.
#' @return param: A list of model parameters, including lambda_gs, theta_gs, beta_ks, psi_gks, and cov_eta.
#' @return logLik: A list of log-likelihood values: final model (only Step 2), random starts, and full model (Step 1 + Step 2).
#' @return model_sel: Model selection metrics such as BIC, AIC, and ICL.
#' @return sample.stats: Observed sample covariance matrices.
#' @return NrPar: A list containing the number of parameters.
#' @return N_gs: The sample size per group.
#' @return nstarts: The number of random starts used in step 2.
#' @return ngroups: Total number of groups.
#' @return iter: The number of iterations needed to reach convergence (from the best random start).
#' @return R2: A matrix containing the explained variance of each endogenous latent variable per group.
#'
#' @export
mmgsem <- function(dat, S1 = NULL, S2 = NULL, s1_type = "lavaan",
                   group, nclus, seed = NULL, userStart = NULL, s1_fit = NULL,
                   max_it = 10000L, nstarts = 20L, printing = FALSE,
                   partition = "hard", endogenous_cov = TRUE,
                   endo_group_specific = TRUE,
                   sam_method = "local", meanstr = FALSE,
                   rescaling = F, only_slow = FALSE,
                   ...) {

  # Get arguments in ...
  # Such arguments are the ones that will pass on lavaan's functions
  dots_args <- list(...)
  constraints <- dots_args$group.equal
  noninv      <- dots_args$group.partial
  ordered     <- dots_args$ordered; if(is.null(ordered)){ordered <- F}
  std.lv      <- dots_args$std.lv;  if(is.null(std.lv)){std.lv <- F}
  missing     <- dots_args$missing; if(is.null(missing)){missing <- "listwise"}
  if(s1_type == "blavaan"){
    wiggle      <- dots_args$wiggle # approximate metric invariance
    wiggle.sd   <- dots_args$wiggle.sd # size of prior sd
    bcontrol    <- dots_args$bcontrol # parallelizing the chains
  }
  # Add a warning in case there is a pre-defined start and the user also requires a multi-start
  if (!(is.null(userStart)) && nstarts > 1) {
    warning("If a start is defined by the user, no multi-start is performed. The results correspond to the one start used an input")
    nstarts <- 1
  }

  # Get several values relevant for future steps
  g_name  <- as.character(unique(dat[, group]))
  if(s1_type != "mplus"){
    vars    <- lavaan::lavNames(lavaan::lavaanify(S1, auto = TRUE))
    lat_var <- lavaan::lavNames(lavaan::lavaanify(S1, auto = TRUE), "lv")
  }
  # n_var   <- length(vars)

  # Add an error in case of incompatibility in the arguments regarding the scale of the latent variables
  if(std.lv == T & rescaling == T){
    warning("std.lv = T and rescaling = T arguments set the factor variances to different scales. Please choose one scaling method.")
  }

  if(ordered == F & rescaling == T){
    stop("rescaling = T only works when ordered = T. When ordered = T, the scale of some of the factor variances are set to 1 (correlations). Rescaling = T effectively turns them back to covariances.")
  }

  # Change the syntax of the model in step 1 if the data is ordered
  s1ori <- NULL # Initialize s1ori object (necessary as input for Step2 function)

  if(ordered == T){
    # Save original syntax for later code
    s1ori <- S1

    # Get new syntax
    S1 <- as.character(
      semTools::measEq.syntax(configural.model = S1,
                              dat              = dat,
                              parameterization = "delta",
                              ordered          = vars,
                              ID.fac           = "std.lv",
                              ID.cat           = "Wu",
                              group            = group,
                              group.equal      = constraints,
                              group.partial    = noninv)
    )

    # When ordered = T, by default, measEq.syntax standardizes the lv following Wu&Estabrook(2016).
    # MMG-SEM does not work with standardized lv by default. Thus, a rescaling is needed
    rescaling <- T # Set to TRUE, it will come later in the code

    # It is possible to work with standardized lv by setting std.lv = T. This means that rescaling must be set to F
    if (std.lv == T){
      rescaling <- F
    }
  }

  # # Center the data per group (so that the mean for all variables in each group is 0)
  centered <- dat

  # if the mean structure is not required, then remove the mean structure of the data (i.e., center the data)
  # but, if the intercepts are required, then meanstr changes to TRUE
  if(isFALSE(meanstr) & "intercepts" %in% constraints){
    warning("If the intercepts are included in the constraints, then meanstr automatically changes to TRUE to include the mean structure.")
    meanstr <- T
  }

  # Only center if data is not categorical and the mean structure is not required
  if(s1_type != "mplus"){
    if(ordered == F){
      if(isFALSE(meanstr)){
        group.idx <- match(dat[, group], g_name)
        group.sizes <- tabulate(group.idx)
        group.means <- rowsum.default(as.matrix(dat[, vars]),
                                      group = group.idx, reorder = FALSE,
                                      na.rm = TRUE # For listwise deletion
        ) / group.sizes
        centered[, vars] <- dat[, vars] - group.means[group.idx, , drop = FALSE]
      }
    }

    if(missing == "fiml"){centered <- dat} # If we want to deal with the missing data using fiml, we cannot center the data
    # Centering the data requires the group means, which are dependent on possible NAs

    # Get sample covariance matrix per group (used later)
    # S_unbiased <- lapply(X = unique(centered[, group]), FUN = function(x) {
    #   cov(centered[centered[, group] == x, vars])
    # })
    }

  ## STEP 1 - MMG-SEM ----------------------------------------------------------------------------------------
  # Save the measurement model results
  # Call function to run Step 1 of MMG-SEM (estimates CFA)
  Step1_args <- list(S1         = S1,
                     s1_fit     = s1_fit,
                     centered   = centered,
                     group      = group)#,
  # S_unbiased = S_unbiased)
  if(s1_type == "lavaan"){
    Step1_args <- c(dots_args, Step1_args)
    MM <- do.call(what = Step1, args = Step1_args)
  } else if(s1_type == "mplus"){
    MM <- do.call(what = Step1_ML, args = Step1_args)
  } else if(s1_type =="blavaan"){
    Step1_args <- c(dots_args, Step1_args, seed = seed)
    MM <- do.call(what = Step1_BSEM, args = Step1_args)
  }

  # Extract necessary objects
  ngroups   <- MM$ngroups
  S1output  <- MM$S1output
  lambda_gs <- MM$lambda_gs
  theta_gs  <- MM$theta_gs
  cov_eta   <- MM$cov_eta
  N_gs      <- MM$N_gs
  S_biased  <- MM$S_biased

  gro_clu   <- ngroups * nclus

  if(s1_type != "lavaan"){ #mplus or blavaan: single-indicator approach - factor scores as observed variables
    vars    <- colnames(theta_gs[[1]])
    lat_var <- colnames(cov_eta[[1]])
    S_unbiased <- MM$S_unbiased
    dat <- MM$single_data #factor scores as observed data
  }

  # Only happens when ordinal = T
  # Rescale covariance matrices when ordinal
  if (rescaling == T){
    # browser()
    for (g in 1:ngroups) {
      # Extract the first loading of each item (the one that would be 1 if unstandardized)
      loadings <- apply(lambda_gs[[g]], 2, function(col) {col[which(col != 0)]}[1])
      # Multiply standardized variances with squared corresponding loading
      sds <- sqrt(diag(cov_eta[[g]]) * loadings^2)
      # Re-scale everything to correlations first (depending on the constraints, only the first group may have a correlation in cov_eta)
      cov_eta[[g]] <- stats::cov2cor(cov_eta[[g]])
      # Use lavaan's cor2cov to go back the covariances
      cov_eta[[g]] <- lavaan::cor2cov(R = cov_eta[[g]], sds = sds)
    }
  }

  # STEP 2 (EM algorithm for model estimation) -----------------------------------------------------
  SM <- Step2(ngroups             = ngroups,
              nclus               = nclus,
              nstarts             = nstarts,
              N_gs                = N_gs,
              seed                = seed,
              max_it              = max_it,
              cov_eta             = cov_eta,
              dat                 = dat,
              S2                  = S2,
              lat_var             = lat_var,
              ordered             = ordered,
              endo_group_specific = endo_group_specific,
              endogenous_cov      = endogenous_cov,
              lambda_gs           = lambda_gs,
              theta_gs            = theta_gs,
              S_unbiased          = S_unbiased,
              S1                  = S1,
              std.lv              = std.lv,
              partition           = partition,
              userStart           = userStart,
              printing            = printing,
              s1ori               = s1ori,
              only_slow           = only_slow,
              s1_type            = s1_type)

  iter            <- SM$iter    # Best start number of iterations
  z_gks           <- SM$z_gks   # Best start posteriors
  LL              <- SM$LL      # Best start loglikelihood
  s2out           <- SM$s2out   # Best start SM fit
  endog           <- SM$endog
  exog            <- SM$exog
  endog1          <- SM$endog1
  endog2          <- SM$endog2
  beta_ks         <- SM$beta_ks
  psi_gks         <- SM$psi_gks
  loglik_nstarts  <- SM$loglik_nstarts
  I               <- diag(length(lat_var))


  # MODEL SELECTION
  # Get observed data log-likelihood for model selection purposes)
  Sigma_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
  Obs.loglik_gks <- matrix(data = 0, nrow = ngroups, ncol = nclus)
  Obs.loglik_gksw <- matrix(data = 0, nrow = ngroups, ncol = nclus)
  pi_ks <- colMeans(z_gks)

  for (k in 1:nclus) {
    ifelse(test = (nclus == 1), yes = (beta <- beta_ks), no = (beta <- beta_ks[[k]]))
    for (g in 1:ngroups) {
      # S_biased <- S_unbiased[[g]] * (N_gs[[g]] - 1) / N_gs[[g]] # Deprecated, we already have S_biased
      var_eta <- solve(I - beta) %*% psi_gks[[g, k]] %*% t(solve(I - beta))
      Sigma_gks[[g, k]] <- lambda_gs[[g]] %*% var_eta %*% t(lambda_gs[[g]]) + theta_gs[[g]]
      # Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]]))
      # Sigma[[g, k]][lower.tri(Sigma[[g, k]])] <- t(Sigma[[g, k]])[lower.tri(Sigma[[g, k]])]
      Obs.loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
        sample.mean = rep(0, length(vars)),
        sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
        # sample.nobs = N_gks[g, k],
        sample.cov  = S_biased[[g]], # Item (observed) covariance matrix from step 1
        Mu          = rep(0, length(vars)),
        Sigma       = Sigma_gks[[g, k]] # Item (observed) model-implied covariance matrix including step 2
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
  if(s1_type == "lavaan"){
    n_free <- 0
    if (is.list(S1)) {
      M <- length(S1)
      for (m in 1:M) {
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
    if (endo_group_specific == F) { # Is endogenous covariance group-specific?
      nr_par_factors <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus)
      nr_pars <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus) + (n_res * ngroups) + (n_load - Q - n_free) + (n_free * ngroups)
    } else if (endo_group_specific == T) {
      nr_par_factors <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * ngroups) + (n_cov_endo2 * ngroups)
      nr_pars <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * ngroups) + (n_cov_endo2 * ngroups) + (n_res * ngroups) + (n_load - Q - n_free) + (n_free * ngroups)
    }} else {
      if (endo_group_specific == F) { # Is endogenous covariance group-specific?
        nr_par_factors <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus)
        nr_pars <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus) + n_load*ngroups + n_res*ngroups #number of loadings and residual variances
      } else if (endo_group_specific == T) {
        nr_par_factors <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * ngroups) + (n_cov_endo2 * ngroups)
        nr_pars <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * ngroups) + (n_cov_endo2 * ngroups) + n_load*ngroups + n_res*ngroups #number of loadings and residual variances
      }
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
    reorder(psi_gks[[x]], exog = exog, endog = endog)
  }), dim = c(ngroups, nclus))
  if (nclus == 1) {
    beta_ks <- reorder(beta_ks, exog = exog, endog = endog)
  } else if (nclus != 1) {
    beta_ks <- lapply(1:nclus, function(x) {
      reorder(beta_ks[[x]], exog = exog, endog = endog)
    }) # Does not work with only one cluster
  }

  names(beta_ks) <- paste("Cluster", seq_len(nclus))

  # Add the group name to the posterior matrix
  z_gks       <- as.data.frame(z_gks)
  z_gks       <- round(z_gks, 4)
  z_gks$Group <- g_name
  g_col_idx   <- ncol(z_gks)
  z_gks       <- z_gks[,c(g_col_idx, 1:(g_col_idx-1))] # Reorder with Group column as the first column

  # Last step, compute the R2 of each endogenous variable per group
  # Compute the modal clustering to extract the correct psi_gks
  post_no_group <- as.matrix(z_gks[, 2:ncol(z_gks)])
  if(nclus == 1){
    Modal_posteriors <- post_no_group
  } else {
    Modal_posteriors <- t(apply(post_no_group, 1, function(x) as.numeric(x == max(x))))
  }
  Modal_posteriors <- as.data.frame(Modal_posteriors)
  # browser()
  # Extract relevant psi_gks
  psi_idx <- which(x = Modal_posteriors == 1, arr.ind = T) # array indices of the relevant group-cluster combinations

  # Reorder based on group number insted of column number
  correct_idx <- order(psi_idx[,1])
  psi_idx <- psi_idx[correct_idx,]

  # Get total and residual variance objects (in array form)
  res_var <- array(data = NA, dim = c(Q, Q, ngroups))
  tot_var <- array(unlist(cov_eta), dim = c(Q, Q, ngroups))
  for(g in 1:ngroups){
    res_var[, , g] <- psi_gks[[psi_idx[g, 1], psi_idx[g, 2]]]
  }

  colnames(res_var) <- rownames(res_var) <- lat_var
  colnames(tot_var) <- rownames(tot_var) <- lat_var

  # Compute R2
  R2 <- array(data = NA, dim = c(Q, Q, ngroups))
  R2 <- 1 - (res_var/tot_var)
  R2 <- apply(R2, 3, diag) # Get only the explained variances
  R2 <- R2[endog, ]
  if(length(endog) == 1){R2 <- t(R2)} # Ensure we have a matrix instead of a vector, even if we only have one endog factor
  colnames(R2) <- g_name

  output <- (list(
    posteriors    = z_gks, # posterior probabilities
    modal_post    = Modal_posteriors, # hard cluster memberships
    final_fit     = s2out, # Final fit of step 2 (contains all group-cluster combinations)
    MM            = S1output, # Output of step 1 (measurement model)
    param         = list(psi_gks = psi_gks, lambda = lambda_gs, # Lambda is invariant across all groups
                         theta = theta_gs, beta_ks = beta_ks, cov_eta = cov_eta), # Factor covariance matrix from first step
    logLik        = list(loglik        = LL, # Final logLik of the model (its meaning depends on argument "sam_method")
                         global_loglik = ifelse(test = sam_method == "global", global_LL, NA), # Only valid if sam_method = "global"
                         #     loglik_gksw   = loglik_gksw, # Weighted logLik per group-cluster combinations
                         runs_loglik   = loglik_nstarts, # loglik for each start
                         obs_loglik    = Obs.LL), # Only useful if fit = "local"
    model_sel     = list(BIC        = list(observed = list(BIC_N = Obs.BIC_N, BIC_G = Obs.BIC_G),
                                           Factors = list(BIC_N = BIC_N, BIC_G = BIC_G)),
                         AIC        = list(observed = Obs.AIC, Factors = AIC),
                         AIC3       = list(observed = Obs.AIC3, Factors = AIC3),
                         R2_entropy = R2_entropy,
                         ICL        = list(observed = Obs.ICL, Factors = ICL)),
    sample.stats  = list(S = S_biased, n_cov_exo = n_cov_exo),
    NrPar         = list(Obs.nrpar = nr_pars, Fac.nrpar = nr_par_factors),
    N_gs          = N_gs,
    nstarts       = nstarts,
    ngroups       = ngroups,
    sam_method    = sam_method,
    iterations    = iter,
    R2            = R2
  ))

  class(output) <- "mmgsem"

  return(output)
}


