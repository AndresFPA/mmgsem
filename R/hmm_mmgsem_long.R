hmm_mmgsem_long <- function(cov_eta, S2, ngroups, nstates, ntimes, N_gs, max_it, s2_fast, nstarts = 20,
                    printing = F){

  # Set up the random starts
  # Initialize all relevant objects
  logliks <- vector(mode = "list", length = nstarts)
  params  <- vector(mode = "list", length = nstarts)
  As      <- vector(mode = "list", length = nstarts)
  pis     <- vector(mode = "list", length = nstarts)
  posts   <- vector(mode = "list", length = nstarts)
  iterations <- integer(nstarts)
  converged <- logical(nstarts)

  for(s in 1:nstarts){
    if(isTRUE(printing)){print(paste("Random Start #", s, "------------------------------------------------------"))}
    # Get a random start for the posteriors per time point (marginals)
    post <- array(data = 0, dim = c(ngroups, nstates, ntimes))
    for(t in 1:ntimes){
      post[, , t] <- matrix(
        data = stats::runif(n = nstates * ngroups),
        ncol = nstates,
        nrow = ngroups
      )
      post[, , t] <- post[, , t] / rowSums(post[, , t])
    }

    # Swap dimensions for easier manipulation later on
    gamma <- aperm(a = post, perm = c(1,3,2))

    # Get random start for the transition probabilities (joint probabilities)
    ntrans <- ntimes - 1
    xi <- array(data = 0, dim = c(ngroups, ntrans, nstates, nstates))

    # Per group and transition, generate random probabilities of going from cluster k to cluster s
    # The sum of joint probabilities must be 1 (per group and transition)
    for(g in 1:ngroups){
      for(t in 1:ntrans){
        xi[g, t, , ] <- matrix(
          data = stats::runif(n = nstates * nstates),
          ncol = nstates,
          nrow = nstates
        )
        xi[g, t, , ] <- xi[g, t, , ]/sum(xi[g, t, , ]) # Ensure sum results in 1
      }
    }

    # ----------------------------------------------------------------------------
    # Start while loop for convergence
    # ----------------------------------------------------------------------------

    # Start convergence loop
    i              <- 0 # Initialize iterations
    prev_LL        <- 0 # previous loglikelihood initialization
    diff_LL        <- 1 # Set a diff of 1 just to start the while loop
    sem_parameters <- NULL

    while(abs(diff_LL) > 1e-6 && i < max_it){
      i <- i + 1
      # --------------------------------------------------------------------------------------------------------------------
      # --------------------------------------------------------------------------------------------------------------------
      # EM ALGORITHM (Backward-Forward version)
      # Why start with the M-Step?
      # Our random start is the probabilities, not the parameters.
      # Thus, our first step is to use those posteriors to *M*aximize the parameter estimates
      # --------------------------------------------------------------------------------------------------------------------
      # --------------------------------------------------------------------------------------------------------------------

      # ---------------------------------------------
      # Do M-Step.
      # ---------------------------------------------
      # 1. Maximize the posterior and transition probabilities
      # What is necessary?
      # Gamma: Marginal probabilities (posteriors) -> Probability of group g belonging to cluster k at timepoint t
      ### Dimensions -> [ngroups, ntimes, nstates]
      # Xi: Joint probabilities -> Probability of group g moving from cluster k to cluster j at transition tr
      ### Dimensions -> [Groups, Time-1, nstates(from), nstates(to)]
      M_results <- MStep_HMM_Transitions_long(gamma = gamma, xi = xi)
      # browser()
      # 2. Based on the resulting probabilities, maximize the SEM parameters
      if(isTRUE(s2_fast)){
        sem_parameters <- sem_est_fast_long(S2      = S2,
                                       post    = gamma,
                                       cov_eta = cov_eta,
                                       nstates = nstates,
                                       ntimes  = ntimes,
                                       ngroups = ngroups,
                                       N_gs    = N_gs,
                                       ite     = i,
                                       sem_par = sem_parameters)
      } else {
        sem_parameters <- sem_est_slow_long(S2      = S2,
                                       post    = gamma,
                                       cov_eta = cov_eta,
                                       nstates = nstates,
                                       ntimes  = ntimes,
                                       ngroups = ngroups,
                                       N_gs    = N_gs)
      }


      # ---------------------------------------------
      # Do E-Step.
      # ---------------------------------------------
      # What is necessary?
      # Inputs:
      # log_liks:  List of length T. Each element is a matrix [Groups x Clusters]
      #            representing log(P(Data | Cluster)) at that time point.
      # log_pi:    Vector of length K. log(Initial Probabilities).
      # log_A:     Matrix [K x K]. log(Transition Probabilities).
      E_results <- EStep_HMM_General_long(log_liks = sem_parameters$log_liks,
                                     log_pi   = log(M_results$pi),
                                     log_A    = log(M_results$A))

      gamma   <- E_results$gamma
      xi      <- E_results$xi
      LL      <- E_results$LL

      # Check convergence
      diff_LL <- LL - prev_LL # Check difference between previous and current loglikelihood
      prev_LL <- LL           # Update previous loglikelihood
      if (i == 1L) diff_LL <- Inf
      if(isTRUE(printing)){print(paste("Iteration", i, "; LL =", LL, "; Diff =", diff_LL))}
    }

    posts[[s]]   <- aperm(gamma, c(1,3,2))
    As[[s]]      <- M_results$A
    pis[[s]]     <- M_results$pi
    params[[s]]  <- sem_parameters
    logliks[[s]] <- LL
    iterations[s] <- i
    converged[s] <- abs(diff_LL) <= 1e-6
  }

  loglik_values <- unlist(logliks, use.names = FALSE)
  if (!any(is.finite(loglik_values))) {
    stop("All longitudinal random starts produced non-finite likelihoods.",
         call. = FALSE)
  }
  loglik_values[!is.finite(loglik_values)] <- -Inf
  best_start <- which.max(loglik_values)
  post       <- posts[[best_start]]
  A          <- As[[best_start]]
  pi         <- pis[[best_start]]
  param      <- params[[best_start]]
  loglik     <- logliks[[best_start]]
  s2out      <- param$s2out

  return(list(
    post        = post,
    A           = A,
    pi          = pi,
    logLik      = loglik,
    iterations  = iterations[best_start],
    converged   = converged[best_start],
    s2out       = s2out,
    param       = param,
    prev_starts = list(
      logliks = logliks,
      As      = As,
      pis     = pis,
      posts   = posts,
      params  = params,
      iterations = iterations,
      converged = converged
    )
  ))
}

