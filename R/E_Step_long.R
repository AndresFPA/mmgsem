# -------------------------------------------------------------------------
# GENERALIZED E-STEP (Forward-Backward Algorithm)
# -------------------------------------------------------------------------
# Inputs:
#   log_liks:  List of length T. Each element is a matrix [Groups x Clusters]
#              representing log(P(Data | Cluster)) at that time point.
#   log_pi:    Vector of length K. log(Initial Probabilities).
#   log_A:     Matrix [K x K]. log(Transition Probabilities).
# -------------------------------------------------------------------------

EStep_HMM_General_long <- function(log_liks, log_pi, log_A) {

  # --- 1. SETUP DIMENSIONS ---
  T_points <- length(log_liks)            # Total number of time points
  N_groups <- nrow(log_liks[[1]])         # Number of groups (countries)
  K_clus   <- ncol(log_liks[[1]])         # Number of clusters

  # Initialize Alpha (Forward) matrix: [Groups x Time x Clusters]
  # Stores: log P(Observation_1...t, State_t = k)
  log_alpha <- array(-Inf, dim = c(N_groups, T_points, K_clus))

  # Initialize Beta (Backward) matrix: [Groups x Time x Clusters]
  # Stores: log P(Observation_t+1...T | State_t = k)
  log_beta  <- array(-Inf, dim = c(N_groups, T_points, K_clus))

  # -----------------------------------------------------------------------
  # 2. FORWARD PASS (The "Filtering" Step)
  # -----------------------------------------------------------------------

  # Initialization at t = 1
  # alpha_1(k) = pi(k) * P(Data_1 | k)
  for (k in 1:K_clus) {
    # Add in log-space is multiplication in normal space
    log_alpha[, 1, k] <- log_pi[k] + log_liks[[1]][, k]
  }

  # Recursion for t = 2 to T
  for (t in 2:T_points) {
    for (k_curr in 1:K_clus) { # Current state (at t)

      # We need to sum over all possible previous states (k_prev)
      # Calculation: sum( alpha[t-1, k_prev] * A[k_prev, k_curr] )

      # Create a temporary matrix to hold values before summing
      # Rows = Groups, Cols = Previous States
      tmp_matrix <- matrix(-Inf, nrow = N_groups, ncol = K_clus)

      for (k_prev in 1:K_clus) {
        # Add alpha from yesterday + transition probability to today
        tmp_matrix[, k_prev] <- log_alpha[, t-1, k_prev] + log_A[k_prev, k_curr]
      }

      # Apply log_sum_exp_long row-wise (summing over previous states) to get single value per group
      summed_prev <- apply(tmp_matrix, 1, log_sum_exp_long)

      # Multiply by likelihood of data today (add in log space)
      log_alpha[, t, k_curr] <- summed_prev + log_liks[[t]][, k_curr]
    }
  }

  # -----------------------------------------------------------------------
  # 3. BACKWARD PASS (The "Smoothing" Step)
  # -----------------------------------------------------------------------

  # Initialization at t = T
  # Beta is 1.0 (or 0 in log space) by definition at the end
  log_beta[, T_points, ] <- 0

  # Recursion for t = T-1 down to 1
  for (t in (T_points-1):1) {
    for (k_curr in 1:K_clus) { # Current state (at t)

      # We need to sum over all possible NEXT states (k_next)
      # Calculation: sum( A[k_curr, k_next] * P(Data_t+1 | k_next) * beta[t+1, k_next] )

      # Create temp matrix: Rows = Groups, Cols = Next States
      tmp_matrix <- matrix(-Inf, nrow = N_groups, ncol = K_clus)

      for (k_next in 1:K_clus) {
        # Transition to next + Likelihood of next data + Beta of next future
        term <- log_A[k_curr, k_next] +           # Transition cost
          log_liks[[t+1]][, k_next] +       # Emission cost (tomorrow)
          log_beta[, t+1, k_next]           # Future cost

        tmp_matrix[, k_next] <- term
      }

      # Sum over all possible next states
      log_beta[, t, k_curr] <- apply(tmp_matrix, 1, log_sum_exp_long)
    }
  }

  # -----------------------------------------------------------------------
  # 4. COMPUTE POSTERIORS (Gamma and Xi)
  # -----------------------------------------------------------------------

  # A. Marginal Probability (Gamma): P(State_t = k | All Data)
  # log_gamma = log_alpha + log_beta - log_total_likelihood

  log_gamma <- array(-Inf, dim = c(N_groups, T_points, K_clus))

  # Calculate Total Likelihood per group (Normalization Constant)
  # We can just sum alpha at time T over all clusters
  log_lik_total <- apply(log_alpha[, T_points, , drop = F], 1, log_sum_exp_long)

  for (t in 1:T_points) {
    for (k in 1:K_clus) {
      # Numerator: alpha * beta
      numer <- log_alpha[, t, k] + log_beta[, t, k]
      # Divide by total likelihood (subtract in log space)
      log_gamma[, t, k] <- numer - log_lik_total
    }
  }

  # B. Joint Transition Probability (Xi): P(State_t = j, State_t+1 = k | Data)
  # Needed to update the Transition Matrix A
  # Dimensions: [Groups x Time-1 x From_Cluster x To_Cluster]

  log_xi <- array(-Inf, dim = c(N_groups, T_points-1, K_clus, K_clus))

  for (t in 1:(T_points-1)) {
    for (j in 1:K_clus) {     # From State j
      for (k in 1:K_clus) {   # To State k

        # Numerator formula:
        # alpha[t, j] + A[j,k] + Likelihood[t+1, k] + beta[t+1, k]

        numer <- log_alpha[, t, j] +    # History up to t
          log_A[j, k] +                 # Jump j -> k
          log_liks[[t+1]][, k] +        # Evidence at t+1
          log_beta[, t+1, k]            # Future after t+1

        # Normalize
        log_xi[, t, j, k] <- numer - log_lik_total
      }
    }
  }

  # Convert back to normal probability space (exp)
  gamma <- exp(log_gamma)
  xi    <- exp(log_xi)

  return(
    list(gamma = gamma,
              xi = xi,
              LL = sum(log_lik_total) # Use the log_lik_total (alpha variable at the last time point) as the loglikelihood value
         )
    )
}
