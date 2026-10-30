# -------------------------------------------------------------------------
# M-STEP: Update HMM Parameters
# -------------------------------------------------------------------------

MStep_HMM_Transitions_long <- function(gamma, xi) {

  # Dimensions of xi: [Groups, Time-1, From, To]
  dim_xi <- dim(xi)
  K <- dim_xi[3] # Number of clusters

  # --- 1. Update Initial Probabilities (Pi) ---
  # Average of Gamma at Time 1 across all groups
  # gamma[, 1, ] is [Groups x Clusters]
  gamma_first <- matrix(gamma[, 1, , drop = FALSE],
                        nrow = dim(gamma)[1],
                        ncol = K)
  pi_new <- colMeans(gamma_first)

  # Ensure it sums to 1 (handle small numerical noise)
  pi_new <- pi_new / sum(pi_new)

  # --- 2. Update Transition Matrix (A) ---
  # We sum Xi over Groups AND Time
  # "How many times did ANY group transition from j to k at ANY time?"

  # Initialize new matrix
  A_new <- matrix(0, nrow=K, ncol=K)

  # Iterate over From (j) and To (k)
  for (j in 1:K) {
    for (k in 1:K) {
      # Extract all transitions j->k (all groups, all times)
      transitions_jk <- xi[, , j, k]

      # Sum them up
      A_new[j, k] <- sum(transitions_jk)
    }
  }

  # Normalize rows to sum to 1
  # If a row sum is 0 (state never visited), keep it uniform or identity
  row_sums <- rowSums(A_new)
  for(j in 1:K) {
    if(row_sums[j] > 0) {
      A_new[j, ] <- A_new[j, ] / row_sums[j]
    } else {
      # Fallback if state is dead
      A_new[j, ] <- 1/K
    }
  }

  # Return the most important objects
  return(list(pi = pi_new, A = A_new))
}
