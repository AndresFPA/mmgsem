# -------------------------------------------------------------------------
# HELPER: Log-Sum-Exp Function
# -------------------------------------------------------------------------
# We need this to add probabilities in log-space: log(A + B)
# formula: log(sum(exp(x))) = max(x) + log(sum(exp(x - max(x))))
# This prevents computer "underflow" (rounding small numbers to zero).

log_sum_exp_long <- function(x) {
  # Find the maximum value in the vector x
  max_val <- max(x)

  # If max is -Inf (all probs are 0), return -Inf
  if (!is.finite(max_val)) return(-Inf)

  # Perform the shift to ensure exponentiation is stable
  # Then take the log of the sum, and add the max back
  return(max_val + log(sum(exp(x - max_val))))
}

# -------------------------------------------------------------------------
# HELPER: Reordering functions
# -------------------------------------------------------------------------
# Create a function to reorder_long matrices (used later). Done due to:
# (1) It allows us to organize the matrices in an easier to understand order. Exogenous variables first, and endogeonus later.
# (2) Used to make sure we are multiplying matrices in the correct way

# Re-order for factors
reorder_long <- function(x, exog = exog, endog = endog) {
  x <- x[c(exog, endog), c(exog, endog)]
  return(x)
}
