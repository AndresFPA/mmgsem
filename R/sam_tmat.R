# Helper function to construct the T matrix (see Wall & Amemiya (2000) eq 7)
# The T-transformation can handle situations that the normal computations for the M-matrix cannot.
# E.g., when one value in the diagonal of theta is 0
# This code is partially taken from the lavaan Github repository (function: lav_sam_tmat)

sam_tmat <- function(lambda, theta){
  # Get necessary objects
  nvar <- nrow(lambda)
  nfac <- ncol(lambda)

  # Extract marker variable indices
  marker.idx <- lavaan:::lav_utils_get_marker(lambda)

  # Get C matrix
  C2 <- diag(nvar)
  C2[, marker.idx] <- -1 * lambda
  C <- C2[-marker.idx, , drop = FALSE] # Remove the rows corresponding to the marker variables

  # compute Sigma.ve and Sigma.vv
  Sigma.ve <- C %*% theta
  # Sigma.vv <- C %*% theta %*% t(C)
  Sigma.vv <- Sigma.ve %*% t(C)

  # construct 'Gamma' (and Gamma2) matrix
  # Gamma <- (t(Sigma.ve) %*% solve(Sigma.vv))[marker.idx,, drop = FALSE]
  Gamma <- try(t(solve(Sigma.vv, Sigma.ve)[, marker.idx, drop = FALSE]), # Normal inverse is not always possible
               silent = TRUE
  )
  if (inherits(Gamma, "try-error")) {
    tmp <- t(Sigma.ve) %*% MASS::ginv(Sigma.vv)
    Gamma <- tmp[marker.idx, , drop = FALSE]
  }
  Gamma2 <- matrix(0, nfac, nvar)
  Gamma2[, -marker.idx] <- Gamma
  Gamma2[, marker.idx] <- diag(nfac)

  # transformation matrix 'T' (we call it here 'Tmat')
  Tmat <- matrix(0, nvar, nvar)
  Tmat[-marker.idx, ] <- C
  Tmat[marker.idx, ] <- -Gamma2 %*% C2

  return(Tmat)
}
