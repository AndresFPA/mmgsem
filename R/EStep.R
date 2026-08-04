#' Internal helper: compute the E-step posterior membership probabilities.
#'
#' Computes the posterior cluster membership probabilities from the current
#' log-likelihood values and cluster weights.
#'
#' @param pi_ks numeric vector of cluster prior probabilities.
#' @param ngroup number of groups.
#' @param nclus number of clusters.
#' @param loglik matrix of log-likelihood values for each group-cluster combination.
#' @return Matrix of posterior membership probabilities.
#' @keywords internal
EStep <- function(pi_ks, ngroup, nclus, loglik) {
  z_gks <- matrix(NA_real_, nrow = ngroup, ncol = nclus)
  for (g in seq_len(ngroup)) {
    z_gks[g, ] <- log(pi_ks) + loglik[g, ]
    z_gks[g, ] <- exp(z_gks[g, ] - max(z_gks[g, ]))
  }
  z_gks / rowSums(z_gks)
}
