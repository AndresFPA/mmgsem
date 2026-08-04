#' Internal helper: estimate the measurement model for MMGSEM
#'
#' Runs lavaan CFA for one or more measurement blocks and returns the
#' measurement model parameters, latent covariance estimates, and sample
#' covariance matrices needed for Step 2.
#'
#' @param S1 lavaan syntax string or list of syntax strings for the measurement model.
#' @param s1_fit optional pre-fitted lavaan object or list of objects.
#' @param centered data frame with centered observed variables.
#' @param group grouping variable name.
#' @param S_unbiased unused internal argument kept for compatibility.
#' @param ... additional arguments passed to lavaan::cfa.
#' @return A list with S1output, lambda_gs, theta_gs, cov_eta, ngroups, N_gs, and S_biased.
#' @keywords internal
Step1 <- function(S1, s1_fit = NULL, centered, group, S_unbiased = NULL, ...) {
  # Create a dummy lavaan model to extract biased sample covariances
  model_dummy <- if (is.list(S1)) unlist(S1) else S1

  s1_dummy <- lavaan::cfa(
    model = model_dummy,
    data = centered,
    group = group,
    test = "none",
    baseline = FALSE,
    loglik = FALSE,
    do.fit = FALSE,
    ...
  )

  S_biased <- lapply(lavaan::lavInspect(s1_dummy, "samplestats"), "[[", "cov")

  if (is.list(S1)) {
    M <- length(S1)

    if (!is.null(s1_fit)) {
      S1output <- s1_fit
    } else {
      S1output <- vector(mode = "list", length = M)
      for (m in seq_len(M)) {
        S1output[[m]] <- lavaan::cfa(
          model = S1[[m]],
          data = centered,
          group = group,
          test = "none",
          baseline = FALSE,
          loglik = FALSE,
          ...
        )
      }
    }

    ngroups <- lavaan::lavInspect(S1output[[1]], "ngroups")
    vars <- lavaan::lavNames(lavaan::lavaanify(S1, auto = TRUE))
    lat_var <- lavaan::lavNames(lavaan::lavaanify(S1, auto = TRUE), "lv")

    # Extract Lambda and Theta per group for each measurement block
    lambda_block <- vector(mode = "list", length = M)
    theta_block <- vector(mode = "list", length = M)
    for (m in seq_len(M)) {
      EST_block <- lavaan::lavInspect(S1output[[m]], "est")
      lambda_block[[m]] <- lapply(EST_block, "[[", "lambda")
      theta_block[[m]] <- lapply(EST_block, "[[", "theta")
    }

    lambda_group <- vector(mode = "list", length = ngroups)
    theta_group <- vector(mode = "list", length = ngroups)
    for (g in seq_len(ngroups)) {
      for (m in seq_len(M)) {
        lambda_group[[g]][[m]] <- lambda_block[[m]][[g]]
        theta_group[[g]][[m]] <- theta_block[[m]][[g]]
      }
      lambda_group[[g]] <- lavaan::lav_matrix_bdiag(lambda_group[[g]])
      theta_group[[g]] <- lavaan::lav_matrix_bdiag(theta_group[[g]])
      rownames(lambda_group[[g]]) <- vars
      colnames(lambda_group[[g]]) <- lat_var
      rownames(theta_group[[g]]) <- vars
      colnames(theta_group[[g]]) <- vars
    }

    lambda_gs <- lambda_group
    theta_gs <- theta_group
    N_gs <- lavaan::lavInspect(S1output[[1]], "nobs")

    cov_eta <- vector(mode = "list", length = ngroups)
    for (g in seq_len(ngroups)) {
      lambda_g <- lambda_gs[[g]]
      theta_g <- theta_gs[[g]]
      marker.idx <- lavaan:::lav_utils_get_marker(lambda_g)

      if (any(diag(theta_g) == 0)) {
        tmat <- sam_tmat(lambda = lambda_g, theta = theta_g)
        M_mat <- tmat[marker.idx, , drop = FALSE]
        rownames(M_mat) <- colnames(lambda_g)
        colnames(M_mat) <- colnames(theta_g)
      } else {
        M_mat <- solve(t(lambda_g) %*% solve(theta_g) %*% lambda_g) %*% t(lambda_g) %*% solve(theta_g)
      }

      cov_eta[[g]] <- M_mat %*% (S_biased[[g]] - theta_g) %*% t(M_mat)
    }
  } else {
    if (!is.null(s1_fit)) {
      S1output <- s1_fit
    } else {
      S1output <- lavaan::cfa(
        model = S1,
        data = centered,
        group = group,
        test = "none",
        baseline = FALSE,
        h1 = FALSE,
        implied = FALSE,
        loglik = FALSE,
        ...
      )
    }

    ngroups <- lavaan::lavInspect(S1output, "ngroups")
    N_gs <- lavaan::lavInspect(S1output, "nobs")
    EST <- lavaan::lavInspect(S1output, "est", add.class = FALSE, add.labels = TRUE)
    theta_gs <- lapply(EST, "[[", "theta")
    lambda_gs <- lapply(EST, "[[", "lambda")
    cov_eta <- lapply(EST, "[[", "psi")
  }

  return(list(
    S1output = S1output,
    lambda_gs = lambda_gs,
    theta_gs = theta_gs,
    cov_eta = cov_eta,
    ngroups = ngroups,
    N_gs = N_gs,
    S_biased = S_biased
  ))
}
