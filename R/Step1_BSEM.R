#' Internal helper: estimate Step 1 using Bayesian CFA with blavaan.
#'
#' Runs a blockwise BCFA for the measurement model, extracts latent score
#' information, and computes the sample and factor covariance matrices required
#' for Step 2 estimation.
#'
#' @param S1 list of lavaan syntax character strings defining measurement blocks.
#' @param s1_fit optional list of pre-fitted blavaan models.
#' @param centered data frame with centered observed variables.
#' @param group grouping variable name.
#' @param group.equal invariance constraints passed to blavaan::bcfa.
#' @param wiggle approximate invariance targets.
#' @param wiggle.sd prior standard deviation for wiggle constraints.
#' @param bcontrol parallelization settings for blavaan.
#' @param seed random seed for blavaan estimation.
#' @param ... additional arguments forwarded to blavaan::bcfa.
#' @return A list with S1output, lambda_gs, theta_gs, cov_eta, ngroups, N_gs, single_data, S_biased, and S_unbiased.
#' @keywords internal
Step1_BSEM <- function(S1, s1_fit = NULL, centered, group,
                       group.equal = c("loadings"),
                       wiggle = c("loadings"),
                       wiggle.sd = sqrt(0.1),
                       bcontrol = list(cores = 3),
                       seed = 100, ...) {
  start_time_step1 <- Sys.time()
  lat_var <- lavNames(lavaanify(S1, auto = TRUE), "lv")

  M <- length(S1)
  nfactors <- length(lat_var)

  if (is.null(s1_fit)) {
    s1_fit <- vector(mode = "list", length = M)
  }

  extract_results <- function(fit) {
    list(
      est = blavInspect(fit, "est", add.class = FALSE, add.labels = TRUE),
      lvs = blavInspect(fit, "lvs"),
      lvmeans = blavInspect(fit, "lvmeans"),
      fit = fitMeasures(fit),
      psrf = blavInspect(fit, "psrf")
    )
  }

  bfit_MM <- EST <- lvs <- lvmeans <- fit_blavaan <- psrf_list <- vector(mode = "list", length = M)

  for (i in seq_len(M)) {
    if (is.null(s1_fit[[i]])) {
      s1_fit[[i]] <- blavaan::bcfa(
        data = centered,
        model = S1[[i]],
        group = group,
        group.equal = group.equal,
        wiggle = wiggle,
        wiggle.sd = wiggle.sd,
        save.lvs = TRUE,
        bcontrol = bcontrol,
        seed = seed
      )
    }

    bfit_MM[[i]] <- s1_fit[[i]]
    res <- extract_results(bfit_MM[[i]])
    EST[[i]] <- res$est
    lvs[[i]] <- res$lvs
    lvmeans[[i]] <- res$lvmeans
    fit_blavaan[[i]] <- res$fit
    psrf_list[[i]] <- res$psrf
  }

  lvmeans <- do.call(cbind, lvmeans)
  sample_size <- nrow(centered)
  g_name <- as.character(unique(centered[, group]))
  group.idx <- match(centered[, group], g_name)
  group.sizes <- tabulate(group.idx)
  N_gs <- group.sizes
  ngroups <- length(group.sizes)

  postmeanfs <- data.frame(cbind(lvmeans, centered[, group]))
  vars_single <- paste0("single", seq_len(nfactors))
  colnames(postmeanfs) <- c(vars_single, "group")
  postmeanfs_list <- split(postmeanfs[, vars_single], postmeanfs$group)

  fcovdat <- array(unlist(lapply(postmeanfs_list, cov)), dim = c(nfactors, nfactors, ngroups))

  postvarfs <- do.call(cbind, lapply(lvs, function(lv) {
    lvs_factor <- do.call("rbind", lv)
    Mfactors <- ncol(lvs_factor) / sample_size
    factor_vars <- lapply(seq_len(Mfactors), function(f) {
      factor_cols <- ((f - 1) * sample_size + 1):(f * sample_size)
      apply(lvs_factor[, factor_cols], 2, var)
    })
    do.call(cbind, factor_vars)
  }))
  postvarfs <- data.frame(cbind(postvarfs, centered[, group]))
  colnames(postvarfs) <- c(vars_single, "group")
  postvarfs_list <- split(postvarfs[, vars_single], postvarfs$group)

  fvar_compute <- replicate(nfactors, numeric(ngroups), simplify = FALSE)
  for (f in seq_len(nfactors)) {
    for (g in seq_len(ngroups)) {
      fvar_compute[[f]][g] <- var(postmeanfs_list[[g]][, f]) + mean(postvarfs_list[[g]][, f])
    }
  }

  rho <- lapply(seq_len(nfactors), function(i) {
    sapply(seq_len(ngroups), function(j) var(postmeanfs_list[[j]][, i]) / fvar_compute[[i]][j])
  })
  lambda <- array(0, dim = c(nfactors, nfactors, ngroups))
  for (j in seq_len(ngroups)) {
    lambda[, , j] <- diag(sapply(rho, `[`, j))
  }
  dimnames(lambda) <- list(vars_single, lat_var, NULL)

  resvar <- lapply(seq_len(ngroups), function(j) {
    sapply(seq_len(nfactors), function(i) fvar_compute[[i]][j] * rho[[i]][j] * (1 - rho[[i]][j]))
  })
  theta <- array(unlist(lapply(resvar, function(x) diag(x))), dim = c(nfactors, nfactors, ngroups))
  dimnames(theta) <- list(vars_single, vars_single, NULL)

  cov_eta <- array(0, dim = c(nfactors, nfactors, ngroups))
  for (j in seq_len(ngroups)) {
    cov_eta[, , j] <- solve(lambda[, , j]) %*% (fcovdat[, , j] - theta[, , j]) %*% solve(t(lambda[, , j]))
  }
  dimnames(cov_eta) <- list(lat_var, lat_var, NULL)

  centered_fs <- postmeanfs
  g_name <- as.character(unique(postmeanfs$group))
  group.idx <- match(postmeanfs$group, g_name)
  group.sizes <- tabulate(group.idx)
  group.means <- rowsum.default(as.matrix(postmeanfs[, vars_single]), group = group.idx, reorder = FALSE, na.rm = FALSE) / group.sizes
  centered_fs[, vars_single] <- postmeanfs[, vars_single] - group.means[group.idx, , drop = FALSE]

  S_unbiased <- lapply(unique(centered_fs$group), function(x) {
    cov(centered_fs[centered_fs$group == x, vars_single])
  })
  S_biased <- lapply(seq_len(ngroups), function(g) {
    S_unbiased[[g]] * (N_gs[g] - 1) / N_gs[g]
  })

  lambda_gs <- lapply(seq_len(dim(lambda)[3]), function(x) lambda[, , x])
  theta_gs <- lapply(seq_len(dim(theta)[3]), function(x) theta[, , x])
  cov_eta <- lapply(seq_len(dim(cov_eta)[3]), function(i) cov_eta[, , i])

  end_time_step1 <- Sys.time()
  step1_duration <- difftime(end_time_step1, start_time_step1, units = "mins")
  cat("Step 1 completed in", step1_duration, "minutes\n")

  return(list(
    S1output = s1_fit,
    lambda_gs = lambda_gs,
    theta_gs = theta_gs,
    cov_eta = cov_eta,
    ngroups = ngroups,
    N_gs = N_gs,
    single_data = centered_fs,
    S_biased = S_biased,
    S_unbiased = S_unbiased
  ))
}
