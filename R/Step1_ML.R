#' Internal helper: estimate Step 1 using Mplus ML-CFA.
#'
#' Runs Mplus ML-CFA for each measurement block, extracts factor-score-based
#' single-indicator representations, and computes group-specific factor
#' covariance matrices.
#'
#' @param S1 list containing Mplus model components and latent variable names.
#' @param s1_fit optional list of pre-computed Mplus results.
#' @param centered data frame with centered observed variables.
#' @param group grouping variable name.
#' @return A list with S1output, lambda_gs, theta_gs, cov_eta, ngroups, N_gs, single_data, S_biased, and S_unbiased.
#' @keywords internal
Step1_ML <- function(S1, s1_fit = NULL, centered, group) {
  start_time_step1 <- Sys.time()

  M <- length(S1$mplus_models)
  nfactors <- length(S1$lat_var)
  latvar_mplus <- toupper(S1$lat_var)
  var_mplus <- S1$vars

  res_temp <- vector(mode = "list", length = M)
  if (is.null(s1_fit)) {
    s1_fit <- vector(mode = "list", length = M)
  }
  post_est_ordered <- vector(mode = "list", length = M)
  matching_columns <- matching_columns_sd <- list()

  for (i in seq_len(M)) {
    if (is.null(s1_fit[[i]])) {
      mod <- mplusObject(
        TITLE = "Multilevel CFA;",
        VARIABLE = paste("USEVARIABLES =", group, var_mplus[[i]], ";",
                         "\nWITHIN =", var_mplus[[i]], ";",
                         "\nCLUSTER =", group, ";"),
        DEFINE = paste0("CENTER ", paste(var_mplus[[i]], "(GROUPMEAN);")),
        ANALYSIS = "TYPE = RANDOM TWOLEVEL;\nESTIMATOR = BAYES;\nBITERATIONS = 10000;\nCONVERGENCE = 0.05;",
        MODEL = S1$mplus_models[[i]],
        SAVEDATA = paste0("\nFILE IS ", paste0("Fscores", i, ".dat"), ";\nSAVE = FSCORES(100 10);\nFORMAT IS FREE;"),
        OUTPUT = "FSCOMPARISON; ",
        rdata = centered
      )
      res_temp[[i]] <- mplusModeler(mod, modelout = paste0("cfa", i, ".inp"), run = 1L)
      s1_fit[[i]] <- res_temp[[i]][["results"]]
    }

    col_names <- colnames(s1_fit[[i]][["savedata"]])
    matched_col <- col_names[toupper(group) == col_names]
    post_est_ordered[[i]] <- s1_fit[[i]][["savedata"]][order(s1_fit[[i]][["savedata"]][[matched_col]]), ]
    matching_columns[[i]] <- grep(paste0("^( ", paste(latvar_mplus, collapse = "|") , ")( %W)?\\s*Mean$"), names(post_est_ordered[[i]]), value = TRUE)
    matching_columns_sd[[i]] <- grep(paste0("^( ", paste(latvar_mplus, collapse = "|") , ")( %W)?\\s*Standard D*"), names(post_est_ordered[[i]]), value = TRUE)
  }

  postmeanfs <- do.call(cbind, lapply(seq_along(post_est_ordered), function(i) {
    post_est_ordered[[i]][, matching_columns[[i]], drop = FALSE]
  }))
  fscores <- data.frame(cbind(postmeanfs, post_est_ordered[[1]][[matched_col]]))
  vars <- paste0("single", seq_len(nfactors))
  colnames(fscores) <- c(vars, "group")
  fcovdat <- lapply(split(fscores[, vars], fscores$group), cov)
  fcovdat <- array(unlist(fcovdat), dim = c(nfactors, nfactors, length(fcovdat)))
  postmeanfs_list <- split(fscores[, vars], fscores$group)

  postvarfs <- cbind((do.call(cbind, lapply(seq_along(post_est_ordered), function(i) {
    post_est_ordered[[i]][, matching_columns_sd[[i]], drop = FALSE]
  })))^2, post_est_ordered[[1]][[matched_col]])
  colnames(postvarfs) <- c(vars, "group")
  postvarfs_list <- split(postvarfs[, vars], postvarfs$group)

  fvar_compute <- replicate(nfactors, numeric(length(postmeanfs_list)), simplify = FALSE)
  for (i in seq_len(nfactors)) {
    for (j in seq_along(postmeanfs_list)) {
      fvar_compute[[i]][j] <- var(postmeanfs_list[[j]][, i]) + mean(postvarfs_list[[j]][, i])
    }
  }

  rho <- lapply(seq_len(nfactors), function(i) {
    sapply(seq_along(postmeanfs_list), function(j) var(postmeanfs_list[[j]][, i]) / fvar_compute[[i]][j])
  })

  lambda <- array(0, dim = c(nfactors, nfactors, length(postmeanfs_list)))
  for (j in seq_along(postmeanfs_list)) {
    lambda[, , j] <- diag(sapply(rho, `[`, j))
  }
  dimnames(lambda) <- list(vars, S1$lat_var, NULL)

  resvar <- lapply(seq_along(postmeanfs_list), function(j) {
    sapply(seq_len(nfactors), function(i) fvar_compute[[i]][j] * rho[[i]][j] * (1 - rho[[i]][j]))
  })
  theta <- array(unlist(lapply(resvar, function(x) diag(x))), dim = c(nfactors, nfactors, length(resvar)))
  dimnames(theta) <- list(vars, vars, NULL)

  cov_eta <- array(0, dim = c(nfactors, nfactors, length(resvar)))
  for (j in seq_along(resvar)) {
    cov_eta[, , j] <- solve(lambda[, , j]) %*% (fcovdat[, , j] - theta[, , j]) %*% solve(t(lambda[, , j]))
  }
  dimnames(cov_eta) <- list(S1$lat_var, S1$lat_var, NULL)

  centered_fs <- fscores
  g_name <- as.character(unique(fscores$group))
  group.idx <- match(fscores$group, g_name)
  group.sizes <- tabulate(group.idx)
  group.means <- rowsum.default(as.matrix(fscores[, vars]), group = group.idx, reorder = FALSE, na.rm = FALSE) / group.sizes
  centered_fs[, vars] <- fscores[, vars] - group.means[group.idx, , drop = FALSE]

  S_unbiased <- lapply(unique(centered_fs$group), function(x) {
    cov(centered_fs[centered_fs$group == x, vars])
  })

  N_gs <- as.numeric(table(centered_fs$group))
  ngroups <- length(N_gs)
  S_biased <- lapply(seq_along(S_unbiased), function(g) {
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
