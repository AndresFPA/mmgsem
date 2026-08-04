#' Internal helper: slow EM Step 2 estimation for MMGSEM.
#'
#' Performs the slow, fully group-cluster-specific EM estimation for Step 2.
#' This is used when endogenous covariance among dependent latent variables
#' must be estimated and the fast approximation is not sufficient.
#'
#' @param ngroups number of groups.
#' @param nclus number of clusters.
#' @param nstarts number of random starts.
#' @param N_gs group sample sizes.
#' @param seed random seed.
#' @param max_it maximum number of EM iterations.
#' @param cov_eta list of group-specific factor covariance matrices.
#' @param dat data frame used for step 2.
#' @param S2 lavaan syntax for the structural model.
#' @param lat_var latent variable names.
#' @param ordered logical indicating ordered observed data.
#' @param endo_group_specific logical for group-specific endogenous covariances.
#' @param endogenous_cov logical for endogenous covariance estimation.
#' @param lambda_gs measurement model loading matrices.
#' @param theta_gs measurement model theta matrices.
#' @param S_unbiased unbiased observed covariance matrices.
#' @param S1 step 1 model syntax.
#' @param std.lv logical for latent variable scaling.
#' @param partition initialization method.
#' @param userStart optional user-defined start matrix.
#' @param printing logical controlling verbose output.
#' @param s1ori original step 1 syntax before ordering transformation.
#' @param only_slow logical forcing slow estimation path.
#' @param beta_ks optional starting beta parameter matrices.
#' @param psi_gks optional starting psi parameter matrices.
#' @param z_gks optional starting posterior matrix.
#' @return A list containing slow EM step results and fit details.
#' @keywords internal
Step2_slow <- function(ngroups, nclus, nstarts, N_gs, seed, max_it,
                       cov_eta, dat, S2, lat_var, ordered,
                       endo_group_specific, endogenous_cov, lambda_gs,
                       theta_gs, S_unbiased, S1, std.lv,
                       partition, userStart, printing, s1ori,
                       only_slow, beta_ks = NULL, psi_gks = NULL, z_gks = NULL) {
  gro_clu <- nclus * ngroups

  fake_cov <- rep(cov_eta, nclus)
  names(fake_cov) <- paste("group", seq_len(gro_clu))
  fake_model <- lavaan::parTable(lavaan::sem(
    model = S2,
    sample.cov = fake_cov,
    sample.nobs = rep(N_gs, nclus),
    do.fit = FALSE,
    meanstructure = FALSE,
    h1 = FALSE,
    check.post = FALSE,
    loglik = FALSE,
    sample.cov.rescale = FALSE,
    fixed.x = TRUE
  ))

  endog1 <- lat_var[(lat_var %in% fake_model$rhs[fake_model$op == "~"]) &
                     (lat_var %in% fake_model$lhs[fake_model$op == "~"])]
  endog2 <- lat_var[!c(lat_var %in% fake_model$rhs[fake_model$op == "~"]) &
                     (lat_var %in% fake_model$lhs[fake_model$op == "~"])]
  endog <- c(endog1, endog2)
  exog <- lat_var[!c(lat_var %in% endog)]

  fake_model$par <- paste0(fake_model$lhs, fake_model$op, fake_model$rhs, ".g", fake_model$group)
  fake_model$cluster <- rep(1:nclus, each = length(fake_model$id[fake_model$group %in% 1:ngroups]))
  fake_model$se <- NULL
  fake_model$cluster <- NULL
  fake_model$par <- NULL
  fake_model$est <- NULL
  fake_model$start <- NULL

  constraints_row <- data.frame(
    id = "", lhs = "", op = "==", rhs = "",
    user = 2, block = 0, group = 0, free = 0,
    ustart = NA, exo = 0, label = "", plabel = "",
    cluster = NA,
    stringsAsFactors = FALSE
  )

  constraints <- fake_model$plabel[fake_model$op == "~"]
  n_reg <- length(fake_model$plabel[fake_model$op == "~" & fake_model$group == 1])
  cons_exo <- fake_model$plabel[fake_model$op == "~~" & fake_model$lhs %in% exog]
  n_exo <- length(fake_model$plabel[fake_model$op == "~~" & fake_model$lhs %in% exog & fake_model$group == 1])

  constraints_matrix <- constraints_row[rep(seq_len(nrow(constraints_row)), length(constraints)), ]
  cons_exo_matrix <- constraints_row[rep(seq_len(nrow(constraints_row)), length(cons_exo)), ]
  rownames(constraints_matrix) <- NULL
  rownames(cons_exo_matrix) <- NULL

  clus_label <- rep(1:nclus, each = ngroups)
  group_label <- rep(1:ngroups, times = nclus)

  for (j in seq_along(clus_label)) {
    fake_model$cluster[fake_model$group == j] <- clus_label[j]
  }

  reg_labels <- rep(clus_label, each = n_reg)
  exo_labels <- rep(group_label, each = n_exo)

  for (k in seq_len(nclus)) {
    cluster_par <- constraints[reg_labels == k]
    constraints_matrix[reg_labels == k, "lhs"] <- cluster_par[seq_len(n_reg)]
    constraints_matrix[reg_labels == k, "rhs"] <- cluster_par
  }

  for (g in seq_len(ngroups)) {
    group_par_exo <- cons_exo[exo_labels == g]
    cons_exo_matrix[exo_labels == g, "lhs"] <- group_par_exo[seq_len(n_exo)]
    cons_exo_matrix[exo_labels == g, "rhs"] <- group_par_exo
  }

  constraints_total <- rbind(constraints_matrix, cons_exo_matrix)
  redundant <- which(constraints_total$lhs == constraints_total$rhs)
  if (length(redundant) > 0) {
    constraints_total <- constraints_total[-redundant, ]
  }
  rownames(constraints_total) <- NULL

  fake_model <- rbind(fake_model, constraints_total)
  fake_model$free <- seq_len(nrow(fake_model))

  if (!only_slow && !is.null(beta_ks) && !is.null(psi_gks) && !is.null(z_gks)) {
    fake_model$parK <- paste0(fake_model$lhs, fake_model$op, fake_model$rhs, ".k", fake_model$cluster)
    beta_vec <- c(); beta_nam <- c()
    for (k in seq_len(nclus)) {
      beta <- if (nclus == 1) beta_ks else beta_ks[[k]]
      non_zer.idx <- which(unlist(beta) != 0)
      beta_vec <- c(beta_vec, unlist(beta)[non_zer.idx])
      beta_nam <- c(beta_nam, as.vector(outer(rownames(beta), colnames(beta), function(x, y) paste0(x, "~", y, ".k", k)))[non_zer.idx])
    }
    beta_vec <- setNames(beta_vec, beta_nam)
    beta.idx <- match(fake_model$parK, names(beta_vec))
    fake_model$ustart <- ifelse(!is.na(beta.idx), beta_vec[beta.idx], fake_model$ustart)

    cov_vec <- c(); cov_nam <- c(); gk <- 0
    for (k in seq_len(nclus)) {
      for (g in seq_len(ngroups)) {
        gk <- gk + 1
        tmp.lower.tri <- psi_gks[[g, k]]
        tmp.lower.tri[upper.tri(tmp.lower.tri)] <- 0
        non_zer.idx <- which(unlist(tmp.lower.tri) != 0)
        unique.idx <- !duplicated(unlist(tmp.lower.tri)[non_zer.idx])
        cov_vec <- c(cov_vec, unlist(tmp.lower.tri)[non_zer.idx][unique.idx])
        cov_nam <- c(cov_nam, as.vector(outer(rownames(psi_gks[[g, k]]), colnames(psi_gks[[g, k]]), function(x, y) paste0(y, "~~", x, ".g", gk)))[non_zer.idx][unique.idx])
      }
    }
    cov_vec <- setNames(cov_vec, cov_nam)
    cov.idx <- match(fake_model$par, names(cov_vec))
    fake_model$ustart <- ifelse(!is.na(cov.idx), cov_vec[cov.idx], fake_model$ustart)
    fake_model$ustart[is.na(fake_model$ustart)] <- 0.01
  }

  fake_model$se <- NULL
  fake_model$cluster <- NULL
  fake_model$par <- NULL
  fake_model$parK <- NULL
  fake_model$est <- NULL
  fake_model$start <- NULL

  if (!only_slow) {
    fake_model$par <- paste0(fake_model$lhs, fake_model$op, fake_model$rhs, ".g", fake_model$group)
    free_idx <- which(fake_model$rhs %in% lat_var)
    fake_model$free[free_idx] <- seq_along(free_idx)

    constraints_row2 <- constraints_row
    constraints <- fake_model$plabel[fake_model$op == "~"]
    n_reg <- length(fake_model$plabel[fake_model$op == "~" & fake_model$group == 1])
    cons_exo <- fake_model$plabel[fake_model$op == "~~" & fake_model$lhs %in% exog]
    n_exo <- length(fake_model$plabel[fake_model$op == "~~" & fake_model$lhs %in% exog & fake_model$group == 1])

    constraints_matrix <- constraints_row2[rep(seq_len(nrow(constraints_row2)), length(constraints)), ]
    cons_exo_matrix <- constraints_row2[rep(seq_len(nrow(constraints_row2)), length(cons_exo)), ]
    rownames(constraints_matrix) <- NULL
    rownames(cons_exo_matrix) <- NULL

    clus_label <- rep(seq_len(nclus), each = ngroups)
    group_label <- rep(seq_len(ngroups), times = nclus)
    for (j in seq_along(clus_label)) fake_model$cluster[fake_model$group == j] <- clus_label[j]

    reg_labels <- rep(clus_label, each = n_reg)
    exo_labels <- rep(group_label, each = n_exo)
    for (k in seq_len(nclus)) {
      cluster_par <- constraints[reg_labels == k]
      constraints_matrix[reg_labels == k, "lhs"] <- cluster_par[seq_len(n_reg)]
      constraints_matrix[reg_labels == k, "rhs"] <- cluster_par
    }
    for (g in seq_len(ngroups)) {
      group_par_exo <- cons_exo[exo_labels == g]
      cons_exo_matrix[exo_labels == g, "lhs"] <- group_par_exo[seq_len(n_exo)]
      cons_exo_matrix[exo_labels == g, "rhs"] <- group_par_exo
    }
    constraints_total <- rbind(constraints_matrix, cons_exo_matrix)
    redundant <- which(constraints_total$lhs == constraints_total$rhs)
    if (length(redundant) > 0) constraints_total <- constraints_total[-redundant, ]
    rownames(constraints_total) <- NULL
    fake_model <- rbind(fake_model, constraints_total)
  }

  if (!only_slow) {
    if (!is.null(seed)) set.seed(seed)
    pi_ks <- colMeans(z_gks)
    N_gks <- z_gks * N_gs
    N_gks <- c(N_gks)
    s2out <- lavaan::sem(
      model = fake_model,
      sample.cov = fake_cov,
      sample.nobs = N_gks,
      baseline = FALSE,
      se = "none",
      h1 = FALSE,
      check.post = FALSE,
      control = list(rel.tol = 1e-09),
      sample.cov.rescale = FALSE,
      fixed.x = FALSE
    )

    loglik_gks <- matrix(0, nrow = ngroups, ncol = nclus)
    loglik_gksw <- matrix(0, nrow = ngroups, ncol = nclus)
    for (k in seq_len(nclus)) {
      for (g in seq_len(ngroups)) {
        gk <- (k - 1) * ngroups + g
        loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
          sample.mean = s2out@SampleStats@mean[[gk]],
          sample.nobs = N_gs[g],
          sample.cov = s2out@SampleStats@cov[[gk]],
          Mu = s2out@SampleStats@mean[[gk]],
          Sigma = s2out@implied$cov[[gk]]
        )
        loglik_gks[g, k] <- loglik_gk
        loglik_gksw[g, k] <- log(pi_ks[k]) + loglik_gk
      }
    }

    max_gs <- apply(loglik_gksw, 1, max)
    minus_max <- sweep(loglik_gksw, 1, max_gs, "-")
    exp_loglik <- exp(minus_max)
    loglik_gsw <- log(rowSums(exp_loglik))
    LL <- sum(loglik_gsw + max_gs)

    z_gks <- EStep(pi_ks = pi_ks, ngroup = ngroups, nclus = nclus, loglik = loglik_gks)
    colnames(z_gks) <- paste("Cluster", seq_len(nclus))

    EST_s2 <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
    beta_gks <- lapply(EST_s2, "[[", "beta")
    psi_gks_tmp <- lapply(EST_s2, "[[", "psi")
    k.idx <- (seq_len(nclus) - 1) * ngroups + 1L
    beta_ks <- beta_gks[k.idx]

    if (nclus == 1) {
      beta_ks <- reorder(beta_ks[[1]], exog = exog, endog = endog)
    } else {
      beta_ks <- lapply(seq_len(nclus), function(x) reorder(beta_ks[[x]], exog = exog, endog = endog))
    }

    psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
    for (gk in seq_len(gro_clu)) {
      psi_gks_tmp[[gk]] <- reorder(psi_gks_tmp[[gk]], exog = exog, endog = endog)
      psi_gks[[gk]] <- psi_gks_tmp[[gk]]
    }

    return(list(
      z_gks = z_gks,
      LL = LL,
      s2out = s2out,
      beta_ks = beta_ks,
      psi_gks = psi_gks
    ))
  }

  if (!is.null(seed)) set.seed(seed)
  results_nstarts <- vector(mode = "list", length = nstarts)
  z_gks_nstarts <- vector(mode = "list", length = nstarts)
  loglik_nstarts <- numeric(nstarts)
  iter_nstarts <- numeric(nstarts)

  for (s in seq_len(nstarts)) {
    if (printing) message("Start ", s, " -----------------")

    if (!is.null(userStart)) {
      z_gks <- userStart
    } else if (partition == "hard") {
      cl <- 0
      while (cl < 1) {
        z_gks <- t(replicate(ngroups, sample(x = c(rep(0, nclus - 1), 1))))
        cl <- min(colSums(z_gks))
      }
    } else {
      z_gks <- matrix(runif(nclus * ngroups), ncol = nclus, nrow = ngroups)
      z_gks <- z_gks / rowSums(z_gks)
    }

    iter <- 0; prev_LL <- 0; diff_LL <- 1; log_test <- TRUE
    while (diff_LL > 1e-6 && iter < max_it && isTRUE(log_test)) {
      iter <- iter + 1
      pi_ks <- colMeans(z_gks)
      N_gks <- z_gks * N_gs
      N_gks <- c(N_gks)

      if (iter == 1) {
        s2out <- lavaan::sem(
          model = fake_model,
          sample.cov = fake_cov,
          sample.nobs = N_gks,
          baseline = FALSE,
          se = "none",
          h1 = FALSE,
          check.post = FALSE,
          control = list(rel.tol = 1e-09),
          sample.cov.rescale = FALSE,
          fixed.x = FALSE
        )
      } else {
        s2out <- lavaan::sem(
          model = fake_model,
          sample.cov = fake_cov,
          sample.nobs = N_gks,
          start = start,
          baseline = FALSE,
          se = "none",
          h1 = FALSE,
          check.post = FALSE,
          control = list(rel.tol = 1e-06),
          sample.cov.rescale = FALSE,
          fixed.x = FALSE
        )
      }

      start <- partable(s2out)$est
      loglik_gks <- matrix(0, nrow = ngroups, ncol = nclus)
      loglik_gksw <- matrix(0, nrow = ngroups, ncol = nclus)
      for (k in seq_len(nclus)) {
        for (g in seq_len(ngroups)) {
          gk <- (k - 1) * ngroups + g
          loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
            sample.mean = s2out@SampleStats@mean[[gk]],
            sample.nobs = N_gs[g],
            sample.cov = s2out@SampleStats@cov[[gk]],
            Mu = s2out@SampleStats@mean[[gk]],
            Sigma = s2out@implied$cov[[gk]]
          )
          loglik_gks[g, k] <- loglik_gk
          loglik_gksw[g, k] <- log(pi_ks[k]) + loglik_gk
        }
      }

      max_gs <- apply(loglik_gksw, 1, max)
      minus_max <- sweep(loglik_gksw, 1, max_gs, "-")
      exp_loglik <- exp(minus_max)
      loglik_gsw <- log(rowSums(exp_loglik))
      LL <- sum(loglik_gsw + max_gs)

      z_gks <- EStep(pi_ks = pi_ks, ngroup = ngroups, nclus = nclus, loglik = loglik_gks)
      diff_LL <- abs(LL - prev_LL)
      log_test <- prev_LL < LL || isTRUE(all.equal(prev_LL, LL))
      if (iter == 1) log_test <- TRUE
      prev_LL <- LL
      if (printing) message(iter, " ", LL)
    }

    results_nstarts[[s]] <- s2out
    z_gks_nstarts[[s]] <- z_gks
    loglik_nstarts[s] <- LL
    iter_nstarts[s] <- iter
  }

  best_idx <- which.max(loglik_nstarts)
  iter <- iter_nstarts[best_idx]
  s2out <- results_nstarts[[best_idx]]
  LL <- loglik_nstarts[best_idx]
  z_gks <- z_gks_nstarts[[best_idx]]
  colnames(z_gks) <- paste("Cluster", seq_len(nclus))

  EST_s2 <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
  beta_gks <- lapply(EST_s2, "[[", "beta")
  psi_gks_tmp <- lapply(EST_s2, "[[", "psi")

  k.idx <- (seq_len(nclus) - 1) * ngroups + 1L
  beta_ks <- beta_gks[k.idx]

  if (nclus == 1) {
    beta_ks <- reorder(beta_ks[[1]], exog = exog, endog = endog)
  } else {
    beta_ks <- lapply(seq_len(nclus), function(x) reorder(beta_ks[[x]], exog = exog, endog = endog))
  }

  psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
  for (gk in seq_len(gro_clu)) {
    psi_gks_tmp[[gk]] <- reorder(psi_gks_tmp[[gk]], exog = exog, endog = endog)
    psi_gks[[gk]] <- psi_gks_tmp[[gk]]
  }

  return(list(
    results_nstarts = results_nstarts,
    z_gks_nstarts = z_gks_nstarts,
    loglik_nstarts = loglik_nstarts,
    iter_nstarts = iter_nstarts,
    iter = iter,
    z_gks = z_gks,
    LL = LL,
    s2out = s2out,
    endog = endog,
    exog = exog,
    endog1 = endog1,
    endog2 = endog2,
    beta_ks = beta_ks,
    psi_gks = psi_gks
  ))
}
