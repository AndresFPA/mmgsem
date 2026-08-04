#' Internal helper: fast EM Step 2 estimation for MMGSEM.
#'
#' Performs the fast EM algorithm for Step 2, estimating structural parameters
#' and cluster posterior probabilities for the MMGSEM model.
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
#' @param s1ori original step 1 syntax before any ordering transformation.
#' @param s1_type step 1 estimation type.
#' @return A list containing the EM results, including fitted models, loglikelihoods, and posterior probabilities.
#' @keywords internal
Step2_fast <- function(ngroups, nclus, nstarts, N_gs, seed, max_it,
                       cov_eta, dat, S2, lat_var, ordered,
                       endo_group_specific, endogenous_cov, lambda_gs,
                       theta_gs, S_unbiased, S1, std.lv,
                       partition, userStart, printing, s1ori, s1_type) {
  results_nstarts <- vector(mode = "list", length = nstarts)
  z_gks_nstarts <- vector(mode = "list", length = nstarts)
  loglik_nstarts <- numeric(nstarts)
  iter_nstarts <- numeric(nstarts)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (ordered) {
    S1 <- s1ori
  }

  if (endogenous_cov) {
    fake <- lavaan::sem(
      model = S2,
      sample.cov = rep(cov_eta[1], nclus),
      sample.nobs = rep(nrow(dat), nclus),
      do.fit = FALSE,
      baseline = FALSE,
      h1 = FALSE,
      check.post = FALSE,
      loglik = FALSE,
      sample.cov.rescale = FALSE,
      fixed.x = TRUE,
      information = "observed"
    )
  } else {
    fake <- lavaan::sem(
      model = S2,
      sample.cov = rep(cov_eta[1], nclus),
      sample.nobs = rep(nrow(dat), nclus),
      do.fit = FALSE,
      baseline = FALSE,
      h1 = FALSE,
      check.post = FALSE,
      loglik = FALSE,
      sample.cov.rescale = FALSE,
      fixed.x = TRUE,
      auto.cov.y = FALSE,
      information = "observed"
    )
  }

  FakeprTbl <- lavaan::parTable(fake)
  fake@Options$do.fit <- TRUE
  fake@Options$se <- "none"
  fake@ParTable$start <- NULL
  fake@ParTable$est <- NULL
  fake@ParTable$se <- NULL
  fake@Options$start <- "default"

  endog1 <- lat_var[(lat_var %in% FakeprTbl$rhs[FakeprTbl$op == "~"]) &
                     (lat_var %in% FakeprTbl$lhs[FakeprTbl$op == "~"])]
  endog2 <- lat_var[!c(lat_var %in% FakeprTbl$rhs[FakeprTbl$op == "~"]) &
                     (lat_var %in% FakeprTbl$lhs[FakeprTbl$op == "~"])]
  endog <- c(endog1, endog2)
  exog <- lat_var[!c(lat_var %in% endog)]

  fake_lv <- vector(mode = "list", length = length(endog))
  prTbl_lv <- vector(mode = "list", length = length(endog))

  for (lv in seq_along(endog)) {
    this_lv <- endog[lv]
    var_not_this_lv <- which(FakeprTbl$lhs != this_lv & FakeprTbl$op == "~~")
    fac_load <- which(FakeprTbl$lhs != this_lv & FakeprTbl$op == "=~")
    prTbl_idx <- sort(c(which(FakeprTbl$lhs == this_lv), fac_load, var_not_this_lv))
    prTbl_lv[[lv]] <- FakeprTbl[prTbl_idx, ]

    fake_lv[[lv]] <- lavaan::sem(
      model = prTbl_lv[[lv]],
      sample.cov = rep(cov_eta[1], nclus),
      sample.nobs = rep(nrow(dat), nclus),
      do.fit = FALSE,
      baseline = FALSE,
      h1 = FALSE,
      check.post = FALSE,
      loglik = FALSE,
      sample.cov.rescale = FALSE,
      fixed.x = TRUE,
      information = "observed"
    )
    fake_lv[[lv]]@Options$do.fit <- TRUE
    fake_lv[[lv]]@Options$se <- "none"
    fake_lv[[lv]]@ParTable$start <- NULL
    fake_lv[[lv]]@ParTable$est <- NULL
    fake_lv[[lv]]@ParTable$se <- NULL
    fake_lv[[lv]]@Options$start <- "default"
  }

  cov_eta <- lapply(seq_len(ngroups), function(x) {
    reorder(cov_eta[[x]], exog = exog, endog = endog)
  })

  if (s1_type == "lavaan") {
    lambda_gs <- lapply(seq_len(ngroups), function(x) {
      reorder_obs(lambda_gs[[x]], matrix = "lambda", exog = exog, endog = endog,
                  endog1 = endog1, endog2 = endog2, S1 = S1, dat = dat)
    })
    theta_gs <- lapply(seq_len(ngroups), function(x) {
      reorder_obs(theta_gs[[x]], matrix = "theta", exog = exog, endog = endog,
                  endog1 = endog1, endog2 = endog2, S1 = S1, dat = dat)
    })
  }

  for (s in seq_len(nstarts)) {
    if (printing) {
      message("Start ", s, " -----------------")
    }

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

    psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
    z_gks_lv <- lapply(seq_along(endog), function(x) matrix(NA, nrow = ngroups, ncol = nclus))

    i <- 0
    prev_LL <- 0
    diff_LL <- 1
    log_test <- TRUE

    while (diff_LL > 1e-6 && i < max_it && isTRUE(log_test)) {
      i <- i + 1
      pi_ks <- colMeans(z_gks)
      N_gks <- z_gks * N_gs

      if (endo_group_specific && i > 1) {
        for (lv in seq_along(endog)) {
          for (k in seq_len(nclus)) {
            for (g in seq_len(ngroups)) {
              z_gks_lv[[lv]][g, k] <- z_gks[g, k] / psi_gks[[g, k]][endog[lv], endog[lv]]
            }
          }
        }
      }

      if (!endo_group_specific || i == 1) {
        COV <- vector("list", length = nclus)
        for (k in seq_len(nclus)) {
          this_nobs <- z_gks[, k] * N_gs
          this_w <- this_nobs / sum(this_nobs)
          tmp <- lapply(seq_along(cov_eta), function(g) cov_eta[[g]] * this_w[g])
          COV[[k]] <- Reduce("+", tmp)
        }
      } else {
        COV_lv <- lapply(seq_along(endog), function(x) vector("list", length = nclus))
        for (lv in seq_along(endog)) {
          for (k in seq_len(nclus)) {
            this_nobs <- z_gks_lv[[lv]][, k] * N_gs
            this_w <- this_nobs / sum(this_nobs)
            tmp <- lapply(seq_along(cov_eta), function(g) cov_eta[[g]] * this_w[g])
            COV_lv[[lv]][[k]] <- Reduce("+", tmp)
          }
        }
      }

      if (!endo_group_specific || i == 1) {
        s2out <- lavaan::lavaan(
          slotOptions = fake@Options,
          slotParTable = fake@ParTable,
          sample.cov = COV,
          sample.nobs = rep(nrow(dat), nclus)
        )
      } else {
        s2out <- lapply(seq_along(endog), function(lv) {
          lavaan::lavaan(
            slotOptions = fake_lv[[lv]]@Options,
            slotParTable = fake_lv[[lv]]@ParTable,
            sample.cov = COV_lv[[lv]],
            sample.nobs = rep(nrow(dat), nclus)
          )
        })
      }

      loglik_gks <- matrix(0, nrow = ngroups, ncol = nclus)
      loglik_gksw <- matrix(0, nrow = ngroups, ncol = nclus)
      Sigma <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
      I <- diag(length(lat_var))

      if (!endo_group_specific || i == 1) {
        if (nclus == 1) {
          EST_s2 <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
          beta_ks <- EST_s2[["beta"]]
          psi_ks <- EST_s2[["psi"]]
        } else {
          EST_s2 <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
          beta_ks <- lapply(EST_s2, "[[", "beta")
          psi_ks <- lapply(EST_s2, "[[", "psi")
        }
        if (nclus == 1) {
          beta_ks <- reorder(beta_ks, exog = exog, endog = endog)
          psi_ks <- reorder(psi_ks, exog = exog, endog = endog)
        } else {
          beta_ks <- lapply(seq_len(nclus), function(x) reorder(beta_ks[[x]], exog = exog, endog = endog))
          psi_ks <- lapply(seq_len(nclus), function(x) reorder(psi_ks[[x]], exog = exog, endog = endog))
        }
      } else {
        EST_s2_lv <- vector(mode = "list", length = length(endog))
        beta_ks_lv <- vector(mode = "list", length = length(endog))
        psi_ks_lv <- vector(mode = "list", length = length(endog))
        for (lv in seq_along(endog)) {
          if (nclus == 1) {
            EST_s2_lv[[lv]] <- lavaan::lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
            beta_ks_lv[[lv]] <- EST_s2_lv[[lv]][["beta"]]
            psi_ks_lv[[lv]] <- EST_s2_lv[[lv]][["psi"]]
          } else {
            EST_s2_lv[[lv]] <- lavaan::lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
            beta_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "beta")
            psi_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "psi")
          }
        }

        for (k in seq_len(nclus)) {
          for (lv in seq_along(endog)) {
            this_lv <- endog[lv]
            col.idx <- colnames(beta_ks_lv[[lv]][[k]])
            beta_ks[[k]][this_lv, col.idx] <- beta_ks_lv[[lv]][[k]][this_lv, col.idx]
          }
        }

        if (nclus == 1) {
          beta_ks <- reorder(beta_ks, exog = exog, endog = endog)
        } else {
          beta_ks <- lapply(seq_len(nclus), function(x) reorder(beta_ks[[x]], exog = exog, endog = endog))
        }
      }

      for (k in seq_len(nclus)) {
        psi_k <- if (nclus == 1) psi_ks else psi_ks[[k]]
        beta <- if (nclus == 1) beta_ks else beta_ks[[k]]
        for (g in seq_len(ngroups)) {
          psi <- psi_k
          psi[exog, exog] <- cov_eta[[g]][exog, exog]
          if (endo_group_specific) {
            solved_psi <- (I - beta) %*% cov_eta[[g]] %*% t(I - beta)
            g_endog1_cov <- solved_psi[endog1, endog1]
            if (length(endog1) > 1) {
              g_endog1_cov[row(g_endog1_cov) != col(g_endog1_cov)] <- 0
            }
            psi[endog1, endog1] <- g_endog1_cov
            psi[endog2, endog2] <- solved_psi[endog2, endog2]
          }
          if (!endogenous_cov) {
            offdiag <- row(psi[endog2, endog2]) != col(psi[endog2, endog2])
            psi[endog2, endog2][offdiag] <- 0
          }

          psi_gks[[g, k]] <- psi
          Sigma[[g, k]] <- solve(I - beta) %*% psi %*% t(solve(I - beta))
          Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]]))

          loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
            sample.mean = rep(0, length(lat_var)),
            sample.nobs = N_gs[g],
            sample.cov = cov_eta[[g]],
            Mu = rep(0, length(lat_var)),
            Sigma = Sigma[[g, k]]
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
      if (i == 1) log_test <- TRUE
      prev_LL <- LL
      if (printing) {
        message(i, " ", LL)
      }
    }

    results_nstarts[[s]] <- s2out
    z_gks_nstarts[[s]] <- z_gks
    loglik_nstarts[s] <- LL
    iter_nstarts[s] <- i
  }

  best_idx <- which.max(loglik_nstarts)
  iter <- iter_nstarts[best_idx]
  s2out <- results_nstarts[[best_idx]]
  LL <- loglik_nstarts[best_idx]
  z_gks <- z_gks_nstarts[[best_idx]]
  colnames(z_gks) <- paste("Cluster", seq_len(nclus))

  if (!endo_group_specific) {
    if (nclus == 1) {
      EST_s2 <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
      beta_ks <- EST_s2[["beta"]]
      psi_ks <- EST_s2[["psi"]]
    } else {
      EST_s2 <- lavaan::lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
      beta_ks <- lapply(EST_s2, "[[", "beta")
      psi_ks <- lapply(EST_s2, "[[", "psi")
    }
  } else {
    EST_s2_lv <- vector(mode = "list", length = length(endog))
    beta_ks_lv <- vector(mode = "list", length = length(endog))
    psi_ks_lv <- vector(mode = "list", length = length(endog))
    for (lv in seq_along(endog)) {
      if (nclus == 1) {
        EST_s2_lv[[lv]] <- lavaan::lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
        beta_ks_lv[[lv]] <- EST_s2_lv[[lv]][["beta"]]
        psi_ks_lv[[lv]] <- EST_s2_lv[[lv]][["psi"]]
      } else {
        EST_s2_lv[[lv]] <- lavaan::lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
        beta_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "beta")
        psi_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "psi")
      }
    }

    for (k in seq_len(nclus)) {
      for (lv in seq_along(endog)) {
        this_lv <- endog[lv]
        col.idx <- colnames(beta_ks_lv[[lv]][[k]])
        if (nclus == 1) {
          beta_ks[this_lv, col.idx] <- beta_ks_lv[[lv]][this_lv, col.idx]
        } else {
          beta_ks[[k]][this_lv, col.idx] <- beta_ks_lv[[lv]][[k]][this_lv, col.idx]
        }
      }
    }
  }

  if (nclus == 1) {
    beta_ks <- reorder(beta_ks, exog = exog, endog = endog)
  } else {
    beta_ks <- lapply(seq_len(nclus), function(x) reorder(beta_ks[[x]], exog = exog, endog = endog))
  }

  psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
  for (k in seq_len(nclus)) {
    psi_k <- if (nclus == 1) psi_ks else psi_ks[[k]]
    beta <- if (nclus == 1) beta_ks else beta_ks[[k]]
    for (g in seq_len(ngroups)) {
      psi <- psi_k
      psi[exog, exog] <- cov_eta[[g]][exog, exog]
      if (endo_group_specific) {
        g_endog1_cov <- ((I - beta) %*% cov_eta[[g]] %*% t(I - beta))[endog1, endog1]
        if (length(endog1) > 1) {
          g_endog1_cov[row(g_endog1_cov) != col(g_endog1_cov)] <- 0
        }
        psi[endog1, endog1] <- g_endog1_cov
        psi[endog2, endog2] <- ((I - beta) %*% cov_eta[[g]] %*% t(I - beta))[endog2, endog2]
      }
      if (!endogenous_cov) {
        offdiag <- row(psi[endog2, endog2]) != col(psi[endog2, endog2])
        psi[endog2, endog2][offdiag] <- 0
      }
      psi_gks[[g, k]] <- psi
    }
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
