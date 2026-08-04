#' Internal helper: wrapper for Step 2 estimation that selects the fast or slow EM algorithm.
#'
#' Selects the appropriate Step 2 routine for MMGSEM based on model structure
#' and the presence of endogenous covariance among dependent latent variables.
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
#' @param only_slow logical forcing slow estimation.
#' @param s1_type step 1 estimation type.
#' @return A list containing the selected step 2 results and convergence details.
#' @keywords internal
Step2 <- function(ngroups, nclus, nstarts, N_gs, seed, max_it,
                  cov_eta, dat, S2, lat_var, ordered,
                  endo_group_specific, endogenous_cov, lambda_gs,
                  theta_gs, S_unbiased, S1, std.lv, partition,
                  userStart, printing, s1ori, only_slow, s1_type) {
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

  if (length(endog2) > 1 && endogenous_cov) {
    if (!only_slow) {
      SM_1 <- Step2_fast(
        ngroups = ngroups,
        nclus = nclus,
        nstarts = nstarts,
        N_gs = N_gs,
        seed = seed,
        max_it = max_it,
        cov_eta = cov_eta,
        dat = dat,
        S2 = S2,
        lat_var = lat_var,
        ordered = ordered,
        endo_group_specific = endo_group_specific,
        endogenous_cov = FALSE,
        lambda_gs = lambda_gs,
        theta_gs = theta_gs,
        S_unbiased = S_unbiased,
        S1 = S1,
        std.lv = std.lv,
        partition = partition,
        userStart = userStart,
        printing = printing,
        s1ori = s1ori,
        s1_type = s1_type
      )

      SM_2 <- Step2_slow(
        ngroups = ngroups,
        nclus = nclus,
        nstarts = nstarts,
        N_gs = N_gs,
        seed = seed,
        max_it = max_it,
        cov_eta = cov_eta,
        dat = dat,
        S2 = S2,
        lat_var = lat_var,
        ordered = ordered,
        endo_group_specific = endo_group_specific,
        endogenous_cov = endogenous_cov,
        lambda_gs = lambda_gs,
        theta_gs = theta_gs,
        S_unbiased = S_unbiased,
        S1 = S1,
        std.lv = std.lv,
        partition = partition,
        userStart = userStart,
        printing = printing,
        s1ori = s1ori,
        only_slow = FALSE,
        beta_ks = SM_1$beta_ks,
        psi_gks = SM_1$psi_gks,
        z_gks = SM_1$z_gks
      )

      return(list(
        iter = SM_2$iter,
        z_gks = SM_2$z_gks,
        LL = SM_2$LL,
        s2out = SM_2$s2out,
        loglik_nstarts = SM_1$loglik_nstarts,
        endog = SM_2$endog,
        exog = SM_2$exog,
        endog1 = SM_2$endog1,
        endog2 = SM_2$endog2,
        beta_ks = SM_2$beta_ks,
        psi_gks = SM_2$psi_gks
      ))
    }

    SM <- Step2_slow(
      ngroups = ngroups,
      nclus = nclus,
      nstarts = nstarts,
      N_gs = N_gs,
      seed = seed,
      max_it = max_it,
      cov_eta = cov_eta,
      dat = dat,
      S2 = S2,
      lat_var = lat_var,
      ordered = ordered,
      endo_group_specific = endo_group_specific,
      endogenous_cov = endogenous_cov,
      lambda_gs = lambda_gs,
      theta_gs = theta_gs,
      S_unbiased = S_unbiased,
      S1 = S1,
      std.lv = std.lv,
      partition = partition,
      userStart = userStart,
      printing = printing,
      s1ori = s1ori,
      only_slow = TRUE
    )

    return(list(
      iter = SM$iter,
      z_gks = SM$z_gks,
      LL = SM$LL,
      s2out = SM$s2out,
      loglik_nstarts = SM$loglik_nstarts,
      endog = SM$endog,
      exog = SM$exog,
      endog1 = SM$endog1,
      endog2 = SM$endog2,
      beta_ks = SM$beta_ks,
      psi_gks = SM$psi_gks
    ))
  }

  SM <- Step2_fast(
    ngroups = ngroups,
    nclus = nclus,
    nstarts = nstarts,
    N_gs = N_gs,
    seed = seed,
    max_it = max_it,
    cov_eta = cov_eta,
    dat = dat,
    S2 = S2,
    lat_var = lat_var,
    ordered = ordered,
    endo_group_specific = endo_group_specific,
    endogenous_cov = endogenous_cov,
    lambda_gs = lambda_gs,
    theta_gs = theta_gs,
    S_unbiased = S_unbiased,
    S1 = S1,
    std.lv = std.lv,
    partition = partition,
    userStart = userStart,
    printing = printing,
    s1ori = s1ori,
    s1_type = s1_type
  )

  return(list(
    iter = SM$iter,
    z_gks = SM$z_gks,
    LL = SM$LL,
    s2out = SM$s2out,
    loglik_nstarts = SM$loglik_nstarts,
    endog = SM$endog,
    exog = SM$exog,
    endog1 = SM$endog1,
    endog2 = SM$endog2,
    beta_ks = SM$beta_ks,
    psi_gks = SM$psi_gks
  ))
}
