#' Longitudinal mixture multigroup structural equation modeling
#'
#' Fits a stepwise measurement-invariance model and then an HMM whose latent
#' states represent clusters of groups over time. Each group must be observed
#' at every time point, with at least two time points.
#'
#' @param dat Data frame containing the observed variables, group identifier,
#'   and time identifier.
#' @param S1 Measurement model in lavaan syntax. Required when `s1_fit` is not
#'   supplied.
#' @param S2 Structural model in lavaan syntax.
#' @param group Name of the grouping variable in `dat`.
#' @param time Name of the time variable in `dat`.
#' @param nclus Number of latent states (clusters).
#' @param seed Optional random seed for reproducibility.
#' @param max_it Maximum number of EM iterations per random start.
#' @param nstarts Number of random starts.
#' @param printing If `TRUE`, print iteration progress.
#' @param s2_fast If `TRUE`, use the fast Step 2 estimator; otherwise use the
#'   slow estimator.
#' @param s1_fit Optional lavaan measurement-model fit with one group per
#'   group-time cell. Its group labels must match
#'   `interaction(dat[[group]], dat[[time]], sep = ".Time")`. If omitted,
#'   measurement invariance is estimated using
#'   [compute_mi_long_stepwise_score()].
#' @param mi_method Method for selecting invariance violations:
#'   `"group_violations"` or `"test_statistic"`.
#' @param mi_alpha Significance level for the score tests.
#' @param mi_d_cfi_threshold Maximum configural-to-metric CFI difference.
#' @param mi_max_iter Maximum iterations in each stepwise invariance analysis.
#' @param mi_fit_indices Fit indices returned by the invariance analysis.
#' @param mi_specific_noninv If `TRUE`, allow non-invariance for specific
#'   group-time pairs; otherwise free loadings globally.
#' @param ... Additional arguments passed to lavaan's CFA fitting functions
#'   during the measurement-invariance analysis.
#'
#' @returns A `lmmgsem` object containing time-specific posterior state
#'   probabilities, transition probabilities, model parameters, likelihoods,
#'   and model-selection criteria. Information criteria use the sequence-level
#'   factor likelihood and marginal posterior state entropy. Standard errors
#'   are not currently implemented. Parameter counts include free measurement
#'   parameters, per-state regression parameters, and initial/transition
#'   probabilities.
#' @export
lmmgsem <- function(dat, S1 = NULL, S2, group, time, nclus, seed = NULL,
                    max_it = 10000L, nstarts = 20L, printing = FALSE,
                    s2_fast = TRUE, s1_fit = NULL,
                    mi_method = c("group_violations", "test_statistic"),
                    mi_alpha = 0.05, mi_d_cfi_threshold = 0.01,
                    mi_max_iter = 20L,
                    mi_fit_indices = c("chisq", "df", "pvalue", "cfi",
                                       "rmsea", "srmr"),
                    mi_specific_noninv = FALSE, ...) {
  if (!is.data.frame(dat)) {
    stop("'dat' must be a data frame.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1L ||
      !group %in% names(dat)) {
    stop("'group' must name a column in 'dat'.", call. = FALSE)
  }
  if (!is.character(time) || length(time) != 1L ||
      !time %in% names(dat)) {
    stop("'time' must name a column in 'dat'.", call. = FALSE)
  }
  if (identical(group, time)) {
    stop("'group' and 'time' must refer to different columns.", call. = FALSE)
  }
  if (!is.null(S1) &&
      (!is.character(S1) || length(S1) != 1L || !nzchar(S1))) {
    stop("'S1' must be NULL or a single non-empty lavaan model string.",
         call. = FALSE)
  }
  if (is.null(s1_fit) && is.null(S1)) {
    stop("Supply 'S1' or a fitted lavaan object through 's1_fit'.",
         call. = FALSE)
  }
  if (!is.character(S2) || length(S2) != 1L || !nzchar(S2)) {
    stop("'S2' must be a single non-empty lavaan model string.",
         call. = FALSE)
  }
  if (length(nclus) != 1L || !is.numeric(nclus) ||
      !is.finite(nclus) || nclus < 1 || nclus != as.integer(nclus)) {
    stop("'nclus' must be a positive integer.", call. = FALSE)
  }
  if (length(max_it) != 1L || !is.numeric(max_it) ||
      !is.finite(max_it) || max_it < 1 || max_it != as.integer(max_it)) {
    stop("'max_it' must be a positive integer.", call. = FALSE)
  }
  if (length(nstarts) != 1L || !is.numeric(nstarts) ||
      !is.finite(nstarts) || nstarts < 1 || nstarts != as.integer(nstarts)) {
    stop("'nstarts' must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(s2_fast) || length(s2_fast) != 1L || is.na(s2_fast)) {
    stop("'s2_fast' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(printing) || length(printing) != 1L || is.na(printing)) {
    stop("'printing' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
       seed != as.integer(seed))) {
    stop("'seed' must be NULL or one finite integer.", call. = FALSE)
  }
  if (anyNA(dat[[group]]) || anyNA(dat[[time]])) {
    stop("The group and time identifiers cannot contain missing values.",
         call. = FALSE)
  }

  mi_method <- match.arg(mi_method)
  group_levels <- unique(as.character(dat[[group]]))
  time_levels <- unique(as.character(dat[[time]]))
  if (length(time_levels) < 2L) {
    stop("Longitudinal estimation requires at least two time points.",
         call. = FALSE)
  }
  if (nclus > length(group_levels)) {
    stop("'nclus' cannot exceed the number of groups.", call. = FALSE)
  }
  if (length(mi_alpha) != 1L || !is.numeric(mi_alpha) ||
      !is.finite(mi_alpha) || mi_alpha <= 0 || mi_alpha >= 1) {
    stop("'mi_alpha' must be strictly between 0 and 1.", call. = FALSE)
  }
  if (length(mi_d_cfi_threshold) != 1L ||
      !is.numeric(mi_d_cfi_threshold) ||
      !is.finite(mi_d_cfi_threshold) || mi_d_cfi_threshold < 0) {
    stop("'mi_d_cfi_threshold' must be a non-negative finite number.",
         call. = FALSE)
  }
  if (length(mi_max_iter) != 1L || !is.numeric(mi_max_iter) ||
      !is.finite(mi_max_iter) || mi_max_iter < 1 ||
      mi_max_iter != as.integer(mi_max_iter)) {
    stop("'mi_max_iter' must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(mi_specific_noninv) || length(mi_specific_noninv) != 1L ||
      is.na(mi_specific_noninv)) {
    stop("'mi_specific_noninv' must be TRUE or FALSE.", call. = FALSE)
  }

  dat[[group]] <- factor(as.character(dat[[group]]), levels = group_levels)
  dat[[time]] <- factor(as.character(dat[[time]]), levels = time_levels)
  cell_counts <- table(dat[[group]], dat[[time]])
  if (any(cell_counts == 0L)) {
    stop("Every group must have observations at every time point.",
         call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  mi_result <- NULL
  if (is.null(s1_fit)) {
    mi_result <- compute_mi_long_stepwise_score(
      data = dat,
      model = S1,
      group_var = group,
      time_var = time,
      method = mi_method,
      alpha = mi_alpha,
      d_cfi_threshold = mi_d_cfi_threshold,
      max_iter = mi_max_iter,
      fit_indices = mi_fit_indices,
      specific_noninv = mi_specific_noninv,
      ...
    )
    s1_fit <- mi_result$final_fit
  }
  if (!inherits(s1_fit, "lavaan")) {
    stop("'s1_fit' must be a fitted lavaan model object.", call. = FALSE)
  }
  if (!isTRUE(lavaan::lavInspect(s1_fit, "converged"))) {
    stop("The Step 1 lavaan model did not converge.",
         call. = FALSE)
  }

  expected_cells <- levels(interaction(
    dat[[group]], dat[[time]], sep = ".Time", drop = TRUE
  ))
  fit_cells <- as.character(lavaan::lavInspect(s1_fit, "group.label"))
  if (length(fit_cells) != length(expected_cells) ||
      anyDuplicated(fit_cells) > 0L ||
      !setequal(fit_cells, expected_cells)) {
    stop(
      "'s1_fit' must be a lavaan multigroup model fitted to every ",
      "group-time cell. Its group labels must match ",
      "interaction(dat[[group]], dat[[time]], sep = '.Time').",
      call. = FALSE
    )
  }

  cov_eta <- lavaan::lavInspect(s1_fit, "cov.lv")
  if (is.array(cov_eta) && length(dim(cov_eta)) == 3L) {
    cov_eta <- lapply(seq_len(dim(cov_eta)[3]), function(i) cov_eta[, , i])
  }
  if (is.matrix(cov_eta) && length(expected_cells) == 1L) {
    cov_eta <- list(cov_eta)
  }
  if (!is.list(cov_eta) || length(cov_eta) != length(expected_cells)) {
    stop("The invariance fit must provide one latent covariance matrix for ",
         "each group-time cell.", call. = FALSE)
  }
  cov_eta <- cov_eta[match(expected_cells, fit_cells)]
  if (!all(vapply(cov_eta, is.matrix, logical(1))) ||
      !all(vapply(cov_eta, function(x) all(is.finite(x)), logical(1)))) {
    stop("The invariance fit contains invalid latent covariance matrices.",
         call. = FALSE)
  }
  covariance_dimensions <- lapply(cov_eta, dim)
  if (!all(vapply(covariance_dimensions, identical, logical(1),
                  covariance_dimensions[[1]])) ||
      length(covariance_dimensions[[1]]) != 2L ||
      covariance_dimensions[[1]][1] != covariance_dimensions[[1]][2]) {
    stop("All latent covariance matrices must be square and have the same ",
         "dimensions.", call. = FALSE)
  }

  nobs <- unlist(lavaan::lavInspect(s1_fit, "nobs"), use.names = TRUE)
  nobs_order <- if (!is.null(names(nobs)) &&
                    setequal(names(nobs), fit_cells)) {
    names(nobs)
  } else {
    fit_cells
  }
  N_gs <- as.numeric(nobs[match(expected_cells, nobs_order)])
  if (length(N_gs) != length(expected_cells) || any(!is.finite(N_gs)) ||
      any(N_gs <= 0)) {
    stop("The measurement model must have positive sample sizes in every ",
         "group-time cell.", call. = FALSE)
  }

  s2_parameters <- lavaan::lavaanify(S2, auto = TRUE)
  regression_parameters <- s2_parameters$free[
    s2_parameters$op == "~" & s2_parameters$free > 0
  ]
  n_regressions <- length(unique(regression_parameters))
  if (n_regressions == 0L) {
    stop("'S2' must contain at least one regression.", call. = FALSE)
  }

  ngroups <- length(group_levels)
  ntimes <- length(time_levels)
  hmm_fit <- hmm_mmgsem_long(
    cov_eta = cov_eta,
    S2 = S2,
    ngroups = ngroups,
    nstates = as.integer(nclus),
    ntimes = ntimes,
    N_gs = N_gs,
    max_it = as.integer(max_it),
    s2_fast = s2_fast,
    nstarts = as.integer(nstarts),
    printing = printing
  )

  post <- hmm_fit$post
  dimnames(post) <- list(group_levels, as.character(seq_len(nclus)), time_levels)
  modal_post <- matrix(
    0L,
    nrow = ngroups,
    ncol = ntimes,
    dimnames = list(group_levels, time_levels)
  )
  for (g in seq_len(ngroups)) {
    for (t in seq_len(ntimes)) {
      modal_post[g, t] <- which.max(post[g, , t])
    }
  }

  posterior_entropy <- -sum(post * log(pmax(post, .Machine$double.xmin)))
  state_proportions <- apply(post, 2L, mean)
  state_entropy <- -sum(state_proportions *
                          log(pmax(state_proportions, .Machine$double.xmin)))
  posterior_entropy_mean <- posterior_entropy / (ngroups * ntimes)
  entropy_r2 <- if (state_entropy == 0) {
    0
  } else {
    max(0, min(1, (state_entropy - posterior_entropy_mean) / state_entropy))
  }

  measurement_parameters <- lavaan::parTable(s1_fit)$free
  n_measurement_parameters <- length(unique(
    measurement_parameters[measurement_parameters > 0]
  ))
  n_parameters <- n_measurement_parameters +
    as.integer(nclus) * n_regressions + as.integer(nclus)^2 - 1L
  loglik <- as.numeric(hmm_fit$logLik)
  total_n <- sum(N_gs)
  bic_n <- -2 * loglik + n_parameters * log(total_n)
  bic_g <- -2 * loglik + n_parameters * log(ngroups)
  aic <- -2 * loglik + 2 * n_parameters
  aic3 <- -2 * loglik + 3 * n_parameters
  icl <- bic_g + 2 * posterior_entropy

  output <- list(
    posteriors = post,
    modal_post = modal_post,
    transition = hmm_fit$A,
    initial_probabilities = hmm_fit$pi,
    final_fit = hmm_fit$s2out,
    MM = s1_fit,
    MI = mi_result,
    param = list(
      beta_ks = hmm_fit$param$beta_ks,
      psi_gks = hmm_fit$param$psi_gks,
      cov_eta = cov_eta
    ),
    logLik = list(
      loglik = loglik,
      obs_loglik = loglik,
      runs_loglik = unlist(hmm_fit$prev_starts$logliks)
    ),
    model_sel = list(
      BIC = list(
        observed = list(BIC_N = bic_n, BIC_G = bic_g),
        Factors = list(BIC_N = bic_n, BIC_G = bic_g)
      ),
      AIC = list(observed = aic, Factors = aic),
      AIC3 = list(observed = aic3, Factors = aic3),
      ICL = list(observed = icl, Factors = icl),
      R2_entropy = entropy_r2
    ),
    NrPar = list(Obs.nrpar = n_parameters, Fac.nrpar = n_parameters),
    N_gs = N_gs,
    ngroups = ngroups,
    ntimes = ntimes,
    nclus = as.integer(nclus),
    nstarts = as.integer(nstarts),
    iterations = hmm_fit$iterations,
    converged = hmm_fit$converged,
    group_levels = group_levels,
    time_levels = time_levels
  )
  class(output) <- c("lmmgsem", "mmgsem")
  output
}
