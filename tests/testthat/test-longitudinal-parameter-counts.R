test_that("longitudinal parameter counts include group-time covariances", {
  s2_parameters <- lavaan::lavaanify(
    "y ~ m + x\nm ~ x",
    auto = TRUE
  )

  counts <- mmgsem:::count_lmmgsem_parameters(
    s2_parameters = s2_parameters,
    ngroups = 2L,
    ntimes = 3L,
    nstates = 4L,
    n_regressions = 3L,
    n_measurement_parameters = 17L
  )

  n_cells <- 2L * 3L
  n_exogenous_covariances <- 1L * n_cells
  n_endog1_variances <- 1L * n_cells
  n_endog2_covariances <- 1L * n_cells
  n_regression_parameters <- 3L * 4L
  n_mixture_parameters <- (4L - 1L) + 4L * (4L - 1L)
  expected_structural <- n_exogenous_covariances +
    n_endog1_variances + n_endog2_covariances +
    n_regression_parameters + n_mixture_parameters

  expect_equal(counts$structural, expected_structural)
  expect_equal(counts$observed, 17L + expected_structural)
})
