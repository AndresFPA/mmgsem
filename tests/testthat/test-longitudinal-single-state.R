test_that("longitudinal HMM supports a single state", {
  cov_eta <- replicate(
    4,
    matrix(
      c(1, 0.3, 0.3, 1),
      nrow = 2,
      dimnames = list(c("x", "y"), c("x", "y"))
    ),
    simplify = FALSE
  )

  for (s2_fast in c(TRUE, FALSE)) {
    fit <- mmgsem:::hmm_mmgsem_long(
      cov_eta = cov_eta,
      S2 = "y ~ x",
      ngroups = 2L,
      nstates = 1L,
      ntimes = 2L,
      N_gs = rep(20, 4),
      max_it = 1L,
      s2_fast = s2_fast,
      nstarts = 1L
    )

    expect_identical(dim(fit$post), c(2L, 1L, 2L))
    expect_equal(unname(apply(fit$post, c(1L, 3L), sum)), matrix(1, 2, 2))
    expect_true(is.finite(fit$logLik))
  }
})
