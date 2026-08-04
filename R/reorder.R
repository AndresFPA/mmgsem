#' Internal helper: reorder latent variable matrices for MMGSEM.
#'
#' Reorders latent variable matrices so exogenous factors appear before
#' endogenous factors, matching the expected internal SEM ordering.
#'
#' @param x matrix to reorder.
#' @param exog character vector of exogenous latent variable names.
#' @param endog character vector of endogenous latent variable names.
#' @return Reordered matrix.
#' @keywords internal
reorder <- function(x, exog, endog) {
  x[c(exog, endog), c(exog, endog)]
}

#' Internal helper: reorder measurement model matrices for MMGSEM.
#'
#' Reorders Lambda and Theta matrices so the observed variable ordering
#' matches the structural model ordering used in Step 2.
#'
#' @param x matrix to reorder.
#' @param matrix character string, either "lambda" or "theta".
#' @param exog exogenous latent variable names.
#' @param endog endogenous latent variable names.
#' @param endog1 endogenous variables that are also predictors.
#' @param endog2 purely dependent endogenous variables.
#' @param S1 step 1 model syntax.
#' @param dat original data frame.
#' @return Reordered matrix.
#' @keywords internal
reorder_obs <- function(x, matrix, exog, endog, endog1, endog2, S1, dat) {
  lines_model <- unlist(strsplit(unlist(S1), "\n"))
  rewritten <- character(0)

  for (lbl in exog) {
    rewritten <- c(rewritten, lines_model[grepl(lbl, lines_model)])
  }
  if (length(endog1) > 0) {
    for (lbl in endog1) {
      rewritten <- c(rewritten, lines_model[grepl(lbl, lines_model)])
    }
  }
  for (lbl in endog2) {
    rewritten <- c(rewritten, lines_model[grepl(lbl, lines_model)])
  }

  fake_measur <- lavaan::cfa(model = rewritten, data = dat, do.fit = FALSE)
  correct_vars <- lavaan::lavNames(fake_measur)

  if (matrix == "lambda") {
    x[correct_vars, c(exog, endog)]
  } else {
    x[correct_vars, correct_vars]
  }
}
