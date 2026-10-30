CHull <- function(loglik, nrpar, nsclust) {
  clusters <- seq.int(nsclust[1], nsclust[2])
  fit_mat <- cbind(clusters, loglik, nrpar)
  nrows <- nrow(fit_mat)
  scores <- rep(NA_real_, nrows)
  if (nrows < 3L) {
    return(cbind(fit_mat, scores))
  }

  calculate_score <- function(left, middle, right) {
    ll_left <- fit_mat[left, 2]
    ll_middle <- fit_mat[middle, 2]
    ll_right <- fit_mat[right, 2]
    par_left <- fit_mat[left, 3]
    par_middle <- fit_mat[middle, 3]
    par_right <- fit_mat[right, 3]
    par_span <- par_right - par_left
    left_span <- par_middle - par_left
    right_span <- par_right - par_middle
    right_gain <- ll_right - ll_middle

    if (any(!is.finite(c(ll_left, ll_middle, ll_right,
                         par_left, par_middle, par_right))) ||
        par_span == 0 || left_span == 0 || right_span == 0 ||
        right_gain == 0) {
      return(NA_real_)
    }

    point_line <- ll_left +
      ((ll_right - ll_left) / par_span) * left_span
    if (ll_middle < point_line - 0.01) {
      return(NA_real_)
    }
    score <- ((ll_middle - ll_left) / left_span) /
      (right_gain / right_span)
    if (is.finite(score)) score else NA_real_
  }

  for (middle in 2:(nrows - 1L)) {
    scores[middle] <- calculate_score(middle - 1L, middle, middle + 1L)
  }

  hull <- c(1L, which(!is.na(scores)), nrows)
  changed <- length(hull) < nrows
  while (length(hull) > 2L && changed) {
    changed <- FALSE
    for (i in 2:(length(hull) - 1L)) {
      left <- hull[i - 1L]
      middle <- hull[i]
      right <- hull[i + 1L]
      if (!identical(c(left, middle, right), c(middle - 1L, middle,
                                               middle + 1L))) {
        score <- calculate_score(left, middle, right)
        if (is.na(score)) {
          scores[middle] <- NA_real_
          changed <- TRUE
        } else {
          scores[middle] <- score
        }
      }
    }
    hull <- c(1L, which(!is.na(scores)), nrows)
  }

  cbind(fit_mat, abs(scores))
}
