# lavaan keeps these helpers internal; preserve their behavior without ::: calls.
lavaan_loglik_samplestats <- function(...) {
  loglik_function <- utils::getFromNamespace(
    "lav_mvnorm_loglik_samplestats",
    "lavaan"
  )
  loglik_function(...)
}

lavaan_marker_indices <- function(lambda) {
  marker_function <- utils::getFromNamespace(
    "lav_utils_get_marker",
    "lavaan"
  )
  marker_function(lambda)
}
