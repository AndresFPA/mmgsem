#' @import ggplot2
#' @import lavaan
#' @import semTools
NULL
#> NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "mmgsem is currently in development and may contain bugs. ",
    "Please report issues at ",
    "https://github.com/AndresFPA/mmgsem/issues."
  )
}
