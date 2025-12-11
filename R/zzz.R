# for compatibility with R <4.4.0
if (!exists("%||%", envir = baseenv())) {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
}