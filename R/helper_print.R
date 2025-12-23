#' @export
print.RSTr <- function(x, ...) {
  cat("RSTr object:\n\n")
  cat("Model name:", x$params$name, "\n")
  cat("Model type:", toupper(x$params$model), "\n")
  cat("Data likelihood:", x$params$method, "\n")
  ci_str <- paste0(round(x$params$perc_ci * 100, 2), "%")
  cat("Estimate Credible Interval:", ci_str, "\n")
  cat("Number of geographic units:", dim(x$data$Y)[1], "\n")
  cat("Number of samples:", x$params$total, "\n")
  cat(
    "Estimates age-standardized:",
    ifelse(x$params$age_standardized, "Yes", "No"),
    "\n"
  )
  if (x$params$age_standardized) {
    cat("Age-standardized groups:", colnames(x$medians_as), "\n")
  }
  cat("Estimates suppressed:", ifelse(x$params$suppressed, "Yes", "No"), "\n")
  if (x$params$suppressed) {
    if (x$params$age_standardized) {
      n_rel_str <- "Number of reliable age-standardized rates:"
      tot_reliable <- sum(!is.na(x$medians_suppressed_as))
      tot_estimate <- length(x$medians_suppressed_as)
    } else {
      n_rel_str <- "Number of reliable rates:"
      tot_reliable <- sum(!is.na(x$medians_suppressed))
      tot_estimate <- length(x$medians_suppressed)
    }
    pct <- round(tot_reliable / tot_estimate * 100, 1)
    cat(n_rel_str, tot_reliable, "/", tot_estimate, paste0("(", pct, "%)\n"))
  }
}
