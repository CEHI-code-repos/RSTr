#' Check data
#'
#' @noRd
check_data <- function(data) {
  message("Checking data...")
  Y <- data$Y
  n <- data$n
  chk <- c("Y", "n")
  miss <- sapply(1:length(chk), \(x) !any(names(data) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'data': ", paste(chk[miss], collapse = ", "))
  }

  # Check for warnings
  warnout <- character()
  # Check for unused elements in 'data'
  chk_elem <- which(!(names(data) %in% c("Y", "n")))
  if (length(chk_elem)) {
    errout <- c(errout, paste("Unused elements of list 'data':", paste(names(data)[chk_elem], collapse = ", ")))
  }
  if (length(warnout) != 0) {
    warnout <- paste(seq_along(warnout), ":", warnout)
    warning(paste(length(warnout), "warning(s) found in list 'data':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout <- character()
  # Dimensions of Y and n are not the same
  dimtest <- NULL
  if (is.null(dim(Y))) {
    dimtest <- length(Y) != length(n)
  } else {
    dimtest <- any(dim(Y) != dim(n))
  }
  if (dimtest) {
    errout <- c(errout, "Data not same dimensions. Ensure dim(Y) == dim(n)")
  }
  # Values of Y are either negative or infinite
  Ychk <- Y[which(!is.na(Y) & !is.null(Y))]
  if (any((Ychk < 0) | is.infinite(Ychk))) {
    errout <- c(errout, "Invalid Y values. Check that all Y's are at least 0 and finite")
  }
  # Sum of Ys is zero
  if (sum(Y) == 0) {
    errout <- c(errout, "No events in Y. Ensure that Y has at least one event")
  }
  # Values of n are either negative or infinite
  if (any((n < 0) | is.infinite(n))) {
    errout <- c(errout, "Invalid n values. Check that all n's are at least 0 and finite")
  }
  if (length(errout) != 0) {
    errout <- paste(seq_along(errout), ":", errout)
    stop(paste(length(errout), "error(s) found in list 'data':\n", paste(errout, collapse = "\n ")))
  }
}
