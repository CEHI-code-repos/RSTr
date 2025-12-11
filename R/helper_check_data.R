check_data <- function(RSTr_obj, errout = NULL) {
  data <- RSTr_obj$data
  check_missing_data_objects(data)
  # Check for warnings
  check_unused_data(data)
  # Check for errors
  c(
    errout,
    check_mismatched_dimensions(data),
    check_invalid_Y(data$Y),
    check_invalid_n(data$n),
    check_zero_events(RSTr_obj)
  )
}

check_unused_data <- function(data) {
  chk_elem <- which(!(names(data) %in% c("Y", "n")))
  if (length(chk_elem)) {
    warning("Unused elements of list 'data':", toString(names(data)[chk_elem]))
  }
}

check_missing_data_objects <- function(data) {
  chk <- c("Y", "n")
  miss <- sapply(chk, \(x) !any(names(data) == x))
  if (any(miss)) {
    stop("One or more objects missing from list 'data': ", toString(chk[miss]))
  }
}

check_mismatched_dimensions <- function(data) {
  if (any(dim(data$Y) != dim(data$n))) {
    "Data not same dimensions. Ensure dim(Y) == dim(n)"
  }
}

check_invalid_Y <- function(Y) {
  Ychk <- Y[(!is.na(Y)) & (!is.null(Y))]
  if (any((Ychk < 0) | is.infinite(Ychk))) {
    "Invalid Y values. Check that all Y's are at least 0 and finite"
  }
}

check_invalid_n <- function(n) {
  if (any((n < 0) | is.infinite(n))) {
    "Invalid n values. Check that all n's are at least 0 and finite"
  }
}

check_zero_events <- function(RSTr_obj) {
  UseMethod("check_zero_events")
}

#' @export
check_zero_events.ucar <- function(RSTr_obj) {
  if (any(apply(RSTr_obj$data$Y, 2:3, sum) == 0)) {
    "At least one set of regions has no events. Ensure that Y has at least one event for each set of regions"
  }
}

#' @export
check_zero_events.mcar <- function(RSTr_obj) {
  if (any(apply(RSTr_obj$data$Y, 3, sum) == 0)) {
    "No events in Y for at least one time period. Ensure that Y has at least one event for each time period"
  }
}

#' @export
check_zero_events.mstcar <- function(RSTr_obj) {
  if (sum(RSTr_obj$data$Y) == 0) {
    "No events in Y. Ensure that Y has at least one event"
  }
}
