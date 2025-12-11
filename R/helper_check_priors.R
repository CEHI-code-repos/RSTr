check_priors <- function(RSTr_obj) {
  UseMethod("check_priors")
}

#' @export
check_priors.ucar <- function(RSTr_obj) {
  priors <- RSTr_obj$priors
  chk <- c("tau_a", "tau_b", "lambda_sd", "lambda_accept", "sig_a", "sig_b")
  check_missing_priors(priors, chk)
  # Check for warnings
  check_unused_priors(priors, chk)
  # Check for errors
  c(
    check_tau_a(priors$tau_a),
    check_tau_b(priors$tau_b),
    check_lambda_sd(priors$lambda_sd, RSTr_obj$data$Y),
    check_sig_a(priors$sig_a),
    check_sig_b(priors$sig_b)
  )
}

#' @export
check_priors.mcar <- function(RSTr_obj) {
  priors <- RSTr_obj$priors
  num_group <- dim(RSTr_obj$data$Y)[2]
  chk <- c("tau_a", "tau_b", "lambda_sd", "lambda_accept", "G_scale", "G_df")
  check_missing_priors(priors, chk)
  # Check for warnings
  check_unused_priors(priors, chk)
  # Check for errors
  c(
    check_tau_a(priors$tau_a),
    check_tau_b(priors$tau_b),
    check_lambda_sd(priors$lambda_sd, RSTr_obj$data$Y),
    check_G_scale(priors$G_scale, num_group),
    check_G_df(priors$G_df, num_group)
  )
}

#' @export
check_priors.mstcar <- function(RSTr_obj) {
  priors <- RSTr_obj$priors
  num_group <- dim(RSTr_obj$data$Y)[2]
  chk <- c("tau_a", "tau_b", "lambda_sd", "lambda_accept", "G_scale", "G_df", "Ag_scale", "Ag_df")
  check_missing_priors(priors, chk)
  # Check for warnings
  check_unused_priors(priors, chk)
  # Check for errors
  c(
    check_tau_a(priors$tau_a),
    check_tau_b(priors$tau_b),
    check_lambda_sd(priors$lambda_sd, RSTr_obj$data$Y),
    check_G_scale(priors$G_scale, num_group),
    check_G_df(priors$G_df, num_group),
    check_Ag_scale(priors$Ag_scale, num_group),
    check_Ag_df(priors$Ag_df, num_group)
  )
}

#' @export
check_priors.mstcar_update_rho <- function(RSTr_obj, errout) {
  priors <- RSTr_obj$priors
  num_group <- dim(RSTr_obj$data$Y)[2]
  chk <- c("tau_a", "tau_b", "lambda_sd", "lambda_accept", "G_scale", "G_df", "Ag_scale", "Ag_df", "rho_a", "rho_b", "rho_sd", "rho_accept")
  check_missing_priors(priors, chk)
  # Check for warnings
  check_unused_priors(priors, chk)
  # Check for errors
  c(
    check_tau_a(priors$tau_a),
    check_tau_b(priors$tau_b),
    check_lambda_sd(priors$lambda_sd, RSTr_obj$data$Y),
    check_G_scale(priors$G_scale, num_group),
    check_G_df(priors$G_df, num_group),
    check_Ag_scale(priors$Ag_scale, num_group),
    check_Ag_df(priors$Ag_df, num_group),
    check_rho_a(priors$rho_a),
    check_rho_b(priors$rho_b),
    check_rho_sd(priors$rho_sd, num_group)
  )
}

check_missing_priors <- function(priors, chk) {
  miss <- sapply(chk, \(x) !any(names(priors) == x))
  if (any(miss)) {
    stop("One or more objects missing from list 'priors': ", toString(chk[miss]))
  }
}

check_unused_priors <- function(priors, chk) {
  chk_elem <- !(names(priors) %in% chk)
  if (any(chk_elem)) {
    warning("Unused elements of list 'priors':", toString(names(priors)[chk_elem]))
  }
}

check_tau_a <- function(tau_a) {
  # is non-positive or infinite
  if ((tau_a <= 0) || !is.finite(tau_a)) {
    "tau_a is not positive. Ensure tau_a > 0 and not infinite or use default value"
  }
}

check_tau_b <- function(tau_b) {
  # is non-positive or infinite
  if ((tau_b <= 0) || !is.finite(tau_b)) {
    "tau_b is not positive. Ensure tau_b > 0 and not infinite or use default value"
  }
}

check_lambda_sd <- function(lambda_sd, Y) {
  err_messages <- character()
  # dim not num_time num_region
  if (!all(dim(lambda_sd) == dim(Y))) {
    err_messages <- c(
      err_messages,
      "lambda_sd has different length than data. Ensure length(lambda_sd) == length(Y) or use default value"
    )
  }
  # is non-positive or infinite
  if (any((lambda_sd <= 0) | !is.finite(lambda_sd))) {
    err_messages <- c(
      err_messages,
      "lambda_sd contains non-positive values. Ensure all(lambda_sd > 0) and not infinite or use default value"
    )
  }
  err_messages
}

check_sig_a <- function(sig_a) {
  # is non-positive or infinite
  if ((sig_a <= 0) || !is.finite(sig_a)) {
    "sig_a is not positive. Ensure sig_a > 0 and not infinite or use default value"
  }
}

check_sig_b <- function(sig_b) {
  # is non-positive or infinite
  if ((sig_b <= 0) || !is.finite(sig_b)) {
    "sig_b is not positive. Ensure sig_b > 0 and not infinite or use default value"
  }
}

check_G_scale <- function(G_scale, num_group) {
  err_messages <- character()
  # dimensions don't match num_group num_group
  if (!all(dim(G_scale) == c(num_group, num_group))) {
    err_messages <- c(
      err_messages,
      "G_scale is not an num_group x num_group matrix. Ensure dim(G_scale) == num_group x num_group or use default value"
    )
  }
  # matrix is not symmetric
  if (!isSymmetric(G_scale)) {
    err_messages <- c(
      err_messages,
      "G_scale is not symmetric. Ensure G_scale is symmetric or use default value"
    )
  }
  # values are infinite
  if (!all(is.finite(G_scale))) {
    err_messages <- c(
      err_messages,
      "G_scale contains infinite values. Ensure G_scale is finite or use default value"
    )
  }
  # diagonals are not positive
  if (any(diag(G_scale) <= 0)) {
    err_messages <- c(
      err_messages,
      "diag(G_scale) contains non-positive values. Ensure diag(G_scale) is all positive or use default value"
    )
  }
  err_messages
}

check_G_df <- function(G_df, num_group) {
  err_messages <- character()
  # is not a whole number
  if (floor(G_df) != G_df) {
    err_messages <- c(
      err_messages,
      "G_df is not a whole number. Ensure G_df is whole number or use default value"
    )
  }
  # is less than df
  if (G_df <= num_group - 1) {
    err_messages <- c(
      err_messages,
      "G_df too small. Ensure G_df > num_group - 1 or use default value"
    )
  }
  err_messages
}

check_Ag_scale <- function(Ag_scale, num_group) {
  err_messages <- character()
  # Ag_scale
  # dimensions don't match num_group num_group
  if (!all(dim(Ag_scale) == c(num_group, num_group))) {
    err_messages <- c(
      err_messages,
      "Ag_scale is not an num_group x num_group matrix. Ensure dim(Ag_scale) == num_group x num_group or use default value"
    )
  }
  # matrix is not symmetric
  if (!isSymmetric(Ag_scale)) {
    err_messages <- c(
      err_messages,
      "Ag_scale is not symmetric. Ensure Ag_scale is symmetric or use default value"
    )
  }
  # values are infinite
  if (!all(is.finite(Ag_scale))) {
    err_messages <- c(
      err_messages,
      "Ag_scale contains infinite values. Ensure Ag_scale is finite or use default value"
    )
  }
  # diagonals are not positive
  if (any(diag(Ag_scale) <= 0)) {
    err_messages <- c(
      err_messages,
      "diag(Ag_scale) contains non-positive values. Ensure diag(Ag_scale) is positive or use default value"
    )
  }
  err_messages
}

check_Ag_df <- function(Ag_df, num_group) {
  err_messages <- character()
  # is not a whole number
  if (floor(Ag_df) != Ag_df) {
    err_messages <- c(
      err_messages,
      "Ag_df is not a whole number. Ensure Ag_df is whole number or use default value"
    )
  }
  # is less than df
  if (Ag_df <= num_group - 1) {
    err_messages <- c(
      err_messages,
      "Ag_df too small. Ensure Ag_df > num_group - 1 or use default value"
    )
  }
  err_messages
}

check_rho_a <- function(rho_a) {
  # is non-positive or infinite
  if ((rho_a <= 0) || !is.finite(rho_a)) {
    "rho_a is not positive. Ensure rho_a > 0 and not infinite or use default value"
  }
}

check_rho_b <- function(rho_b) {
  # is non-positive or infinite
  if ((rho_b <= 0) || !is.finite(rho_b)) {
    "rho_b is not positive. Ensure rho_b > 0 and not infinite or use default value"
  }
}

check_rho_sd <- function(rho_sd, num_group) {
  err_messages <- character()
  # length not num_group
  if (length(rho_sd) != num_group) {
    err_messages <- c(
      err_messages,
      "rho_sd is not length num_group. Ensure length(rho_sd) == num_group or use default value"
    )
  }
  # is non-positive or infinite
  if (any((rho_sd <= 0) | !is.finite(rho_sd))) {
    err_messages <- c(
      err_messages,
      "rho_sd contains non-positive values. Ensure all(rho_sd > 0) and not infinite or use default value"
    )
  }
  err_messages
}