check_initial_values <- function(RSTr_obj) {
  UseMethod("check_initial_values")
}

#' @export
check_initial_values.ucar <- function(RSTr_obj) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- RSTr_obj$spatial_data$num_island
  chk <- c("beta", "tau2", "lambda", "Z", "sig2")
  check_missing_initial_values(RSTr_obj, chk)
  # Check for warnings
  check_unused_initial_values(RSTr_obj, chk)
  # Check for errors
  c(
    check_beta(RSTr_obj$initial_values$beta, num_island, num_group, num_time),
    check_lambda(RSTr_obj$initial_values$lambda, Y, method),
    check_sig2(RSTr_obj$initial_values$sig2),
    check_tau2(RSTr_obj$initial_values$tau2),
    check_Z(RSTr_obj$initial_values$Z, Y)
  )
}

#' @export
check_initial_values.mcar <- function(RSTr_obj) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- RSTr_obj$spatial_data$num_island
  chk <- c("beta", "tau2", "lambda", "Z", "G")
  check_missing_initial_values(RSTr_obj, chk)
  # Check for warnings
  check_unused_initial_values(RSTr_obj, chk)

  # Check for errors
  c(
    check_beta(RSTr_obj$initial_values$beta, num_island, num_group, num_time),
    check_lambda(RSTr_obj$initial_values$lambda, Y, method),
    check_G(RSTr_obj$initial_values$G),
    check_tau2(RSTr_obj$initial_values$tau2),
    check_Z(RSTr_obj$initial_values$Z, Y)
  )
}

#' @export
check_initial_values.mstcar <- function(RSTr_obj) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- RSTr_obj$spatial_data$num_island
  chk <- c("beta", "tau2", "lambda", "Z", "G", "rho", "Ag")
  check_missing_initial_values(RSTr_obj, chk)
  # Check for warnings
  check_unused_initial_values(RSTr_obj, chk)
  # Check for errors
  c(
    check_beta(RSTr_obj$initial_values$beta, num_island, num_group, num_time),
    check_lambda(RSTr_obj$initial_values$lambda, Y, method),
    check_G(RSTr_obj$initial_values$G),
    check_rho(RSTr_obj$initial_values$rho),
    check_tau2(RSTr_obj$initial_values$tau2),
    check_Z(RSTr_obj$initial_values$Z, Y),
    check_Ag(RSTr_obj$initial_values$Ag)
  )
}

check_missing_initial_values <- function(RSTr_obj, chk) {
  miss <- sapply(seq_along(chk), \(x) !any(names(RSTr_obj$initial_values) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'initial_values': ", toString(chk[miss]))
  }
}

check_unused_initial_values <- function(RSTr_obj, chk) {
  chk_elem <- which(!(names(RSTr_obj$initial_values) %in% chk))
  if (length(chk_elem)) {
    warning(paste("Unused elements of list 'initial_values':", toString(names(RSTr_obj$initial_values)[chk_elem])))
  }
}

check_beta <- function(beta, num_island, num_group, num_time) {
  err_messages <- character()
  # dimensions don't match num_island num_group num_time
  if (!all(dim(beta) == c(num_island, num_group, num_time))) {
    err_messages <- c(
      err_messages,
      "beta is not an num_island x num_group x num_time array. Ensure dim(beta) == num_island x num_group x num_time or use default value"
    )
  }
  # values are infinite
  if (!all(is.finite(beta))) {
    err_messages <- c(
      err_messages,
      "beta contains infinite values. Ensure all(is.finite(beta)) or use default value"
    )
  }
  err_messages
}

check_lambda <- function(lambda, Y, method) {
  err_messages <- character()
  # dimensions don't match num_region num_group num_time
  if (!all(dim(lambda) == dim(Y))) {
    err_messages <- c(
      err_messages,
      "lambda is not a num_region x num_group x num_time array. Ensure dim(lambda) == dim(Y) or use default value"
    )
  }
  # values are unsupported
  lower_lim <- 0
  upper_lim <- ifelse(method == "binomial", 1, Inf)
  if (any((lambda <= lower_lim) | (lambda >= upper_lim))) {
    err_messages <- c(
      err_messages,
      "lambda contains unsupported values. Ensure lambdas are within range (0, 1) for `method = binomial` or (0, Inf) for `method = poisson` or use default value"
    )
  }
  err_messages
}

check_sig2 <- function(sig2) {
  # is non-positive or infinite
  if (any(sig2 <= 0) || !all(is.finite(sig2))) {
      "sig2 contains non-positive or infinite values. Ensure all sig2 > 0 and not infinite or use default value"
  }
}

check_tau2 <- function(tau2) {
  # is non-positive or infinite
  if (any(tau2 <= 0) || !all(is.finite(tau2))) {
    "Some or all tau2 are non-positive or infinite. Ensure all tau2 > 0 and not infinite or use default value"
  }
}

check_Z <- function(Z, Y) {
  err_messages <- character()
  # dimensions don't match num_region num_group num_time
  if (!all(dim(Z) == dim(Y))) {
    err_messages <- c(
      err_messages,
      "Z is not an num_region x num_group x num_time array. Ensure dim(Z) == dim(Y) or use default value"
    )
  }
  # values are infinite
  if (!all(is.finite(Z))) {
    err_messages <- c(
      err_messages,
      "Z contains infinite values. Ensure all(is.finite(Z)) or use default value"
    )
  }
  err_messages
}

check_G <- function(G) {
  err_messages <- character()
  sig2 <- apply(G, 3, diag)
  gcor <- apply(G, 3, \(G) G[lower.tri(G)])
  # diagonals non-positive or infinite
  if (any((sig2 <= 0) | !is.finite(sig2))) {
    err_messages <- c(
      err_messages,
      "Diagonals of G contain non-positive values. Ensure all diag(G) > 0 and not infinite or use default value"
    )
  }
  # off-diagonal values are infinite
  if (!all(is.finite(gcor))) {
    err_messages <- c(
      err_messages,
      "Off-diagonals of G contain infinite values. Ensure all(is.finite(G)) or use default value"
    )
  }
  err_messages
}

check_rho <- function(rho) {
  # is non-positive or infinite
  if (any((rho <= 0) | !is.finite(rho))) {
    "rho contains non-positive values. Ensure all(rho > 0) and not infinite or use default value"
  }
}

check_Ag <- function(Ag) {
  err_messages <- character()
  # matrix is not symmetric
  if (!isSymmetric(Ag)) {
    err_messages <- c(
      err_messages,
      "Ag is not symmetric. Ensure Ag is symmetric or use default value"
    )
  }
  # values are infinite
  if (!all(is.finite(Ag))) {
    err_messages <- c(
      err_messages,
      "Ag contains infinite values. Ensure Ag is finite or use default value"
    )
  }
  # diagonals are not positive
  if (any(diag(Ag) <= 0)) {
    err_messages <- c(
      err_messages,
      "diag(Ag) contains non-positive values. Ensure diag(Ag) is positive or use default value"
    )
  }
  err_messages
}
