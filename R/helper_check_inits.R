check_inits <- function(RSTr_obj) {
  UseMethod("check_inits")
}

#' @export
check_inits.car <- function(RSTr_obj) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  n_group <- dim(Y)[[2]]
  n_time <- dim(Y)[[3]]
  n_island <- RSTr_obj$sp_data$n_island
  chk <- c("beta", "tau2", "lambda", "Z", "sig2")
  check_missing_inits(RSTr_obj, chk)
  # Check for warnings
  check_unused_inits(RSTr_obj, chk)
  # Check for errors
  c(
    check_beta(RSTr_obj$inits$beta, n_island, n_group, n_time),
    check_lambda(RSTr_obj$inits$lambda, Y, method),
    check_sig2(RSTr_obj$inits$sig2),
    check_tau2(RSTr_obj$inits$tau2),
    check_Z(RSTr_obj$inits$Z, Y)
  )
}

#' @export
check_inits.mcar <- function(RSTr_obj) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  n_group <- dim(Y)[[2]]
  n_time <- dim(Y)[[3]]
  n_island <- RSTr_obj$sp_data$n_island
  chk <- c("beta", "tau2", "lambda", "Z", "G")
  check_missing_inits(RSTr_obj, chk)
  # Check for warnings
  check_unused_inits(RSTr_obj, chk)

  # Check for errors
  c(
    check_beta(RSTr_obj$inits$beta, n_island, n_group, n_time),
    check_lambda(RSTr_obj$inits$lambda, Y, method),
    check_G(RSTr_obj$inits$G),
    check_tau2(RSTr_obj$inits$tau2),
    check_Z(RSTr_obj$inits$Z, Y)
  )
}

#' @export
check_inits.mstcar <- function(RSTr_obj) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  n_group <- dim(Y)[[2]]
  n_time <- dim(Y)[[3]]
  n_island <- RSTr_obj$sp_data$n_island
  chk <- c("beta", "tau2", "lambda", "Z", "G", "rho", "Ag")
  check_missing_inits(RSTr_obj, chk)
  # Check for warnings
  check_unused_inits(RSTr_obj, chk)
  # Check for errors
  c(
    check_beta(RSTr_obj$inits$beta, n_island, n_group, n_time),
    check_lambda(RSTr_obj$inits$lambda, Y, method),
    check_G(RSTr_obj$inits$G),
    check_rho(RSTr_obj$inits$rho),
    check_tau2(RSTr_obj$inits$tau2),
    check_Z(RSTr_obj$inits$Z, Y),
    check_Ag(RSTr_obj$inits$Ag)
  )
}

check_missing_inits <- function(RSTr_obj, chk) {
  miss <- sapply(seq_along(chk), \(x) {
    !any(names(RSTr_obj$inits) == chk[x])
  })
  if (sum(miss)) {
    stop(
      "One or more objects missing from list 'inits': ",
      toString(chk[miss])
    )
  }
}

check_unused_inits <- function(RSTr_obj, chk) {
  chk_elem <- which(!(names(RSTr_obj$inits) %in% chk))
  if (length(chk_elem)) {
    warning(paste(
      "Unused elements of list 'inits':",
      toString(names(RSTr_obj$inits)[chk_elem])
    ))
  }
}

check_beta <- function(beta, n_island, n_group, n_time) {
  err_messages <- character()
  # dimensions don't match n_island n_group n_time
  if (!all(dim(beta) == c(n_island, n_group, n_time))) {
    err_messages <- c(
      err_messages,
      "beta is not an n_island x n_group x n_time array. Ensure dim(beta) == n_island x n_group x n_time or use default value"
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
  # dimensions don't match n_region n_group n_time
  if (!all(dim(lambda) == dim(Y))) {
    err_messages <- c(
      err_messages,
      "lambda is not a n_region x n_group x n_time array. Ensure dim(lambda) == dim(Y) or use default value"
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
  # dimensions don't match n_region n_group n_time
  if (!all(dim(Z) == dim(Y))) {
    err_messages <- c(
      err_messages,
      "Z is not an n_region x n_group x n_time array. Ensure dim(Z) == dim(Y) or use default value"
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
