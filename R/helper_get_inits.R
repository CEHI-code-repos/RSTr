get_inits <- function(RSTr_obj, inits, method) {
  UseMethod("get_inits")
}

#' @export
get_inits.car <- function(RSTr_obj, inits, method) {
  Y <- RSTr_obj$data$Y
  n <- RSTr_obj$data$n
  isl_id <- RSTr_obj$sp_data$isl_id
  n_group <- dim(Y)[[2]]
  n_time <- dim(Y)[[3]]
  n_island <- length(unique(isl_id))
  inits$beta <- inits$beta %||%
    get_inits_beta(Y, n, n_group, n_time, n_island, method)
  inits$lambda <- inits$lambda %||%
    get_inits_lambda(inits, Y, n, method, isl_id)

  inits$Z <- inits$Z %||%
    log_logit(inits$lambda, method) -
    inits$beta[isl_id, , , drop = FALSE]
  inits$tau2 <- inits$tau2 %||% matrix(1 / 100, n_group, n_time)
  inits$sig2 <- inits$sig2 %||% matrix(1 / 100, n_group, n_time)
  inits
}

#' @export
get_inits.mcar <- function(RSTr_obj, inits, method) {
  Y <- RSTr_obj$data$Y
  n <- RSTr_obj$data$n
  isl_id <- RSTr_obj$sp_data$isl_id
  n_group <- dim(Y)[[2]]
  n_time <- dim(Y)[[3]]
  n_island <- length(unique(isl_id))
  inits$beta <- inits$beta %||%
    get_inits_beta(Y, n, n_group, n_time, n_island, method)
  inits$lambda <- inits$lambda %||%
    get_inits_lambda(inits, Y, n, method, isl_id)
  inits$Z <- inits$Z %||%
    log_logit(inits$lambda, method) -
    inits$beta[isl_id, , , drop = FALSE]
  inits$tau2 <- inits$tau2 %||% matrix(1 / 100, n_group, n_time)
  inits$G <- inits$G %||%
    array(diag(n_group) / 7, dim = c(n_group, n_group, n_time))
  inits
}

#' @export
get_inits.mstcar <- function(RSTr_obj, inits, method) {
  Y <- RSTr_obj$data$Y
  n <- RSTr_obj$data$n
  isl_id <- RSTr_obj$sp_data$isl_id
  n_group <- dim(Y)[[2]]
  n_time <- dim(Y)[[3]]
  n_island <- length(unique(isl_id))
  inits$beta <- inits$beta %||%
    get_inits_beta(Y, n, n_group, n_time, n_island, method)
  inits$lambda <- inits$lambda %||%
    get_inits_lambda(inits, Y, n, method, isl_id)
  inits$Z <- inits$Z %||%
    log_logit(inits$lambda, method) -
    inits$beta[isl_id, , , drop = FALSE]
  inits$tau2 <- inits$tau2 %||% matrix(1 / 100, n_group, 1)
  inits$G <- inits$G %||%
    array(diag(n_group) / 7, dim = c(n_group, n_group, n_time))
  inits$rho <- inits$rho %||% matrix(0.95, 1, n_group)
  inits$Ag <- inits$Ag %||% diag(1 / 7, n_group)
  inits
}

get_inits_beta <- function(Y, n, n_group, n_time, n_island, method) {
  beta <- apply(Y, 2:3, sum, na.rm = TRUE) / apply(n, 2:3, sum)
  beta <- array(log_logit(beta, method), dim = c(n_group, n_time, n_island))
  beta[!is.finite(beta)] <- log_logit(sum(Y, na.rm = TRUE) / sum(n), method)
  beta <- aperm(beta, c(3, 1, 2))
  beta
}

get_inits_lambda <- function(inits, Y, n, method, isl_id) {
  limit_lo <- 0
  limit_hi <- ifelse(method == "binomial", 1, Inf)
  lambda <- Y / n
  lambda_na <- (lambda <= limit_lo) | (lambda >= limit_hi) | (is.na(lambda))
  if (any(lambda_na)) {
    lambda[lambda_na] <- exp_expit(inits$beta, method)[isl_id, , ][lambda_na]
  }
  lambda
}
