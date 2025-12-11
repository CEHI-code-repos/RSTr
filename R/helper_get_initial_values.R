get_initial_values <- function(RSTr_obj, initial_values, method) {
  UseMethod("get_initial_values")
}

#' @export
get_initial_values.ucar <- function(RSTr_obj, initial_values, method) {
  Y <- RSTr_obj$data$Y
  n <- RSTr_obj$data$n
  island_id <- RSTr_obj$spatial_data$island_id
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- length(unique(island_id))
  # beta
  if (is.null(initial_values$beta)) {
    initial_values$beta <- get_initial_values_beta(Y, n, num_group, num_time, num_island, method)
  }
  # lambda
  if (is.null(initial_values$lambda)) {
    initial_values$lambda <- get_initial_values_lambda(initial_values, Y, n, method, island_id)
  }
  # Z
  if (is.null(initial_values$Z)) {
    initial_values$Z <- log_logit(initial_values$lambda, method) - initial_values$beta[island_id + 1, , , drop = FALSE]
  }
  # tau2
  if (is.null(initial_values$tau2)) {
    initial_values$tau2 <- matrix(1 / 100, num_group, num_time)
  }
  if (is.null(initial_values$sig2)) {
    initial_values$sig2 <- matrix(1 / 100, num_group, num_time)
  }
  initial_values
}

#' @export
get_initial_values.mcar <- function(RSTr_obj, initial_values, method) {
  Y <- RSTr_obj$data$Y
  n <- RSTr_obj$data$n
  island_id <- RSTr_obj$spatial_data$island_id
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- length(unique(island_id))
  # beta
  if (is.null(initial_values$beta)) {
    initial_values$beta <- get_initial_values_beta(Y, n, num_group, num_time, num_island, method)
  }
  # lambda
  if (is.null(initial_values$lambda)) {
    initial_values$lambda <- get_initial_values_lambda(initial_values, Y, n, method, island_id)
  }
  # Z
  if (is.null(initial_values$Z)) {
    initial_values$Z <- log_logit(initial_values$lambda, method) - initial_values$beta[island_id + 1, , , drop = FALSE]
  }
  # tau2
  if (is.null(initial_values$tau2)) {
    initial_values$tau2 <- matrix(1 / 100, num_group, num_time)
  }
  if (is.null(initial_values$G)) {
    initial_values$G <- array(diag(num_group) / 7, dim = c(num_group, num_group, num_time))
  }
  initial_values
}

#' @export
get_initial_values.mstcar <- function(RSTr_obj, initial_values, method) {
  Y <- RSTr_obj$data$Y
  n <- RSTr_obj$data$n
  island_id <- RSTr_obj$spatial_data$island_id
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- length(unique(island_id))
  # beta
  if (is.null(initial_values$beta)) {
    initial_values$beta <- get_initial_values_beta(Y, n, num_group, num_time, num_island, method)
  }
  # lambda
  if (is.null(initial_values$lambda)) {
    initial_values$lambda <- get_initial_values_lambda(initial_values, Y, n, method, island_id)
  }
  # Z
  if (is.null(initial_values$Z)) {
    initial_values$Z <- log_logit(initial_values$lambda, method) - initial_values$beta[island_id + 1, , , drop = FALSE]
  }
  # tau2
  if (is.null(initial_values$tau2)) {
    initial_values$tau2 <- matrix(1 / 100, num_group, 1)
  }
  if (is.null(initial_values$G)) {
    initial_values$G <- array(diag(num_group) / 7, dim = c(num_group, num_group, num_time))
  }
  # rho
  if (is.null(initial_values$rho)) {
    initial_values$rho <- matrix(0.95, 1, num_group)
  }
  # Ag
  if (is.null(initial_values$Ag)) {
    initial_values$Ag <- diag(1 / 7, num_group)
  }
  initial_values
}

get_initial_values_beta <- function(Y, n, num_group, num_time, num_island, method) {
  beta <- apply(Y, 2:3, sum, na.rm = TRUE) / apply(n, 2:3, sum)
  beta <- array(log_logit(beta, method), dim = c(num_group, num_time, num_island))
  beta[!is.finite(beta)] <- log_logit(sum(Y, na.rm = TRUE) / sum(n), method)
  beta <- aperm(beta, c(3, 1, 2))
  beta
}

get_initial_values_lambda <- function(initial_values, Y, n, method, island_id) {
  lower_limit <- 0
  upper_limit <- ifelse(method == "binomial", 1, Inf)
  lambda <- Y / n
  lambda_unsupported <- (lambda <= lower_limit) | (lambda >= upper_limit)
  if (any(lambda_unsupported)) {
    lambda[lambda_unsupported] <- exp_expit(initial_values$beta, method)[island_id + 1, , ][lambda_unsupported]
  }
  lambda
}