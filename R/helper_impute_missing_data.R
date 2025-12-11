impute_missing_data <- function(RSTr_obj) {
  lambda <- RSTr_obj$current_sample$lambda
  params <- RSTr_obj$params
  method <- params$method
  miss <- params$miss
  impute_lb <- params$impute_lb
  impute_ub <- params$impute_ub
  data <- RSTr_obj$data
  Y <- data$Y
  n <- data$n
  if (method == "binomial") {
    rate <- expit(lambda[miss])
    rp <- stats::runif(
      length(miss),
      stats::pbinom(impute_lb - 0.1, round(n[miss]), rate),
      stats::pbinom(impute_ub + 0.1, round(n[miss]), rate)
    )
    Y[miss] <- stats::qbinom(rp, round(n[miss]), rate)
  }
  if (method == "poisson") {
    rate <- n[miss] * exp(lambda[miss])
    rp <- stats::runif(
      length(miss),
      stats::ppois(impute_lb - 0.1, rate),
      stats::ppois(impute_ub + 0.1, rate)
    )
    Y[miss] <- stats::qpois(rp, rate)
  }
  data$Y <- Y
  RSTr_obj$data <- data
  RSTr_obj
}