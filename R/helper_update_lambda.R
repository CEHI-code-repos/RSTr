clamp_array <- function(arr, min, max) {
  pmax(pmin(arr, max), min)
}

rnorm_array <- function(mean, sd) {
  mean + stats::rnorm(length(mean), sd = sd)
}

update_lambda <- function(RSTr_obj) {
  l_0 <- RSTr_obj$sample$lambda
  l_sd <- RSTr_obj$priors$lambda_sd
  l_star <- rnorm_array(l_0, l_sd) |> clamp_array(-100.0, 15.0)

  method <- RSTr_obj$params$method
  tau2 <- array(RSTr_obj$sample$tau2, dim = dim(l_0)[c(2, 3, 1)]) |>
    aperm(c(3, 1, 2))
  beta <- RSTr_obj$sample$beta[RSTr_obj$sp_data$isl_id + 1, , , drop = FALSE]
  Z <- RSTr_obj$sample$Z
  n <- RSTr_obj$data$n
  Y <- RSTr_obj$data$Y
  delta <- l_star - l_0

  r1 <- Y * delta
  r2 <- switch(
    method,
    "binomial" = n * (log(exp(l_star) + 1) - log(exp(l_0) + 1)),
    "poisson" = n * (exp(l_star) - exp(l_0))
  )
  r3 <- delta * (l_star + l_0 - 2.0 * (beta + Z))
  r <- exp(r1 - r2 - 0.5 * r3 / tau2)

  update <- r >= stats::runif(length(l_0))
  RSTr_obj$sample$lambda[update] <- l_star[update]
  RSTr_obj$priors$lambda_acpt[update] <- RSTr_obj$priors$lambda_acpt[update] + 1
  RSTr_obj
}
