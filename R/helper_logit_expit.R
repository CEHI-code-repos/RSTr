logit <- function(x) {
  log(x / (1 - x))
}

expit <- function(x) {
  1 / (1 + exp(-x))
}

log_logit <- function(x, method) {
  if (method == "binomial") logit(x) else if (method == "poisson") log(x)
}

exp_expit <- function(x, method) {
  if (method == "binomial") expit(x) else if (method == "poisson") exp(x)
}
