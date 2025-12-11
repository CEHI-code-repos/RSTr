get_params <- function(data, seed, method, model, name, dir, perc_ci, restricted, A, m0, update_rho, impute_lb, impute_ub) {
  params <- list(
    batch = 0,
    total = 0,
    method = method,
    model = model,
    dimnames = dimnames(data$Y),
    name = name,
    dir = dir,
    perc_ci = perc_ci,
    missing_Y = FALSE,
    age_standardized = FALSE,
    suppressed = FALSE
  )
  if (!all(is.finite(data$Y))) {
    params$missing_Y <- TRUE
    params$impute_lb <- impute_lb
    params$impute_ub <- impute_ub
    params$miss <- which(!is.finite(data$Y))
  }
  if (!is.null(seed)) {
    set.seed(seed)
    params$seed <- .Random.seed
  }
  if (model %in% c("ucar", "eucar")) {
    params$restricted <- restricted
    if (restricted) {
      if (is.null(A)) A <- array(6 / dim(data$Y)[2], dim = dim(data$Y)[-1])
      if (is.null(m0)) m0 <- 3
      params$A <- array(A, dim = dim(data$Y)[2:3])
      params$m0 <- m0
    }
  } else if (model == "mstcar") {
    params$update_rho <- update_rho
  }
  params
}