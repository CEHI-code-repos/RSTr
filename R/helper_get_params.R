get_params <- function(
  data,
  seed,
  method,
  model,
  name,
  dir,
  perc_ci,
  restricted,
  A,
  m0,
  update_rho,
  impute_bounds
) {
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
    params$impute_bounds <- impute_bounds %||% c(0, 10)
    params$miss <- which(!is.finite(data$Y))
  }
  if (!is.null(seed)) {
    set.seed(seed)
    params$seed <- .Random.seed
  }
  if (model %in% c("car", "rcar")) {
    params$restricted <- restricted
    if (restricted) {
      A <- A %||% array(6, dim = dim(data$Y)[-1])
      A <- array(A, dim = dim(data$Y)[2:3])
      m0 <- m0 %||% 3
      params$A <- A
      params$m0 <- m0
    }
  } else if (model == "mstcar") {
    params$update_rho <- update_rho
  }
  params
}
