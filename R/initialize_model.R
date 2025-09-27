#' Initialize model
#' @noRd
initialize_model <- function(
  name,
  data,
  adjacency,
  dir,
  show_plots,
  ignore_checks,
  method,
  impute_lb,
  impute_ub,
  seed,
  inits,
  priors,
  model = c("mstcar", "ucar", "mcar"),
  restricted = NULL,
  m0 = NULL,
  A = NULL,
  rho_up = NULL
) {
  model <- match.arg(model)
  if (model == "ucar") {
    if (is.null(dim(data$Y))) {
      data <- lapply(data, \(x) array(x, dim = c(length(x), 1, 1), dimnames = list(names(x))))
    } else if (length(dim(data$Y)) == 2) {
      data <- lapply(data, \(x) array(x, dim = c(dim(x), 1), dimnames = dimnames(x)))
    }
  }
  miss <- which(!is.finite(data$Y))
  if (!ignore_checks) {
    check_data(data)
  }
  if (model == "mstcar" & show_plots) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    par(mfrow = c(1, 2))
    plot(dimnames(data$Y)[[3]], apply(data$Y, 3, sum, na.rm = TRUE), xlab = "Year", ylab = "Events")
    plot(dimnames(data$Y)[[3]], apply(data$n, 3, sum), xlab = "Year", ylab = "Population")
  }
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  if (!dir.exists(paste0(dir, name))) {
    dir.create(paste0(dir, name))
  }
  pars <- list(
    "ucar" = c("theta", "beta", "Z", "sig2", "tau2"),
    "mcar" = c("theta", "beta", "Z", "G", "tau2"),
    "mstcar" = c("theta", "beta", "Z", "G", "Ag", "tau2")
  )[[model]]
  if (model == "mstcar") if (rho_up) pars <- c(pars, "rho")
  for (par in pars) {
    if (!dir.exists(paste0(dir, name, "/", par))) {
      dir.create(paste0(dir, name, "/", par))
    }
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  params <- list(
    seed = .Random.seed,
    batch = 0,
    total = 0,
    model = model,
    method = method,
    dimnames = dimnames(data$Y)
  )
  if (model == "ucar") {
    params$restricted <- restricted
    if (restricted) {
      params$A <- array(1 / dim(data$Y)[2], dim = dim(data$Y)[-1])
      params$m0 <- m0
    }
  } else if (model == "mstcar") {
    params$rho_up <- rho_up
  }
  if (length(miss)) {
    params$impute_lb <- impute_lb
    params$impute_ub <- impute_ub
  }

  spatial_data <- get_spatial_data(adjacency, ignore_checks)
  if (model == "ucar") {
    inits <- get_inits_u(inits, data, spatial_data, method, ignore_checks)
  } else if (model == "mcar") {
    inits <- get_inits_m(inits, data, spatial_data$island_id, method, ignore_checks)
  } else if (model == "mstcar") {
    inits <- get_inits_mst(inits, data, spatial_data$island_id, method, ignore_checks)
  }
  priors <- get_priors(priors, data, model, ignore_checks)

  saveRDS(data, file = paste0(dir, name, "/data.Rds"))
  saveRDS(params, file = paste0(dir, name, "/params.Rds"))
  saveRDS(spatial_data, file = paste0(dir, name, "/spatial_data.Rds"))
  saveRDS(priors, file = paste0(dir, name, "/priors.Rds"))
  saveRDS(inits, file = paste0(dir, name, "/inits.Rds"))
  message("Model ready!")
}

#' Initialize CAR model
#'
#' This function performs checks and prepares data for use with either an MSTCAR, MCAR, or UCAR model. This function additionally specifies all of the model parameters, such as model type, event data type, intensity of smoothing in the UCAR model, and more.
#'
#' @param name Name of model and corresponding folder
#' @param data Dataset including mortality (Y) and population (n) information
#' @param adjacency Dataset including adjacency information
#' @param dir Directory where model will live
#' @param show_plots If set to \code{FALSE}, suppresses check plots generated for MSTCAR models
#' @param ignore_checks If set to \code{TRUE}, ignores data checks. Only use if you are certain that your input data is
#' correct and you are encountering bugs during setup
#' @param method Run model with either Binomial data or Poisson data
#' @param impute_lb If counts are suppressed for privacy reasons, \code{impute_lb} is lower bound of suppression, typically 0 or 1
#' @param impute_ub If counts are suppressed for privacy reasons, \code{impute_ub} is upper bound of suppression, typically 10
#' @param seed Set of random seeds to use for data replication
#' @param inits Optional list of initial conditions for each parameter
#' @param priors Optional list of priors for updates
#' @param m0 For restricted UCAR models, baseline neighbor count by region
#' @param A For restricted UCAR models, describes intensity of smoothing between regions
#' @param rho_up For MSTCAR models, controls whether rho update is performed for MSTCAR models
#' @returns No output, only sets up model and saves files to directory
#' @examples
#' # Initialize an MSTCAR model
#' initialize_mstcar(name = "test", data = miheart, adjacency = miadj, dir = tempdir())
#' # Initialize an MCAR model
#' data_m <- lapply(miheart, \(x) x[, , "1979"])
#' initialize_mcar("test", data_m, miadj, tempdir())
#' # Initialize an MCAR model with Poisson-distributed event data
#' initialize_mcar("test", data_m, miadj, tempdir(), method = "poisson")
#' # Initialize a restricted UCAR model
#' data_u <- lapply(miheart, \(x) x[, "65-74", "1979"])
#' initialize_ucar_restricted("test", data_u, miadj, tempdir(), A = 6)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
initialize_ucar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = 1,
  impute_ub = 9,
  seed = NULL,
  inits = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    seed = seed,
    inits = inits,
    priors = priors,
    model = "ucar",
    restricted = FALSE,
  )
}

#' Initialize Restricted UCAR model
#' @rdname initialize_ucar
#' @export
initialize_ucar_restricted <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  A = NULL,
  m0 = NULL,
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = 1,
  impute_ub = 9,
  seed = NULL,
  inits = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  if (is.null(A)) A <- 6
  if (is.null(m0)) m0 <- 3
  initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    seed = seed,
    inits = inits,
    priors = priors,
    model = "ucar",
    restricted = TRUE,
    A = A,
    m0 = m0
  )
}

#' Initialize MCAR model
#' @rdname initialize_ucar
#' @export
initialize_mcar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = 1,
  impute_ub = 9,
  seed = NULL,
  inits = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    seed = seed,
    inits = inits,
    priors = priors,
    model = "mcar"
  )
}

#' Initialize MSTCAR model
#' @rdname initialize_ucar
#' @export
initialize_mstcar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = 1,
  impute_ub = 9,
  seed = NULL,
  inits = NULL,
  priors = NULL,
  rho_up = FALSE
) {
  method <- match.arg(method)
  initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    seed = seed,
    inits = inits,
    priors = priors,
    model = "mstcar",
    rho_up = rho_up
  )
}