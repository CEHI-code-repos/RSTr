#' Create CAR model
#'
#' \code{*car()} generates an \code{RSTr} model object, samples, and estimates for either an MSTCAR, MCAR, RCAR, or CAR model.
#'
#' @param name Name of model and corresponding folder.
#' @param data Dataset including mortality (Y) and population (n) information.
#' @param adjacency Dataset including adjacency information.
#' @param dir Directory where model will live.
#' @param seed Set of random seeds to use for data replication.
#' @param perc_ci The percentage of the desired estimate credible interval. Defaults to 95 percent (0.95).
#' @param iterations The number of iterations to run the model for.
#' @param show_plots If set to \code{FALSE}, suppresses traceplots.
#' @param verbose If set to \code{FALSE}, suppresses model progress messages.
#' @param ignore_checks If set to \code{TRUE}, skips model validation.
#' @param method Run model with either Binomial data or Poisson data.
#' @param impute_bounds If counts are suppressed for privacy reasons, \code{impute_bounds} is the lower/upper bound of suppression, typically 0 or 1 and 10, respectively.
#' @param inits Optional list of initial conditions for each parameter.
#' @param priors Optional list of priors for updates.
#' @param m0 For RCAR models, baseline neighbor count by region.
#' @param A For RCAR models, describes maximum intensity of smoothing between regions.
#' @param update_rho For MSTCAR models, controls whether rho update is performed.
#' @returns An \code{RSTr} model object.
#' @examples
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' # MSTCAR model
#' mod_mst <- mstcar(
#'   name = "test",
#'   data = data_min,
#'   adjacency = adj_min,
#'   dir = tempdir(),
#'   show_plots = FALSE,
#'   verbose = FALSE
#' )
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
car <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_bounds = NULL,
  inits = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "sig2", "tau2")
  RSTr_obj <- initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    seed = seed,
    perc_ci = perc_ci,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_bounds = impute_bounds,
    inits = inits,
    priors = priors,
    model = "car",
    pars = pars,
    restricted = FALSE
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  check_informativeness(RSTr_obj)
  RSTr_obj
}

#' Initialize Restricted CAR model
#' @rdname car
#' @export
rcar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  A = NULL,
  m0 = NULL,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_bounds = NULL,
  inits = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "sig2", "tau2")
  RSTr_obj <- initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    seed = seed,
    perc_ci = perc_ci,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_bounds = impute_bounds,
    inits = inits,
    priors = priors,
    model = "rcar",
    pars = pars,
    restricted = TRUE,
    A = A,
    m0 = m0
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}

#' Initialize MCAR model
#' @rdname car
#' @export
mcar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_bounds = NULL,
  inits = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "G", "tau2")
  RSTr_obj <- initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    seed = seed,
    perc_ci = perc_ci,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_bounds = impute_bounds,
    inits = inits,
    priors = priors,
    model = "mcar",
    pars = pars
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}

#' Initialize MSTCAR model
#' @rdname car
#' @export
mstcar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_bounds = NULL,
  inits = NULL,
  priors = NULL,
  update_rho = FALSE
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "G", "Ag", "tau2")
  if (update_rho) {
    pars <- c(pars, "rho")
  }
  RSTr_obj <- initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    seed = seed,
    perc_ci = perc_ci,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_bounds = impute_bounds,
    inits = inits,
    priors = priors,
    model = "mstcar",
    pars = pars,
    update_rho = update_rho
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}

initialize_model <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = "binomial",
  impute_bounds = NULL,
  inits = NULL,
  priors = NULL,
  model = c("mstcar", "car", "mcar"),
  pars,
  restricted = NULL,
  A = NULL,
  m0 = NULL,
  update_rho = NULL
) {
  RSTr_obj <- new_model(model, data, restricted, update_rho)
  RSTr_obj$params <- get_params(
    RSTr_obj$data,
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
  )
  RSTr_obj$sp_data <- get_sp_data(adjacency)
  RSTr_obj <- get_priors(RSTr_obj, priors)
  RSTr_obj$inits <- get_inits(
    RSTr_obj,
    inits,
    method
  )
  RSTr_obj$sample <- RSTr_obj$inits
  if (!ignore_checks) {
    validate_model(RSTr_obj)
  }
  create_model_directory(name, dir, pars)
  save_model(RSTr_obj)
  RSTr_obj
}
