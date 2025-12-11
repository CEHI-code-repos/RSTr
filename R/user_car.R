#' Create CAR model
#' 
#' \code{*car()} generates an \code{RSTr} model object, samples, and estimates for either an MSTCAR, MCAR, EUCAR, or UCAR model.
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
#' @param impute_lb If counts are suppressed for privacy reasons, \code{impute_lb} is lower bound of suppression, typically 0 or 1.
#' @param impute_ub If counts are suppressed for privacy reasons, \code{impute_ub} is upper bound of suppression, typically 10.
#' @param initial_values Optional list of initial conditions for each parameter.
#' @param priors Optional list of priors for updates.
#' @param m0 For EUCAR models, baseline neighbor count by region.
#' @param A For EUCAR models, describes maximum intensity of smoothing between regions.
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
ucar <- function(
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
  impute_lb = NULL,
  impute_ub = NULL,
  initial_values = NULL,
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
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    initial_values = initial_values,
    priors = priors,
    model = "ucar",
    pars = pars,
    restricted = FALSE
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}

#' Initialize Enhanced UCAR model
#' @rdname ucar
#' @export
eucar <- function(
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
  impute_lb = NULL,
  impute_ub = NULL,
  initial_values = NULL,
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
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    initial_values = initial_values,
    priors = priors,
    model = "eucar",
    pars = pars,
    restricted = TRUE,
    A = A,
    m0 = m0
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}

#' Initialize MCAR model
#' @rdname ucar
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
  impute_lb = NULL,
  impute_ub = NULL,
  initial_values = NULL,
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
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    initial_values = initial_values,
    priors = priors,
    model = "mcar",
    pars = pars
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}

#' Initialize MSTCAR model
#' @rdname ucar
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
  impute_lb = NULL,
  impute_ub = NULL,
  initial_values = NULL,
  priors = NULL,
  update_rho = FALSE
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "G", "Ag", "tau2")
  if (update_rho) pars <- c(pars, "rho")
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
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    initial_values = initial_values,
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
  impute_lb = 0,
  impute_ub = 10,
  initial_values = NULL,
  priors = NULL,
  model = c("mstcar", "ucar", "mcar"),
  pars,
  restricted = NULL,
  A = NULL,
  m0 = NULL,
  update_rho = NULL
) {
  RSTr_obj <- create_new_model(model, data, restricted, update_rho)
  RSTr_obj$params <- get_params(RSTr_obj$data, seed, method, model, name, dir, perc_ci, restricted, A, m0, update_rho, impute_lb, impute_ub)
  RSTr_obj$spatial_data <- get_spatial_data(adjacency)
  RSTr_obj <- get_priors(RSTr_obj, priors)
  RSTr_obj$initial_values <- get_initial_values(RSTr_obj, initial_values, method)
  RSTr_obj$current_sample <- RSTr_obj$initial_values
  if (!ignore_checks) validate_model(RSTr_obj)
  create_model_directory(name, dir, pars)
  save_model(RSTr_obj)
  RSTr_obj
}