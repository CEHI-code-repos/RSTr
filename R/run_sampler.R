#' Run Gibbs sampler
#' @param name Name of model and corresponding folder
#' @param dir Directory where model lives
#' @param iterations Specifies number of iterations to run
#' @param .show_plots Show or hide traceplots as 
#' @param .discard_burnin If set to \code{TRUE}, won't save burn-in samples
#'
#' @export
run_sampler <- function(name, dir = getwd(), iterations = 6000, .show_plots = TRUE, .discard_burnin = FALSE) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  iterations = iterations - iterations %% 100
  model <- readRDS(paste0(dir, name, "/params.Rds"))$model
  if (model == "ucar") {
    gibbs_u(name, dir, iterations, .show_plots, .discard_burnin)
  }
  if (model == "mcar") {
    gibbs_m(name, dir, iterations, .show_plots, .discard_burnin)
  }
  if (model == "mstcar") {
    gibbs_mst(name, dir, iterations, .show_plots, .discard_burnin)
  }
}
