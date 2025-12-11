#' Load MCMC samples
#'
#' \code{load_samples()} gathers samples saved for model \code{RSTr_obj}. By default, loads the rate estimate samples \code{lambda}, but any model parameters can be loaded. Users can also specify a burn-in period.
#' 
#' @param RSTr_obj \code{RSTr} model object to load in samples from.
#' @param param Which parameter samples to load.
#' @param burn Number of burn-in samples to discard.
#' @returns An \code{array} of samples from model \code{RSTr_obj}.
#' @examples
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
#' samples <- load_samples(mod_mst) * 1e5
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
load_samples <- function(RSTr_obj, param = "lambda", burn = 2000) {
  params <- RSTr_obj$params
  name <- RSTr_obj$params$name
  dir <- RSTr_obj$params$dir
  mar <- c(lambda = 4, beta = 4, Z = 4, G = 4, Ag = 3, tau2 = 3, sig2 = 3, rho = 2)
  if (params$model == "mstcar") mar["tau2"] <- 2
  batch <- which(seq_len(params$batch) * 100 > burn)
  files <- paste0(dir, "/", name, "/", param, "/", param, "_out_", batch, ".Rds")
  output <- abind::abind(lapply(files, readRDS), along = mar[param])
  dims <- params$dimnames
  its <- seq(burn + 10, max(batch) * 100, by = 10)
  if (param == "beta") {
    num_island <- readRDS(paste0(dir, name, "/spatial_data.Rds"))$num_island
    dimnames(output) <- list(island = seq_len(num_island), group = dims[[2]], time = dims[[3]], its = its)
  } else if (param %in% c("Z", "lambda")) {
    dimnames(output) <- c(dims, list(its = its))
  } else if (param == "rho") {
    dimnames(output) <- list(group = dims[[2]], its = its)
  } else if (param == "tau2") {
    if (params$model == "mstcar") {
      dimnames(output) <- list(group = dims[[2]], its = its)
    } else if (params$model %in% c("ucar", "mcar")) {
      dimnames(output) <- list(group = dims[[2]], time = dims[[3]], its = its)
    }
  } else if (param == "Ag") {
    dimnames(output) <- list(group1 = dims[[2]], group2 = dims[[2]], its = its)
  } else if (param == "G") {
    dimnames(output) <- list(group1 = dims[[2]], group2 = dims[[2]], time = dims[[3]], its = its)
  } else if (param == "sig2") {
    dimnames(output) <- list(group = dims[[2]], time = dims[[3]], its = its)
  }
  output
}