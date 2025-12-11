#' Run Gibbs sampler
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#' @noRd
run_sampler <- function(RSTr_obj, iterations = 6000, show_plots = TRUE, verbose = TRUE) {
  iterations <- iterations - iterations %% 100
  sampler_start <- Sys.time()
  missing_Y <- RSTr_obj$params$missing_Y
  start_batch <- RSTr_obj$params$batch
  total <- RSTr_obj$params$total
  method <- RSTr_obj$params$method
  batches <- seq(start_batch + 1, start_batch + iterations / 100)
  if (verbose) message("NAs detected in Y. Events will be imputed for missing values")
  if (verbose) message("Starting sampler on Batch ", start_batch + 1, " at ", format(Sys.time(), "%a %b %d %X"))
  for (batch in batches) {
    if (verbose) display_progress(batch, max(batches), total, 0, sampler_start)
    output <- stats::setNames(vector("list", length(RSTr_obj$current_sample)), names(RSTr_obj$current_sample))
    RSTr_obj$current_sample$lambda <- log_logit(RSTr_obj$current_sample$lambda, method)
    RSTr_obj <- convert_index(RSTr_obj, "zero")
    for (it in 1:100) {
      if (missing_Y) RSTr_obj <- impute_missing_data(RSTr_obj)
      RSTr_obj <- update_current_sample(RSTr_obj)
      if (it %% 10 == 0) output <- append_to_output(output, RSTr_obj)
      if (verbose) display_progress(batch, max(batches), total, it, sampler_start)
    }
    RSTr_obj <- convert_index(RSTr_obj, "one")
    output <- prepare_output(output, method)
    RSTr_obj$current_sample$lambda <- exp_expit(RSTr_obj$current_sample$lambda, method)
    RSTr_obj <- update_priors_sd(RSTr_obj)
    RSTr_obj <- update_params(RSTr_obj, batch)
    save_model(RSTr_obj)
    save_output(output, batch, RSTr_obj$params$dir, RSTr_obj$params$name)

    if (show_plots) {
      if (!exists("plots")) plots <- NULL
      plots <- update_plots(plots, output, RSTr_obj$params$batch, start_batch)
      plot(plots, xlab = "Iterations", main = "Traceplots")
    } 
  }
  RSTr_obj
}

#' Append new values to plots
#' @noRd
update_plots <- function(plots, output, batch, start_batch) {
  if (start_batch < 40) {
    start <- min(batch * 100 / 2, 2000) + 10
  } else {
    start <- start_batch * 100 + 10
  }
  plots <- rbind(plots, sapply(output, extract_last_margin))
  if (start < 2000) plots <- plots[-(1:5), ]
  stats::ts(plots, start = start, frequency = 0.1)
}

#' Append new values to output
#' @noRd
append_to_output <- function(output, RSTr_obj) {
  current_sample <- RSTr_obj$current_sample
  along <- sapply(current_sample, \(par) length(dim(par)) + 1) 
  mapply(abind::abind, output, current_sample, along = along)
}

#' Update params list inside of RSTr_obj
#' @noRd
update_params <- function(RSTr_obj, current_batch) {
  params <- RSTr_obj$params
  params$total <- params$total + 100
  params$batch <- current_batch
  RSTr_obj$params <- params
  RSTr_obj
}

#' Prepare output to be saved and used for plots
#' @noRd
prepare_output <- function(output, method) {
  output$lambda <- exp_expit(output$lambda, method)
  # remove parameter from `output` if no changes detected
  difftest <- lapply(output, \(par) diff(extract_last_margin(par)))
  output <- output[!sapply(difftest, \(par) all(par == 0))]
  output
}