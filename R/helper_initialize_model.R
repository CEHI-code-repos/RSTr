validate_model <- function(RSTr_obj) {
  errout <- c(
    check_data(RSTr_obj),
    check_inits(RSTr_obj),
    check_priors(RSTr_obj)
  )
  display_errors(errout)
}

display_errors <- function(errout) {
  if (length(errout)) {
    stop(paste0(
      length(errout),
      " error(s) found:\n",
      paste(errout, collapse = "\n ")
    ))
  }
}

prepare_data <- function(data) {
  if (is.null(dim(data$Y))) {
    data <- lapply(data, \(x) {
      array(x, dim = c(length(x), 1, 1), dimnames = list(names(x)))
    })
  } else if (length(dim(data$Y)) == 2) {
    data <- lapply(data, \(x) {
      array(x, dim = c(dim(x), 1), dimnames = dimnames(x))
    })
  }
  data
}

post_sampler_output <- function(RSTr_obj) {
  samples <- load_samples(RSTr_obj)
  medians <- get_medians(samples)
  ci <- get_credible_interval(samples, RSTr_obj$params$perc_ci)
  rel_prec <- get_relative_precision(medians, ci)
  RSTr_obj$medians <- medians
  RSTr_obj$ci <- ci
  RSTr_obj$rel_prec <- rel_prec
  RSTr_obj
}
