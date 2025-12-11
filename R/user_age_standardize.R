#' Age-standardize model objects
#' 
#' Age-standardizes samples using a standard population for an \code{RSTr} model object.
#' 
#' @param RSTr_obj An \code{RSTr} model object.
#' @param std_pop A vector of standard populations.
#' @param new_name The name to assign to the age-standardized group.
#' @param groups A vector of either indices for each group or a vector of strings for each group name. If set to \code{NULL}, will use all groups in the dataset.
#' @returns An \code{RSTr} object with age-standardized estimates.
#' @examples
#' std_pop <- c(113154, 100640, 95799)
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
#' # age-standardize by all age groups
#' mod_mst <- age_standardize(mod_mst, std_pop, "35-64")
#' # Add onto age-standardized estimates. Age-standardize only by the first two age groups
#' mod_mst <- age_standardize(mod_mst, std_pop[1:2], "35-54", groups = 1:2)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
age_standardize <- function(RSTr_obj, std_pop, new_name, groups = NULL) {
  RSTr_obj$params$age_standardized <- TRUE
  samples <- load_samples(RSTr_obj)
  if (is.null(groups)) groups <- seq_len(dim(samples)[2])
  data <- lapply(RSTr_obj$data, aggregate_count, 2, groups, TRUE, new_name)
  data <- lapply(data, \(x) x[, new_name, , drop = FALSE])
  samples <- subset_array(samples, 2, groups)
  samples <- standardize_samples(samples, std_pop, 2, groups, TRUE, new_name)[, new_name, , , drop = FALSE]
  medians <- get_medians(samples)
  credible_interval <- get_credible_interval(samples)
  relative_precision <- get_relative_precision(medians, credible_interval)
  RSTr_obj$medians_as <- erase_duplicates(abind::abind(RSTr_obj$medians_as, medians, along = 2))
  RSTr_obj$data_as$Y <- erase_duplicates(abind::abind(RSTr_obj$data_as$Y, data$Y, along = 2))
  RSTr_obj$data_as$n <- erase_duplicates(abind::abind(RSTr_obj$data_as$n, data$n, along = 2))
  RSTr_obj$credible_interval_as$lower <- erase_duplicates(abind::abind(RSTr_obj$credible_interval_as$lower, credible_interval$lower, along = 2))
  RSTr_obj$credible_interval_as$upper <- erase_duplicates(abind::abind(RSTr_obj$credible_interval_as$upper, credible_interval$upper, along = 2))
  RSTr_obj$relative_precision_as <- erase_duplicates(abind::abind(RSTr_obj$relative_precision_as, relative_precision, along = 2))
  RSTr_obj$age_metadata$names <- colnames(RSTr_obj$medians_as)
  RSTr_obj$age_metadata$std_pop[[new_name]] <- std_pop
  RSTr_obj$age_metadata$groups[[new_name]] <- groups
  if (RSTr_obj$params$suppressed) RSTr_obj <- suppress_estimates(RSTr_obj, RSTr_obj$params$suppress_threshold)
  RSTr_obj
}

#' @noRd
erase_duplicates <- function(arr) {
  arr_groups <- which(!duplicated(dimnames(arr)[[2]], fromLast = TRUE))
  arr[, arr_groups, , drop = FALSE]
}