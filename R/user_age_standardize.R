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
  if (is.null(groups)) {
    groups <- seq_len(dim(samples)[2])
  }
  data <- RSTr_obj$data |>
    lapply(aggregate_count, 2, groups, TRUE, new_name) |>
    lapply(\(x) x[, new_name, , drop = FALSE])
  samples <- samples |>
    subset_array(2, groups) |>
    standardize_samples(std_pop, 2, groups, TRUE, new_name) |>
    _[, new_name, , , drop = FALSE]
  medians <- get_medians(samples)
  ci <- get_credible_interval(samples)
  rel_prec <- get_relative_precision(medians, ci)
  RSTr_obj$medians_as <- bind_objects(RSTr_obj$medians_as, medians)
  RSTr_obj$data_as$Y <- bind_objects(RSTr_obj$data_as$Y, data$Y)
  RSTr_obj$data_as$n <- bind_objects(RSTr_obj$data_as$n, data$n)
  RSTr_obj$ci_as$lower <- bind_objects(RSTr_obj$ci_as$lower, ci$lower)
  RSTr_obj$ci_as$upper <- bind_objects(RSTr_obj$ci_as$upper, ci$upper)
  RSTr_obj$rel_prec_as <- bind_objects(RSTr_obj$rel_prec_as, rel_prec)
  RSTr_obj$as_data$names <- colnames(RSTr_obj$medians_as)
  RSTr_obj$as_data$std_pop[[new_name]] <- std_pop
  RSTr_obj$as_data$groups[[new_name]] <- groups
  if (RSTr_obj$params$suppressed) {
    RSTr_obj <- suppress_estimates(RSTr_obj, RSTr_obj$params$supp_thres)
  }
  RSTr_obj
}

erase_duplicates <- function(arr) {
  arr_groups <- which(!duplicated(dimnames(arr)[[2]], fromLast = TRUE))
  arr[, arr_groups, , drop = FALSE]
}

bind_objects <- function(obj, obj_new) {
  obj |> abind::abind(obj_new, along = 2) |> erase_duplicates()
}
