#' Age-standardize model objects
#'
#' Generate instructions for age-standardizing samples using a standard population for an \code{RSTr} model object.
#'
#' @param RSTr_obj An \code{RSTr} model object.
#' @param std_pop A vector of standard populations.
#' @param new_name The name to assign to the age-standardized group.
#' @param groups A vector of either indices for each group or a vector of strings for each group name. If set to \code{NULL}, will use all groups in the dataset.
#' @returns An \code{RSTr} object with age-standardized parameters.
#' @examples
#' std_pop <- c(113154, 100640, 95799)
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' on.exit(unlink(file.path(tempdir(), "test"), recursive = TRUE), add = TRUE)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
#' # age-standardize by all age groups
#' mod_mst <- age_standardize(mod_mst, std_pop, "35-64")
#' # Add onto age-standardized estimates. Age-standardize only by the first two age groups
#' mod_mst <- age_standardize(mod_mst, std_pop[1:2], "35-54", groups = 1:2)
#' @export
age_standardize <- function(RSTr_obj, std_pop, new_name, groups = NULL) {
  RSTr_obj$params$age_standardized <- TRUE
  RSTr_obj$as_data$std_pop[[new_name]] <- std_pop
  RSTr_obj$as_data$groups[[new_name]] <- groups %||% seq_along(std_pop)
  RSTr_obj$as_data$names <- c(RSTr_obj$as_data$names, new_name)
  RSTr_obj
}
