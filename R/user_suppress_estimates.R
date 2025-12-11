#' Suppress estimates based on reliability criteria
#' 
#' Generates suppressed estimates for an \code{RSTr} model object with a given relative precision and population/event threshold.
#' 
#' While the \code{threshold} argument is optional, population/event thresholds are necessary for non-enhanced models. Population/event thresholds should only be omitted for enhanced CAR models, such as the EUCAR.
#' 
#' @param RSTr_obj An \code{RSTr} model object.
#' @param threshold The population/event suppression threshold.
#' @param type Determines whether suppression threshold is based on population counts or event counts.
#' @returns An \code{RSTr} model object with suppressed estimates.
#' @examples
#' std_pop <- c(113154, 100640, 95799)
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
#' estimates_table <- get_estimates(mod_mst)
#' mod_mst <- age_standardize(mod_mst, std_pop, "35-64")
#' estimates_table_as <- get_estimates(mod_mst)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
suppress_estimates <- function(RSTr_obj, threshold = 0, type = c("population", "event")) {
  type <- match.arg(type)
  RSTr_obj$params$suppressed <- TRUE
  RSTr_obj$params$suppress_threshold <- threshold
  if (threshold == 0 && !(RSTr_obj$params$model %in% c("eucar"))) {
    warning("Suppressing estimates without a population/event threshold is not recommended for non-enhanced models. Specify `threshold` or re-run with enhanced model")
  }
  if (threshold > 0 && (RSTr_obj$params$model %in% c("eucar"))) {
    warning("Suppressing estimates with a population/event threshold not necessary for EUCAR models")
  }
  medians_suppressed <- RSTr_obj$medians
  medians_suppressed[(RSTr_obj$relative_precision < 1) | (RSTr_obj$data$n < threshold)] <- NA
  RSTr_obj$medians_suppressed <- medians_suppressed
  if (RSTr_obj$params$age_standardized) {
    medians_suppressed_as <- RSTr_obj$medians_as
    medians_suppressed_as[(RSTr_obj$relative_precision_as < 1) | (RSTr_obj$data_as$n < threshold)] <- NA
    RSTr_obj$medians_suppressed_as <- medians_suppressed_as
  }
  RSTr_obj
}