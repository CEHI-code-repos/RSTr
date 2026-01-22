#' Suppress estimates based on reliability criteria
#'
#' Generates suppressed estimates for an \code{RSTr} model object with a given relative precision and population/event threshold.
#'
#' While the \code{threshold} argument is optional, population/event thresholds are necessary for non-restricted models. Population/event thresholds should only be omitted for restricted CAR models, such as the RCAR.
#'
#' @param RSTr_obj An \code{RSTr} model object.
#' @param threshold The population/event suppression threshold.
#' @param type Determines whether suppression threshold is based on population counts or event counts.
#' @returns An \code{RSTr} model object with suppressed estimates.
#' @examples
#' std_pop <- c(113154, 100640, 95799)
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' on.exit(unlink(file.path(tempdir(), "test"), recursive = TRUE), add = TRUE)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
#' mod_mst <- suppress_estimates(mod_mst, threshold = 1000, type = "population")
#' estimates_table <- get_estimates(mod_mst)
#' @export
suppress_estimates <- function(
  RSTr_obj,
  threshold = 0,
  type = c("population", "event")
) {
  type <- match.arg(type)
  RSTr_obj$params$suppressed <- TRUE
  RSTr_obj$params$supp_thres <- threshold
  if (threshold == 0 && !(RSTr_obj$params$model %in% c("rcar"))) {
    warning(
      "Suppressing estimates without a population/event threshold is not recommended for non-restricted models. Specify `threshold` or re-run with restricted model"
    )
  }
  if (threshold > 0 && (RSTr_obj$params$model %in% c("rcar"))) {
    warning(
      "Suppressing estimates with a population/event threshold not necessary for restricted models"
    )
  }
  medians_suppressed <- RSTr_obj$medians
  supp <- (RSTr_obj$rel_prec < 1) | (RSTr_obj$data$n < threshold)
  medians_suppressed[supp] <- NA
  RSTr_obj$medians_suppressed <- medians_suppressed
  if (RSTr_obj$params$age_standardized) {
    medians_suppressed_as <- RSTr_obj$medians_as
    supp <- (RSTr_obj$rel_prec_as < 1) | (RSTr_obj$data_as$n < threshold)
    medians_suppressed_as[supp] <- NA
    RSTr_obj$medians_suppressed_as <- medians_suppressed_as
  }
  RSTr_obj
}
