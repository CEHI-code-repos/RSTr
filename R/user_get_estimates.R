#' Extract estimates from RSTr model object
#' 
#' Gathers model and estimate information for an \code{RSTr} model object, exported as a long table. Estimate rates and their respective credible intervals are displayed by default in rates per 100,000.
#' 
#' @param RSTr_obj An \code{RSTr} model object.
#' @param rates_per The desired scaling for estimate rates.
#' @param standardized If \code{RSTr_obj} contains age-standardized rates, shows the age-standardized rates. If set to \code{FALSE}, always shows the non-age-standardized rates.
#' @returns A long \code{table} containing region/group/time period names, estimates, credible intervals, relative precisions, and the associated event/population counts.
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
get_estimates <- function(RSTr_obj, rates_per = 1e5, standardized = TRUE) {
  marnames <- names(RSTr_obj$params$dimnames)
  if (is.null(marnames)) marnames <- c("region", "group", "time")
  marnames[!nzchar(marnames)] <- c("region", "group", "time")[!nzchar(marnames)]
  if (RSTr_obj$params$age_standardized && standardized) {
    est_table <- stats::setNames(as.data.frame.table(RSTr_obj$medians_as * rates_per), c(marnames, "medians"))
    if (RSTr_obj$params$suppressed) est_table$medians_suppressed <- c(RSTr_obj$medians_suppressed_as * rates_per)
    est_table$credible_interval_lower <- c(RSTr_obj$credible_interval_as$lower * rates_per)
    est_table$credible_interval_upper <- c(RSTr_obj$credible_interval_as$upper * rates_per)
    est_table$relative_precision <- c(RSTr_obj$relative_precision_as)
    est_table$events <- c(RSTr_obj$data_as$Y)
    est_table$population <- c(RSTr_obj$data_as$n)
  } else {
    est_table <- stats::setNames(as.data.frame.table(RSTr_obj$medians * rates_per), c(marnames, "medians"))
    if (RSTr_obj$params$suppressed) est_table$medians_suppressed <- c(RSTr_obj$medians_suppressed * rates_per)
    est_table$credible_interval_lower <- c(RSTr_obj$credible_interval$lower * rates_per)
    est_table$credible_interval_upper <- c(RSTr_obj$credible_interval$upper * rates_per)
    est_table$relative_precision <- c(RSTr_obj$relative_precision)
    est_table$events <- c(RSTr_obj$data$Y)
    est_table$population <- c(RSTr_obj$data$n)
  }
  na_test <- which(apply(est_table[, 1:3], 2, \(col) all(is.na(col))))
  if (length(na_test) > 0) est_table <- est_table[, -na_test]
  est_table
}
