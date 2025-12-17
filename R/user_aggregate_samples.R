#' Aggregate samples by non-age group
#'
#' Consolidates a set of samples over non-age groups using a population array to create weighted-average samples.
#'
#' \code{aggregate_samples()} is only meant for non-age group data, such as spatial regions, time periods, or other sociodemographic groups (race, sex, etc.). If you are interested in consolidating samples by age group, use \code{age_standardize()} instead. Additionally, if you plan on doing age-standardization along with aggregating by other groups, always aggregate groups first before doing age-standardization to ensure that the samples are properly standardized.
#' @inheritParams standardize_samples
#' @param pop The population array to be used for weighted averages.
#' @returns An \code{array} of weighted-average samples.
#' @examples
#' pop <- miheart$n[1:2, 1:3, 1:3]
#' time_margin <- 3
#' # calculate prevalence by aggregating over time periods
#' samples_3564 <- aggregate_samples(minsample, pop, margin = time_margin)
#' # calculate prevalence of only the first two time periods
#' samples_3554 <- aggregate_samples(minsample, pop, time_margin, groups = 1:2)
#' # bind prevalence samples to original samples
#' samples_prev <- aggregate_samples(
#'   minsample,
#'   pop,
#'   time_margin,
#'   bind_new = TRUE,
#'   new_name = "1979-1981"
#' )
#' @export
aggregate_samples <- function(
  sample,
  pop,
  margin,
  groups = NULL,
  bind_new = FALSE,
  new_name = NULL
) {
  mar <- seq_along(dim(sample))[-margin]
  pop_arr <- array(pop, dim = c(dim(pop), rev(dim(sample))[1]))
  sub_sample <- sample
  sub_pop_arr <- pop_arr
  if (!is.null(groups)) {
    sub_sample <- subset_array(sample, margin, groups)
    sub_pop_arr <- subset_array(pop_arr, margin, groups)
  }
  if (bind_new) {
    new_dim <- dim(sample)
    new_dim[margin] <- 1
    agg_sample <- array(
      apply(sub_sample * sub_pop_arr, mar, sum, na.rm = TRUE) /
        apply(sub_pop_arr, mar, sum, na.rm = TRUE),
      dim <- new_dim
    )
    array_new <- abind::abind(sample, agg_sample, along = margin)
    newnames <- c(dimnames(sample)[[margin]], new_name)
    dimnames(array_new)[[margin]] <- newnames
  } else {
    array_new <- apply(sub_sample * sub_pop_arr, mar, sum, na.rm = TRUE) /
      apply(sub_pop_arr, mar, sum, na.rm = TRUE)
  }
  array_new
}
