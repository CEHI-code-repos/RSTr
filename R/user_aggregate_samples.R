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
  sub_sample <- sample
  subpop <- pop
  if (!is.null(groups)) {
    sub_sample <- subset_array(sample, margin, groups)
    subpop <- subset_array(pop, margin, groups)
  }
  perm <- c(margin, setdiff(seq_along(dim(sub_sample)), margin))
  rest <- prod(dim(sub_sample)[-c(margin, length(dim(sub_sample)))])
  its <- dim(sub_sample)[length(dim(sub_sample))]
  ng <- dim(sub_sample)[margin]
  sub_sample <- arr_to_matrix(sub_sample, perm, ng, rest * its)
  subpop <- arr_to_matrix(subpop, perm[-length(perm)], ng, rest)
  num <- (sub_sample * rep(subpop, times = its)) |> colSums(na.rm = TRUE)
  den <- subpop |> colSums(na.rm = TRUE) |> rep(times = its)
  (num / den) |>
    create_array_new(sample, margin, bind_new, new_name)
}
