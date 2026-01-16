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
  perm <- c(margin, setdiff(seq_along(dim(sample)), margin))
  d <- dim(sample)[length(dim(sample))]
  m <- dim(sample)[margin]
  rest <- prod(dim(sample)[-c(margin, length(dim(sample)))])
  x <- sample |> aperm(perm) |> matrix(nrow = m, ncol = rest * d)
  w <- pop |> aperm(perm[-length(perm)]) |> matrix(nrow = m, ncol = rest)
  num <- colSums(x * rep(w, times = d), na.rm = TRUE)
  den <- rep(colSums(w, na.rm = TRUE), times = d)
  wtd_avg <- num / den
  if (bind_new) {
    new_dim <- dim(sample)
    new_dim[margin] <- 1
    new_dimnames <- dimnames(sample)
    new_dimnames[[margin]] <- c(new_dimnames[[margin]], new_name)
    agg_sample <- array(wtd_avg, dim = new_dim)
    array_new <- abind::abind(sample, agg_sample, along = margin)
    dimnames(array_new) <- new_dimnames
  } else {
    new_dim <- dim(sample)[-margin]
    new_dimnames <- dimnames(sample)[-2]
    array_new <- array(wtd_avg, dim = new_dim, dimnames = new_dimnames)
  }
  array_new
}
