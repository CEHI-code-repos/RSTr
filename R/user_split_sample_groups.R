#' Split sample groups
#'
#' Sequesters stratified sociodemographic group margin into individual array margins.
#'
#' When using \code{aggregate_samples()} or \code{standardize_samples()}, the group/age margin must only feature groups of similar type, e.g., you cannot age-standardize with groups that specify both age and race. \code{split_sample_groups()} sequesters each category of group into its own margin to allow group-aggregation and age-standardization of these multiply-stratified groups. Ensure that the delimiter character is only used to split groups. E.g., for an age-sex group named \code{35-64_m}, \code{"_"} will split the margins with names \code{"35-64"} and \code{"m"}, whereas for a group named \code{35_64_m}, \code{split_sample_group()} will fail.
#'
#' @param sample an \code{array} of samples imported with \code{load_samples()}
#' @param new_groups A string vector of names for each new group.
#' @param delimiter A character that specifies the break between group categories.
#' @returns An \code{array} of samples with separate margins for stratified groups.
#' @examples
#' dimnames(minsplit)[2] # Can't age-standardize due to age-sex stratification
#' new_groups = c("age", "sex")
#' delimiter = "_"
#' sample_split <- split_sample_groups(minsplit, new_groups, delimiter)
#' dimnames(sample_split)[2:3] # can now age-standardize
#' std_pop <- c(113154, 100640, 95799)
#' age_margin <- 2
#' sample_as <- standardize_samples(sample_split, std_pop, age_margin)
#' @export
split_sample_groups <- function(sample, new_groups, delimiter = "_") {
  ng <- length(new_groups)
  nm <- length(dim(sample)) + ng - 1
  resid_mar <- (ng + 2):nm
  group_dims <- dimnames(sample)[[2]] |>
    strsplit(delimiter) |>
    simplify2array() |>
    apply(1, unique) |>
    stats::setNames(new_groups) |>
    rev()
  new_dimnames <- c(
    dimnames(sample)[1],
    group_dims,
    dimnames(sample)[resid_mar - 2]
  )

  sample_new <- array(
    sample,
    dim = sapply(new_dimnames, length),
    dimnames = new_dimnames
  ) |>
    aperm(c(1, 2 + rev(0:(ng - 1)), resid_mar))
  sample_new
}
