#' Age-standardize samples
#' 
#' Age-standardizes samples using a standard population.
#' 
#' @param sample an \code{array} of samples imported with \code{load_samples()}
#' @param std_pop A vector of standard populations.
#' @param margin For \code{array}s, The margin on which the groups of interest are stratified.
#' @param groups A vector of either indices for each group or a vector of strings for each group name. If set to \code{NULL}, will use all groups in the dataset.
#' @param bind_new If set to \code{TRUE}, will bind an \code{array} to the original sample dataset. Otherwise, will generate a standalone array of samples.
#' @param new_name The name to assign to the age-standardized group.
#' @returns An \code{array} of age-standardized samples.
#' @examples
#' std_pop <- c(113154, 100640, 95799)
#' age_margin <- 2
#' # age-standardize by all age groups
#' samples_3564 <- standardize_samples(minsample, std_pop, age_margin)
#' # age-standardize only by the first two age groups
#' samples_3554 <- standardize_samples(minsample, std_pop[1:2], age_margin, groups = 1:2)
#' # bind age-standardized samples to original samples
#' samples_as <- standardize_samples(
#'   minsample,
#'   std_pop,
#'   age_margin,
#'   bind_new = TRUE,
#'   new_name = "35-64"
#' )
#' @export
standardize_samples <- function(sample, std_pop, margin, groups = NULL, bind_new = FALSE, new_name = NULL) {
  mar <- seq_along(dim(sample))[-margin]
  wts <- std_pop / sum(std_pop)
  sub_sample <- sample
  if (!is.null(groups)) {
    sub_sample <- subset_array(sample, margin, groups)
  }
  if (bind_new) {
    new_dim <- dim(sample)
    new_dim[margin] <- 1
    agg_sample <- array(
      apply(sweep(sub_sample, margin, wts, "*"), mar, sum, na.rm = TRUE),
      dim <- new_dim
    )
    array_new <- abind::abind(sample, agg_sample, along = margin)
    newnames <- c(dimnames(sample)[[margin]], new_name)
    dimnames(array_new)[[margin]] <- newnames
  } else {
    array_new <- apply(sweep(sub_sample, margin, wts, "*"), mar, sum, na.rm = TRUE)
  }
  array_new
}