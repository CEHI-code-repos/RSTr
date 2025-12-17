#' Aggregate count arrays
#'
#' Sums counts over event/population arrays. Useful when manually generating group-aggregated/age-standardized estimates and a population threshold is needed for suppression.
#'
#' @inheritParams standardize_samples
#' @param count The \code{array} to aggregate.
#' @returns An \code{array} of aggregated count data.
#' @examples
#' margin_time <- 3
#' # aggregate population from all years for each county-group
#' pop_7988 <- aggregate_count(miheart$n, margin_time)
#' # aggregate population from 1980-1984 for each county-group
#' pop_8084 <- aggregate_count(miheart$n, margin_time, groups = as.character(1980:1984))
#' # bind aggregated pop from all years to population data
#' pop_agg <- aggregate_count(miheart$n, margin_time, bind_new = TRUE, new_name = "1979-1988")
#' @export
aggregate_count <- function(
  count,
  margin,
  groups = NULL,
  bind_new = FALSE,
  new_name = NULL
) {
  mar <- seq_along(dim(count))[-margin]
  sub_count <- count
  if (!is.null(groups)) {
    sub_count <- subset_array(count, margin, groups)
  }
  if (bind_new) {
    new_dim <- dim(count)
    new_dim[margin] <- 1
    new_count <- array(apply(sub_count, mar, sum, na.rm = TRUE), dim = new_dim)
    newnames <- c(dimnames(count)[[margin]], new_name)
    array_agg <- abind::abind(count, new_count, along = margin)
    dimnames(array_agg)[[margin]] <- newnames
  } else {
    array_agg <- apply(sub_count, mar, sum, na.rm = TRUE)
  }
  array_agg
}
