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
  perm <- c(margin, setdiff(seq_along(dim(sub_count)), margin))
  rest <- prod(dim(sub_count)[-margin])
  ng <- dim(sub_count)[margin]
  sub_count |>
    arr_to_matrix(perm, ng, rest) |>
    colSums(na.rm = TRUE) |>
    create_array_new(count, margin, bind_new, new_name)
}
