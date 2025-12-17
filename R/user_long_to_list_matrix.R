#' Generate count data for RSTr object
#'
#' \code{long_to_list_matrix()} converts a long \code{table} featuring event counts across regions and other optional margins into a \code{list} that is readable by \code{*car()}.
#'
#' \code{long_to_list_matrix()} will sum along any group/time stratifications that aren't specified; for example, if your dataset contains time periods and time is not specified in \code{long_to_list_matrix()}, the output will be a sum of all time periods. Filter data by desired groups and time periods before running \code{long_to_list_matrix()}.
#'
#' @param table A \code{table} containing event and mortality counts stratified by group/region/time.
#' @param event The column containing event counts.
#' @param population The column containing population counts.
#' @param region The column containing region names.
#' @param group An optional column containing sociodemographic group names.
#' @param time An optional column containing time period names.
#' @returns A \code{list} of mortality and population counts organized into multi-dimensional arrays.
#' @examples
#' ma_data <- maexample[!is.na(maexample$Year), ]
#' # Generates data from 1979-1981 stratified by sex
#' ma_data_mst <- long_to_list_matrix(ma_data, Deaths, Population, County.Code, Sex.Code, Year.Code)
#' ma_data_79 <- ma_data[ma_data$Year == 1979, ]
#' # Generates 1979 data stratified by sex
#' ma_data_m <- long_to_list_matrix(ma_data_79, Deaths, Population, County.Code, Sex.Code)
#' # Generates 1979 data summarized for all sexes
#' ma_data_u <- long_to_list_matrix(ma_data_79, Deaths, Population, County.Code)
#' @export
long_to_list_matrix <- function(
  table,
  event,
  population,
  region,
  group = NULL,
  time = NULL
) {
  ev <- deparse(substitute(event))
  po <- deparse(substitute(population))
  re <- deparse(substitute(region))
  gr <- deparse(substitute(group))
  ti <- deparse(substitute(time))
  formula_event <- stats::reformulate(c(re, gr, ti), response = ev)
  formula_population <- stats::reformulate(c(re, gr, ti), response = po)
  list(
    Y = stats::xtabs(formula_event, table),
    n = stats::xtabs(formula_population, table)
  )
}
