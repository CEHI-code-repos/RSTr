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
#' on.exit(unlink(file.path(tempdir(), "test"), recursive = TRUE), add = TRUE)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
#' estimates_table <- get_estimates(mod_mst)
#' mod_mst <- age_standardize(mod_mst, std_pop, "35-64")
#' estimates_table_as <- get_estimates(mod_mst)
#' @export
get_estimates <- function(RSTr_obj, rates_per = 1e5, standardized = TRUE) {
  samples <- load_samples(RSTr_obj) * rates_per
  sample_data <- RSTr_obj$data
  if (RSTr_obj$params$age_standardized && standardized) {
    samples <- get_as_samples(samples, RSTr_obj)
    sample_data <- get_as_data(sample_data, RSTr_obj)
  }
  medians <- get_medians(samples)
  ci <- get_credible_interval(samples)
  rp <- get_relative_precision(medians, ci)

  medians_table <- get_medians_table(medians, rp, sample_data, RSTr_obj)

  ci_table_lower <- as.data.frame.table(ci$lower, responseName = "ci_lower")
  ci_table_upper <- as.data.frame.table(ci$upper, responseName = "ci_upper")
  ci_table <- merge(ci_table_lower, ci_table_upper)

  rp_table <- as.data.frame.table(rp, responseName = "rel_prec")

  events_table <- as.data.frame.table(sample_data$Y, responseName = "events")
  population_table <- as.data.frame.table(
    sample_data$n,
    responseName = "population"
  )
  data_table <- merge(events_table, population_table)

  medians_table |>
    merge(ci_table) |>
    merge(rp_table) |>
    merge(data_table)
}

get_as_samples <- function(samples, RSTr_obj) {
  m_age <- 2
  for (new_name in RSTr_obj$as_data$names) {
    stdpop <- RSTr_obj$as_data$std_pop[[new_name]]
    grps <- RSTr_obj$as_data$groups[[new_name]] %||% seq_along(stdpop)
    samples <- standardize_samples(
      samples,
      stdpop,
      m_age,
      grps,
      bind_new = TRUE,
      new_name = new_name
    )
  }
  subset_array(samples, m_age, RSTr_obj$as_data$names)
}
get_as_data <- function(sample_data, RSTr_obj) {
  m_age <- 2
  for (new_name in RSTr_obj$as_data$names) {
    grps <- RSTr_obj$as_data$groups
    sample_data <- lapply(
      sample_data,
      \(x) aggregate_count(x, m_age, grps[[new_name]], TRUE, new_name)
    )
  }
  lapply(sample_data, subset_array, m_age, RSTr_obj$as_data$names)
}
get_medians_table <- function(medians, rp, sample_data, RSTr_obj) {
  medians_table <- as.data.frame.table(medians, responseName = "medians")
  if (RSTr_obj$params$suppressed) {
    medians_suppressed <- medians
    type <- RSTr_obj$params$supp_type %||% "population"
    type <- ifelse(type == "event", "Y", "n")
    supp_thres <- sample_data[[type]] < RSTr_obj$params$supp_thres
    medians_suppressed[(rp < 1) | supp_thres] <- NA
    supp_table <- as.data.frame.table(
      medians_suppressed,
      responseName = "medians_suppressed"
    )
    medians_table <- merge(medians_table, supp_table)
  }
  medians_table
}
