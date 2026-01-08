get_sp_data <- function(adjacency) {
  if (!inherits(adjacency, "nb")) {
    adjacency <- lapply(adjacency, as.integer)
    class(adjacency) <- c("nb")
  }
  check_regions_unlinked(adjacency)
  comp <- spdep::n.comp.nb(adjacency)
  n_adj <- spdep::card(adjacency)
  n_island <- comp$nc
  isl_id <- comp$comp.id
  isl_region <- lapply(1:n_island, \(isl) which(isl_id == isl))
  n_isl_region <- lengths(isl_region)
  list(
    adjacency = adjacency,
    n_adj = n_adj,
    isl_region = isl_region,
    n_isl_region = n_isl_region,
    isl_id = isl_id,
    n_island = n_island
  )
}

convert_index <- function(RSTr_obj, index = c("zero", "one")) {
  match.arg(index)
  adjacency <- RSTr_obj$sp_data$adjacency
  isl_id <- RSTr_obj$sp_data$isl_id
  isl_region <- RSTr_obj$sp_data$isl_region
  if (index == "zero") {
    adjacency <- lapply(adjacency, \(x) x - 1)
    isl_region <- lapply(isl_region, \(x) x - 1)
    isl_id <- isl_id - 1
  }
  if (index == "one") {
    adjacency <- lapply(adjacency, \(x) x + 1)
    isl_region <- lapply(isl_region, \(x) x + 1)
    isl_id <- isl_id + 1
    as_nb(adjacency)
  }
  RSTr_obj$sp_data$adjacency <- adjacency
  RSTr_obj$sp_data$isl_id <- isl_id
  RSTr_obj$sp_data$isl_region <- isl_region
  RSTr_obj
}

check_regions_unlinked <- function(adjacency) {
  if (any(spdep::card(adjacency) == 0)) {
    stop(
      "Some regions in 'adjacency' have no neighbors. Ensure all regions have at least 1 neighbor. Check vignette('RSTr-adjacency') for more information"
    )
  }
}
