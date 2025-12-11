get_spatial_data <- function(adjacency) {
  if (!inherits(adjacency, "nb")) {
    adjacency <- lapply(adjacency, as.integer)
    class(adjacency) <- c("nb")
  }
  check_regions_unlinked(adjacency)
  comp <- spdep::n.comp.nb(adjacency)
  num_adj <- spdep::card(adjacency)
  num_island <- comp$nc
  island_id <- comp$comp.id
  island_region <- lapply(1:num_island, \(isl) which(island_id == isl))
  num_island_region <- lengths(island_region)
  list(
    adjacency = adjacency,
    num_adj = num_adj,
    island_region = island_region,
    num_island_region = num_island_region,
    island_id = island_id,
    num_island = num_island
  )
}

convert_index <- function(RSTr_obj, index = c("zero", "one")) {
  if (index == "zero") {
    RSTr_obj$spatial_data$adjacency <- lapply(RSTr_obj$spatial_data$adjacency, \(x) x - 1)
    RSTr_obj$spatial_data$island_region <- lapply(RSTr_obj$spatial_data$island_region, \(x) x - 1)
    RSTr_obj$spatial_data$island_id <- RSTr_obj$spatial_data$island_id - 1
  }
  if (index == "one") {
    RSTr_obj$spatial_data$adjacency <- lapply(RSTr_obj$spatial_data$adjacency, \(x) x + 1)
    RSTr_obj$spatial_data$island_region <- lapply(RSTr_obj$spatial_data$island_region, \(x) x + 1)
    RSTr_obj$spatial_data$island_id <- RSTr_obj$spatial_data$island_id + 1
    as_nb(RSTr_obj$spatial_data$adjacency)
  }
  RSTr_obj
}

check_regions_unlinked <- function(adjacency) {
  if (any(spdep::card(adjacency) == 0)) {
    stop("Some regions in 'adjacency' have no neighbors. Ensure all regions have at least 1 neighbor. Check vignette('RSTr-adjacency') for more information")
  }
}



