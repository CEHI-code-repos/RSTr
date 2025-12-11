get_spatial_data <- function(adjacency) {
  if (!inherits(adjacency, "nb")) {
    adjacency <- lapply(adjacency, as.integer)
    class(adjacency) <- c("nb")
  }
  check_regions_unlinked(adjacency)
  num_adj <- spdep::card(adjacency)
  island_region <- lapply(get_islands(adjacency), \(x) x - 1)
  num_island_region <- lengths(island_region)
  adjacency <- lapply(adjacency, \(x) x - 1)
  num_island <- length(island_region)
  island_id <- rep(NA, length(adjacency))
  for (isl in 1:num_island) {
    island_id[island_region[[isl]] + 1] <- isl - 1
  }
  list(
    adjacency = adjacency,
    num_adj = num_adj,
    island_region = island_region,
    num_island_region = num_island_region,
    island_id = island_id,
    num_island = num_island
  )
}

get_islands <- function(adjacency) {
  f <- seq_along(adjacency)
  island_region <- list()
  group <- 0
  while (length(f) > 0) {
    active_list <- f[1]
    inactive_list <- NULL
    while (length(active_list) > 0) {
      Na <- adjacency[[active_list[1]]]
      active_list <- unique(c(active_list, Na[which(!(Na %in% inactive_list))]))
      inactive_list <- c(inactive_list, active_list[1])
      active_list <- active_list[-1]
    }
    group <- group + 1
    inactive_list <- sort(inactive_list)
    island_region[[group]] <- inactive_list
    f <- setdiff(f, inactive_list)
  }
  island_region
}

check_regions_unlinked <- function(adjacency) {
  if (any(spdep::card(adjacency) == 0)) {
    stop("Some regions in 'adjacency' have no neighbors. Ensure all regions have at least 1 neighbor. Check vignette('RSTr-adjacency') for more information")
  }
}



