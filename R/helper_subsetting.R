#' Extract all samples for first instance of parameter
#' @noRd
extract_last_margin <- function(arr) {
  idx <- c(rep(list(1), length(dim(arr)) - 1), list(seq_len(rev(dim(arr))[1])))
  do.call(`[`, c(list(arr), idx))
}

subset_array <- function(array, margin, groups) {
  index_list <- lapply(dim(array), \(x) 1:x)
  index_list[[margin]] <- groups
  do.call(`[`, c(list(array), index_list, drop = FALSE))
}