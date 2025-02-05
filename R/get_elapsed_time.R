#' Get elapsed time
#' @noRd
get_elapsed_time <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  format(.POSIXct(dt, tz = "GMT"), "%H:%M:%OS")
}

#' Display time
#' @noRd
display_progress <- function(batch, total_batches, total, it, T_inc, sampler_start) {
  cat(
    "Batch", paste0(batch, "/", total_batches, ","),
    "Progress:", paste0("|", paste0(rep("*", floor(it / 2)), collapse = ""), paste0(rep(".", ceiling((T_inc - it) / 2)), collapse = ""), "|"), 
    "Elapsed Time:", get_elapsed_time(sampler_start),
    "\r"
  )
}
