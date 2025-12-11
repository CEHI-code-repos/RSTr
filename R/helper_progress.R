get_elapsed_time <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units = "secs")
  format(.POSIXct(dt, tz = "GMT"), "%H:%M:%OS")
}

display_progress <- function(batch, total_batches, total, it, sampler_start) {
  cat(
    "Batch", paste0(batch, "/", total_batches, ","),
    "Progress:", paste0("|", strrep("*", floor(it / 2)), strrep(".", ceiling((100 - it) / 2)), "|"),
    "Elapsed Time:", get_elapsed_time(sampler_start),
    "\r"
  )
  if (batch == total_batches && it == 100) {
    cat("\n")
  }
}
