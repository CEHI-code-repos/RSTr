get_elapsed_time <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units = "secs")
  format(.POSIXct(dt, tz = "GMT"), "%H:%M:%OS")
}

display_progress <- function(batch, total_batches, total, it, sampler_start) {
  batch_ratio <- paste0(batch, "/", total_batches, ",")
  prog_done <- strrep("*", floor(it / 2))
  prog_notdone <- strrep(".", ceiling((100 - it) / 2))
  progress_bar <- paste0("|", prog_done, prog_notdone, "|")
  t <- get_elapsed_time(sampler_start)
  cat("Batch", batch_ratio, "Progress:", progress_bar, "Elapsed Time:", t, "\r")
  if (batch == total_batches && it == 100) {
    cat("\n")
  }
}
