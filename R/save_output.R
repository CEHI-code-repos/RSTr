#' Save output
#'
#' @noRd
save_output <- function(output, batch, dir, name, .discard_burnin) {
  if (!.discard_burnin | batch > 20) {
    for (par in names(output)) {
      saveRDS(output[[par]], paste0(dir, name, "/", par, "/", par, "_out_", batch, ".Rds"))
    }
  }
}