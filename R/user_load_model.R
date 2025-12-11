#' Load model
#' 
#' \code{load_model()} imports an \code{RSTr} object with name \code{name} in directory \code{dir}.
#' 
#' @param name The name of the model to load.
#' @param dir The directory in which the model lives.
#' @returns An \code{RSTr} model object.
#' @examples
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
#' mod_mst <- load_model(name = "test", dir = tempdir())
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
load_model <- function(name, dir = tempdir()) {
  RSTr_obj <- readRDS(file.path(dir, name, paste0(name, ".Rds")))
  RSTr_obj$params$name <- name
  RSTr_obj$params$dir <- dir
  RSTr_obj
}