save_output <- function(output, batch, dir, name) {
  for (par in names(output)) saveRDS(output[[par]], file.path(dir, name, par, paste0(par, "_out_", batch, ".Rds")))
}

save_model <- function(RSTr_obj) {
  dir <- RSTr_obj$params$dir
  name <- RSTr_obj$params$name
  saveRDS(RSTr_obj, file.path(dir, name, paste0(name, ".Rds")))
}

create_model_directory <- function(name, dir, pars) {
  new_dir <- file.path(dir, name)
  if (!dir.exists(new_dir)) dir.create(new_dir)
  param_dirs <- file.path(dir, name, pars)
  for (par in param_dirs) if(!dir.exists(par)) dir.create(par)
}