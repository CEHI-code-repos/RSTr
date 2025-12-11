create_new_model <- function(model, data, restricted = NULL, update_rho = NULL) {
  data <- prepare_data(data)
  switch(
    model,
    ucar = new_ucar(data),
    eucar = new_ucar_restricted(data),
    mcar = new_mcar(data),
    mstcar = {
      if (update_rho) new_mstcar_update_rho(data) else new_mstcar(data)
    }
  )
}

new_RSTr <- function(data, subclass = character()) {
  structure(list(data = data), class = c(subclass, "RSTr"))
}

new_ucar <- function(data, subclass = character()) {
  new_RSTr(data, subclass = c(subclass, "ucar"))
}

new_ucar_restricted <- function(data) {
  new_ucar(data, subclass = "eucar")
}

new_mcar <- function(data, subclass = character()) {
  new_RSTr(data, subclass = c(subclass, "mcar"))
}

new_mstcar <- function(data, subclass = character()) {
  new_RSTr(data, subclass = c(subclass, "mstcar"))
}

new_mstcar_update_rho <- function(data) {
  new_mstcar(data, subclass = "mstcar_update_rho")
}

