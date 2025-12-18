new_model <- function(model, data, restricted = NULL, update_rho = NULL) {
  data <- prepare_data(data)
  switch(
    model,
    car = new_car(data),
    rcar = new_rcar(data),
    mcar = new_mcar(data),
    mstcar = {
      if (update_rho) new_mstcar_update_rho(data) else new_mstcar(data)
    }
  )
}

new_RSTr <- function(data, subclass = character()) {
  structure(list(data = data), class = c(subclass, "RSTr"))
}

new_car <- function(data, subclass = character()) {
  new_RSTr(data, subclass = c(subclass, "car"))
}

new_rcar <- function(data) {
  new_car(data, subclass = "rcar")
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
