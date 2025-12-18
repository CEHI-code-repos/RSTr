update_beta <- function(RSTr_obj) {
  UseMethod("update_beta")
}

#' @export
update_beta.default <- function(RSTr_obj) {
  update_beta_default(RSTr_obj)
}

#' @export
update_beta.rcar <- function(RSTr_obj) {
  update_beta_rcar(RSTr_obj)
}

update_G <- function(RSTr_obj) {
  UseMethod("update_G")
}

#' @export
update_G.default <- function(RSTr_obj) {
  update_G_default(RSTr_obj)
}

#' @export
update_G.mstcar <- function(RSTr_obj) {
  update_G_mstcar(RSTr_obj)
}

update_sig2 <- function(RSTr_obj) {
  UseMethod("update_sig2")
}

#' @export
update_sig2.default <- function(RSTr_obj) {
  update_sig2_default(RSTr_obj)
}

#' @export
update_sig2.rcar <- function(RSTr_obj) {
  update_sig2_rcar(RSTr_obj)
}

update_tau2 <- function(RSTr_obj) {
  UseMethod("update_tau2")
}

#' @export
update_tau2.default <- function(RSTr_obj) {
  update_tau2_default(RSTr_obj)
}

#' @export
update_tau2.rcar <- function(RSTr_obj) {
  update_tau2_rcar(RSTr_obj)
}

#' @export
update_tau2.mstcar <- function(RSTr_obj) {
  update_tau2_mstcar(RSTr_obj)
}

update_Z <- function(RSTr_obj) {
  UseMethod("update_Z")
}

#' @export
update_Z.car <- function(RSTr_obj) {
  update_Z_car(RSTr_obj)
}

#' @export
update_Z.mcar <- function(RSTr_obj) {
  update_Z_mcar(RSTr_obj)
}

#' @export
update_Z.mstcar <- function(RSTr_obj) {
  update_Z_mstcar(RSTr_obj)
}
