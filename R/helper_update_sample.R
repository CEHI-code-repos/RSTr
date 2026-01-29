update_sample <- function(RSTr_obj) {
  UseMethod("update_sample")
}

#' @export
update_sample.car <- function(RSTr_obj) {
  update_lambda(RSTr_obj)
  update_Z(RSTr_obj)
  update_tau2(RSTr_obj)
  update_beta(RSTr_obj)
  update_sig2(RSTr_obj)
  RSTr_obj
}

#' @export
update_sample.mcar <- function(RSTr_obj) {
  update_lambda(RSTr_obj)
  update_Z(RSTr_obj)
  update_tau2(RSTr_obj)
  update_beta(RSTr_obj)
  update_G(RSTr_obj)
  RSTr_obj
}

#' @export
update_sample.mstcar <- function(RSTr_obj) {
  update_lambda(RSTr_obj)
  update_Z(RSTr_obj)
  update_tau2(RSTr_obj)
  update_beta(RSTr_obj)
  update_G(RSTr_obj)
  update_Ag(RSTr_obj)
  RSTr_obj
}

#' @export
update_sample.mstcar_update_rho <- function(RSTr_obj) {
  update_lambda(RSTr_obj)
  update_Z(RSTr_obj)
  update_tau2(RSTr_obj)
  update_beta(RSTr_obj)
  update_G(RSTr_obj)
  update_Ag(RSTr_obj)
  update_rho(RSTr_obj)
  RSTr_obj
}
