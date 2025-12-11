tune_metropolis_sd <- function(sd, accept) {
  accept <- pmin(pmax(accept, 1 / 6), 0.75)
  sd <- ifelse(accept > 0.5, sd * accept / 0.5, sd)
  sd <- ifelse(accept < 0.35, sd * accept / 0.35, sd)
  sd
}

#' Update metropolis SD inside of priors
#' @noRd
update_priors_sd <- function(RSTr_obj) {
  UseMethod("update_priors_sd")
}

#' @export
update_priors_sd.default <- function(RSTr_obj) {
  priors <- RSTr_obj$priors
  priors$lambda_sd <- tune_metropolis_sd(priors$lambda_sd, priors$lambda_accept / 100)
  priors$lambda_accept[] <- 0
  RSTr_obj$priors <- priors
  RSTr_obj
}

#' @export
update_priors_sd.mstcar_update_rho <- function(RSTr_obj) {
  priors <- RSTr_obj$priors
  priors$lambda_sd <- tune_metropolis_sd(priors$lambda_sd, priors$lambda_accept / 100)
  priors$lambda_accept[] <- 0
  priors$rho_sd <- tune_metropolis_sd(priors$rho_sd, priors$rho_accept / 100)
  priors$rho_accept[] <- 0
  RSTr_obj$priors <- priors
  RSTr_obj
}