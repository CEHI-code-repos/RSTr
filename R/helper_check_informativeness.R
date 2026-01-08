check_informativeness <- function(RSTr_obj) {
  sig2 = get_medians(load_samples(RSTr_obj, "sig2"))
  tau2 = get_medians(load_samples(RSTr_obj, "tau2"))
  a0 = stats::median(1 / (exp(sig2 + (sig2 + tau2) / 3) - 1))
  if (a0 > 6) {
    warning(
      "Model informativeness is too high. Re-run data with restricted CAR model using 'rcar()'. See vignette('rstr-informativeness') for more details"
    )
  }
}
