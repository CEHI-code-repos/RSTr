#' Gibbs sampler
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#'
#' @noRd
gibbs_u <- function(name, dir, iterations, .show_plots, .discard_burnin) {
  sampler_start <- Sys.time()
  data <- readRDS(paste0(dir, name, "/data.Rds"))
  Y    <- data$Y
  n    <- data$n
  miss <- which(!is.finite(Y))

  spatial_data <- readRDS(paste0(dir, name, "/spatial_data.Rds"))
  adjacency     <- spatial_data$adjacency
  num_adj       <- spatial_data$num_adj
  island_region <- spatial_data$island_region
  island_id     <- spatial_data$island_id
  num_island    <- spatial_data$num_island
  num_island_region <- sapply(island_region, length)

  priors <- readRDS(paste0(dir, name, "/priors.Rds"))
  tau_a    <- priors$tau_a
  tau_b    <- priors$tau_b
  sig_a    <- priors$sig_a
  sig_b    <- priors$sig_b
  theta_sd <- priors$theta_sd
  t_accept <- priors$t_accept

  inits <- readRDS(paste0(dir, name, "/inits.Rds"))
  theta <- inits$theta
  beta  <- inits$beta
  Z     <- inits$Z
  sig2  <- inits$sig2
  tau2  <- inits$tau2

  params <- readRDS(paste0(dir, name, "/params.Rds"))
  total       <- params$total
  m0          <- params$m0
  A           <- params$A
  rho_up      <- params$rho_up
  method      <- params$method
  impute_lb   <- params$impute_lb
  impute_ub   <- params$impute_ub
  start_batch <- params$batch
  batches     <- seq(start_batch + 1, start_batch + iterations / 100)
  num_region  <- length(Y)
  num_island  <- length(beta)

  cat("Starting sampler on Batch", start_batch + 1, "at", format(Sys.time(), "%a %b %d %X"), "\n")
  plots <- output <- vector("list", length(inits))
  names(plots) <- names(output) <- par_up <- names(inits)
  plot_its <- NULL
  for (batch in batches) {
    T_inc <- 100
    display_progress(batch, max(batches), total, 0, T_inc, sampler_start)
    output$theta <- array(dim = c(num_region, T_inc / 10))
    output$beta  <- array(dim = c(num_island, T_inc / 10))
    output$sig2  <- array(dim = c(T_inc / 10))
    output$tau2  <- array(dim = c(T_inc / 10))
    output$Z     <- array(dim = c(num_region, T_inc / 10))

    # Metropolis for Yikt
    t_accept <- ifelse(t_accept < 1 / 6, 1 / 6, ifelse(t_accept > 0.75, 0.75, t_accept))
    theta_sd <- ifelse(
      t_accept > 0.5,
      theta_sd * t_accept / 0.5,
      ifelse(t_accept < 0.35, theta_sd * t_accept / 0.35, theta_sd)
    )
    t_accept <- rep(0, num_region)
    for(it in 1:T_inc) {
      #### impute missing Y's ####
      if (length(miss)) {
        if (method == "binom") {
          rate <- expit(theta[miss])
          rp   <- stats::runif(
            length(miss),
            stats::pbinom(impute_lb - 0.1, round(n[miss]), rate),
            stats::pbinom(impute_ub + 0.1, round(n[miss]), rate)
          )
          Y[miss] <- stats::qbinom(rp, round(n[miss]), rate)
        }
        if (method == "pois") {
          rate <- n[miss] * exp(theta[miss])
          rp   <- stats::runif(
            length(miss),
            stats::ppois(impute_lb - 0.1, rate),
            stats::ppois(impute_ub + 0.1, rate)
          )
          Y[miss] <- stats::qpois(rp, rate)
        }
      }

      Z     <- update_Z_u(Z, sig2, theta, beta, tau2, adjacency, num_adj, island_region, island_id)
      sig2  <- update_sig2_u(sig2, Z, beta, tau2, adjacency, num_adj, island_region, num_island_region, A, m0, sig_a, sig_b, method)
      tau2  <- update_tau2_u(tau2, theta, beta, Z, sig2, num_island_region, island_id, A, m0, tau_a, tau_b, method)
      theta <- update_theta_u(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method)
      beta  <- update_beta_u(beta, theta, Z, tau2, sig2, A, m0, island_region, method)

      #### Save outputs ####
      if (it %% 10 == 0) {
        output$beta [, it / 10] <- beta
        output$sig2 [  it / 10] <- sig2
        output$tau2 [  it / 10] <- tau2
        output$theta[, it / 10] <- theta
        output$Z    [, it / 10] <- Z
      }
      display_progress(batch, max(batches), total, it, T_inc, sampler_start)
    }

    # modify meta-parameters, save outputs to respective files
    total <- total + T_inc
    t_accept <- t_accept / T_inc
    inits <- list(
      theta = theta,
      beta  = beta,
      Z     = Z,
      sig2  = sig2,
      tau2  = tau2
    )
    priors$theta_sd <- theta_sd
    priors$t_accept <- t_accept
    params$total    <- total
    params$batch    <- batch
    saveRDS(params, paste0(dir, name, "/params.Rds"))
    saveRDS(priors, paste0(dir, name, "/priors.Rds"))
    saveRDS(inits,  paste0(dir, name, "/inits.Rds"))
    save_output(output, batch, dir, name, .discard_burnin)
    if (.show_plots) {
      output_its  <- seq((batch - 1) * 100 + 10, batch * 100, 10)
      plot_its    <- c(plot_its, output_its)
      plots$beta  <- c(plots$beta,  output$beta [1, ])
      plots$theta <- c(plots$theta, output$theta[1, ])
      plots$Z     <- c(plots$Z,     output$Z    [1, ])
      plots$tau2  <- c(plots$tau2,  output$tau2)
      plots$sig2  <- c(plots$sig2,  output$sig2)
      grid <- c(2, 3)
      graphics::par(mfrow = grid)
      # Gradually remove plots in burn-in, then plot
      if (plot_its[1] < 2000) {
        plots    <- lapply(plots, \(par) par[-(1:5)])
        plot_its <- plot_its[-(1:5)]
      }
      lapply(
        par_up,
        \(par) plot(plot_its, plots[[par]], type = "l", main = par, xlab = "Iteration", ylab = "Value")
      )
    }

  }
  cat("\nModel finished at", format(Sys.time(), "%a %b %d %X"), "\n")
}
