#' Gibbs sampler
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#'
#' @noRd
gibbs_m <- function(name, dir, iterations, .show_plots, .discard_burnin) {
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

  priors <- readRDS(paste0(dir, name, "/priors.Rds"))
  G_scale  <- priors$G_scale
  G_df     <- priors$G_df
  tau_a    <- priors$tau_a
  tau_b    <- priors$tau_b
  theta_sd <- priors$theta_sd
  t_accept <- priors$t_accept

  inits <- readRDS(paste0(dir, name, "/inits.Rds"))
  theta <- inits$theta
  beta  <- inits$beta
  Z     <- inits$Z
  G     <- inits$G
  tau2  <- inits$tau2

  params <- readRDS(paste0(dir, name, "/params.Rds"))
  total       <- params$total
  method      <- params$method
  impute_lb   <- params$impute_lb
  impute_ub   <- params$impute_ub
  start_batch <- params$batch
  batches     <- seq(start_batch + 1, start_batch + iterations / 100)

  cat("Starting sampler on Batch", start_batch + 1, "at", format(Sys.time(), "%a %b %d %X"), "\n")
  plots <- output <- vector("list", length(inits))
  names(plots) <- names(output) <- par_up <- names(inits)
  plot_its <- NULL
  for (batch in batches) {
    T_inc <- 100
    display_progress(batch, max(batches), total, 0, T_inc, sampler_start)
    output$theta <- array(dim = c(dim(theta),   T_inc / 10))
    output$beta  <- array(dim = c(dim(beta),    T_inc / 10))
    output$G     <- array(dim = c(dim(G),       T_inc / 10))
    output$tau2  <- array(dim = c(length(tau2), T_inc / 10))
    output$Z     <- array(dim = c(dim(Z),       T_inc / 10))

    # Metropolis for Yikt
    t_accept <- ifelse(t_accept < 1 / 6, 1 / 6, ifelse(t_accept > 0.75, 0.75, t_accept))
    theta_sd <- ifelse(
      t_accept > 0.5,
      theta_sd * t_accept / 0.5,
      ifelse(t_accept < 0.35, theta_sd * t_accept / 0.35, theta_sd)
    )
    t_accept <- array(0, dim = dim(theta))
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

      #### Update parameters ####
      beta  <- update_beta_m(beta, theta, Z, tau2, island_region)
      Z     <- update_Z_m(Z, G, theta, beta, tau2, adjacency, num_adj, island_region, island_id)
      G     <- update_G_m(G, Z, G_df, G_scale, adjacency, num_island)
      tau2  <- update_tau2_m(tau2, theta, beta, Z, tau_a, tau_b, island_id)
      theta <- update_theta_m(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method)

      #### Save outputs ####
      if (it %% 10 == 0) {
        output$beta [, , it / 10] <- beta
        output$G    [, , it / 10] <- G
        output$tau2 [,   it / 10] <- tau2
        output$theta[, , it / 10] <- theta
        output$Z    [, , it / 10] <- Z
      }
      display_progress(batch, max(batches), total, it, T_inc, sampler_start)
    }

    # modify meta-parameters, save outputs to respective files
    total    <- total + T_inc
    t_accept <- t_accept / T_inc
    inits    <- list(
      theta = theta,
      beta  = beta,
      Z     = Z,
      G     = G,
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
      plots$beta  <- c(plots$beta,  output$beta [1, 1, ])
      plots$theta <- c(plots$theta, output$theta[1, 1, ])
      plots$Z     <- c(plots$Z,     output$Z    [1, 1, ])
      plots$tau2  <- c(plots$tau2,  output$tau2 [1,    ])
      plots$G     <- c(plots$G,     output$G    [1, 1, ])
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
