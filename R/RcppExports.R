# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

get_regs <- function(arr, ind) {
    .Call(`_RSTr_get_regs`, arr, ind)
}

get_grp <- function(arr, reg, time) {
    .Call(`_RSTr_get_grp`, arr, reg, time)
}

Sig_eta_i <- function(G, rho) {
    .Call(`_RSTr_Sig_eta_i`, G, rho)
}

Sig_eta <- function(Sein) {
    .Call(`_RSTr_Sig_eta`, Sein)
}

cpp_rmvnorm <- function(mean, covar) {
    .Call(`_RSTr_cpp_rmvnorm`, mean, covar)
}

geteig <- function(covar) {
    .Call(`_RSTr_geteig`, covar)
}

update_beta_m <- function(beta, theta, Z, tau2, island_region) {
    .Call(`_RSTr_update_beta_m`, beta, theta, Z, tau2, island_region)
}

update_Z_m <- function(Z, G, theta, beta, tau2, adjacency, num_adj, island_region, island_id) {
    .Call(`_RSTr_update_Z_m`, Z, G, theta, beta, tau2, adjacency, num_adj, island_region, island_id)
}

update_G_m <- function(G, Z, G_df, G_scale, adjacency, num_island) {
    .Call(`_RSTr_update_G_m`, G, Z, G_df, G_scale, adjacency, num_island)
}

update_tau2_m <- function(tau2, theta, beta, Z, tau_a, tau_b, island_id) {
    .Call(`_RSTr_update_tau2_m`, tau2, theta, beta, Z, tau_a, tau_b, island_id)
}

update_theta_m <- function(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method) {
    .Call(`_RSTr_update_theta_m`, theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method)
}

update_beta_mst <- function(beta, theta, Z, tau2, island_region) {
    .Call(`_RSTr_update_beta_mst`, beta, theta, Z, tau2, island_region)
}

update_Z_mst <- function(Z, G, theta, beta, rho, tau2, adjacency, num_adj, island_region, island_id) {
    .Call(`_RSTr_update_Z_mst`, Z, G, theta, beta, rho, tau2, adjacency, num_adj, island_region, island_id)
}

update_G_mst <- function(G, Z, Ag, rho, G_df, adjacency, num_island) {
    .Call(`_RSTr_update_G_mst`, G, Z, Ag, rho, G_df, adjacency, num_island)
}

update_Ag_mst <- function(Ag, G, Ag_scale, G_df, Ag_df) {
    .Call(`_RSTr_update_Ag_mst`, Ag, G, Ag_scale, G_df, Ag_df)
}

update_tau2_mst <- function(tau2, theta, beta, Z, tau_a, tau_b, island_id) {
    .Call(`_RSTr_update_tau2_mst`, tau2, theta, beta, Z, tau_a, tau_b, island_id)
}

update_theta_mst <- function(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method) {
    .Call(`_RSTr_update_theta_mst`, theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method)
}

update_rho_mst <- function(rho, r_accept, G, Z, rho_a, rho_b, rho_sd, adjacency, num_island) {
    .Call(`_RSTr_update_rho_mst`, rho, r_accept, G, Z, rho_a, rho_b, rho_sd, adjacency, num_island)
}

update_Z_u <- function(Z, sig2, theta, beta, tau2, adjacency, num_adj, island_region, island_id) {
    .Call(`_RSTr_update_Z_u`, Z, sig2, theta, beta, tau2, adjacency, num_adj, island_region, island_id)
}

update_sig2_u <- function(sig2, Z, beta, tau2, adjacency, num_adj, island_region, num_island_region, A, m0, sig_a, sig_b, method) {
    .Call(`_RSTr_update_sig2_u`, sig2, Z, beta, tau2, adjacency, num_adj, island_region, num_island_region, A, m0, sig_a, sig_b, method)
}

update_tau2_u <- function(tau2, theta, beta, Z, sig2, num_island_region, island_id, A, m0, tau_a, tau_b, method) {
    .Call(`_RSTr_update_tau2_u`, tau2, theta, beta, Z, sig2, num_island_region, island_id, A, m0, tau_a, tau_b, method)
}

update_theta_u <- function(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method) {
    .Call(`_RSTr_update_theta_u`, theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method)
}

update_beta_u <- function(beta, theta, Z, tau2, sig2, A, m0, island_region, method) {
    .Call(`_RSTr_update_beta_u`, beta, theta, Z, tau2, sig2, A, m0, island_region, method)
}

