#' Check initial values MSTCAR
#'
#' @noRd
#'
check_inits_u <- function(inits, data, num_island) {
  message("Checking inits...")
  Y <- data$Y
  num_region <- dim(Y)[[1]]
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  theta <- inits$theta
  beta <- inits$beta
  sig2 <- inits$sig2
  tau2 <- inits$tau2
  Z <- inits$Z
  chk <- c("theta", "beta", "tau2", "sig2", "Z")
  miss <- sapply(1:length(chk), \(x) !any(names(inits) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'inits': ", paste(chk[miss], collapse = ", "))
  }
  # Check for warnings
  warnout <- character()
  # Check for unused elements in 'inits'
  chk_elem <- which(!(names(inits) %in% chk))
  if (length(chk_elem)) {
    warnout <- c(warnout, paste("Unused elements of list 'inits':", paste(names(inits)[chk_elem], collapse = ", ")))
  }
  if (length(warnout) != 0) {
    warnout <- paste(seq_along(warnout), ":", warnout)
    warning(paste(length(warnout), "warning(s) found in list 'inits':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout <- character()
  # theta
  # dimensions don't match num_region num_group num_time
  if (!all(dim(theta) == dim(Y))) {
    errout <- c(errout, "theta is not a num_region x num_group x num_time array. Ensure dim(theta) == dim(Y) or use default value")
  }
  # values are infinite
  if (any(!is.finite(theta))) {
    errout <- c(errout, "theta contains infinite values. Ensure all(is.finite(theta)) or use default value")
  }

  # beta
  # length doesn't match num_island
  if (!all(dim(beta) == c(num_island, num_group, num_time))) {
    errout <- c(errout, "beta is not an num_island x num_group x num_time array. Ensure dim(beta) == num_island x num_group x num_time or use default value")
  }
  # values are infinite
  if (any(!is.finite(beta))) {
    errout <- c(errout, "beta contains infinite values. Ensure all(is.finite(beta)) or use default value")
  }
  # sig2
  # is non-positive or infinite
  if (any(sig2 <= 0) | any(!is.finite(sig2))) {
    errout <- c(errout, "Some or all sig2 are non-positive or infinite. Ensure all sig2 > 0 and not infinite or use default value")
  }

  # tau2
  # is non-positive or infinite
  if (any(tau2 <= 0) | any(!is.finite(tau2))) {
    errout <- c(errout, "Some or all tau2 are non-positive or infinite. Ensure all tau2 > 0 and not infinite or use default value")
  }

  # Z
  # length doesn't match num_region
  if (!all(dim(Z) == dim(Y))) {
    errout <- c(errout, "Z is not an num_region x num_group x num_time array. Ensure dim(Z) == dim(Y) or use default value")
  }
  # values are infinite
  if (any(!is.finite(Z))) {
    errout <- c(errout, "Z contains infinite values. Ensure all(is.finite(Z)) or use default value")
  }

  if (length(errout) != 0) {
    errout <- paste(seq_along(errout), ":", errout)
    stop(paste(length(errout), "error(s) found in list 'inits':\n", paste(errout, collapse = "\n ")))
  }
}

#' Check initial values MCAR
#'
#' @noRd
#'
check_inits_m <- function(inits, num_region, num_group, num_island) {
  message("Checking inits...")
  theta <- inits$theta
  beta <- inits$beta
  G <- inits$G
  tau2 <- inits$tau2
  Z <- inits$Z
  chk <- c("theta", "beta", "tau2", "G", "Z")
  miss <- sapply(1:length(chk), \(x) !any(names(inits) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'inits': ", paste(chk[miss], collapse = ", "))
  }
  # Check for warnings
  warnout <- character()
  # Check for unused elements in 'inits'
  chk_elem <- which(!(names(inits) %in% chk))
  if (length(chk_elem)) {
    warnout <- c(warnout, paste("Unused elements of list 'inits':", paste(names(inits)[chk_elem], collapse = ", ")))
  }
  if (length(warnout) != 0) {
    warnout <- paste(seq_along(warnout), ":", warnout)
    warning(paste(length(warnout), "warning(s) found in list 'inits':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout <- character()
  # theta
  # dimensions don't match num_region num_group
  if (!all(dim(theta) == c(num_region, num_group))) {
    errout <- c(errout, "theta is not an num_region x num_group array. Ensure dim(theta) == dim(Y) or use default value")
  }
  # values are infinite
  if (any(!is.finite(theta))) {
    errout <- c(errout, "theta contains infinite values. Ensure all(is.finite(theta)) or use default value")
  }

  # beta
  # dimensions don't match num_island num_group
  if (!all(dim(beta) == c(num_island, num_group))) {
    errout <- c(errout, "beta is not an num_island x num_group array. Ensure dim(beta) == num_island x num_group or use default value")
  }
  # values are infinite
  if (any(!is.finite(beta))) {
    errout <- c(errout, "beta contains infinite values. Ensure all(is.finite(beta)) or use default value")
  }
  # G
  sig2 <- diag(G)
  gcor <- G[lower.tri(G)]
  # G not symmetric
  if (!isSymmetric(G)) {
    errout <- c(errout, "G is not symmetric. Ensure G is symmetric or use default value")
  }

  # sig2
  # is non-positive or infinite
  if (any((sig2 <= 0) | !is.finite(sig2))) {
    errout <- c(errout, "Diagonals of G contain non-positive values. Ensure all diag(G) > 0 and not infinite or use default value")
  }

  # gcor
  # values are infinite
  if (any(!is.finite(gcor))) {
    errout <- c(errout, "Off-diagonals of G contain infinite values. Ensure all(is.finite(G)) or use default value")
  }

  # tau2
  # length not num_group
  if (length(tau2) != num_group) {
    errout <- c(errout, "tau2 is not length num_group. Ensure length(tau2) == num_group or use default value")
  }
  # is non-positive or infinite
  if (any((tau2 <= 0) | !is.finite(tau2))) {
    errout <- c(errout, "tau2 contains non-positive values. Ensure all(tau2 > 0) and not infinite or use default value")
  }

  # Z
  # dimensions don't match num_region num_group
  if (!all(dim(Z) == c(num_region, num_group))) {
    errout <- c(errout, "Z is not an num_region x num_group array. Ensure dim(Z) == dim(Y) or use default value")
  }
  # values are infinite
  if (any(!is.finite(Z))) {
    errout <- c(errout, "Z contains infinite values. Ensure all(is.finite(Z)) or use default value")
  }

  if (length(errout) != 0) {
    errout <- paste(seq_along(errout), ":", errout)
    stop(paste(length(errout), "error(s) found in list 'inits':\n", paste(errout, collapse = "\n ")))
  }
}


#' Check initial values MSTCAR
#'
#' @noRd
#'
check_inits_mst <- function(inits, num_region, num_group, num_time, num_island) {
  message("Checking inits...")
  theta <- inits$theta
  beta <- inits$beta
  G <- inits$G
  sig2 <- apply(G, 3, diag)
  gcor <- apply(G, 3, \(G) G[lower.tri(G)])
  tau2 <- inits$tau2
  Z <- inits$Z
  rho <- inits$rho
  Ag <- inits$Ag
  chk <- c("theta", "beta", "tau2", "G", "Ag", "Z", "rho")
  miss <- sapply(1:length(chk), \(x) !any(names(inits) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'inits': ", paste(chk[miss], collapse = ", "))
  }
  # Check for warnings
  warnout <- character()
  # Check for unused elements in 'inits'
  chk_elem <- which(!(names(inits) %in% chk))
  if (length(chk_elem)) {
    warnout <- c(warnout, paste("Unused elements of list 'inits':", paste(names(inits)[chk_elem], collapse = ", ")))
  }
  if (length(warnout) != 0) {
    warnout <- paste(seq_along(warnout), ":", warnout)
    warning(paste(length(warnout), "warning(s) found in list 'inits':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout <- character()
  # theta
  # dimensions don't match num_region num_group num_time
  if (!all(dim(theta) == c(num_region, num_group, num_time))) {
    errout <- c(errout, "theta is not an num_region x num_group x num_time array. Ensure dim(theta) == dim(Y) or use default value")
  }
  # values are infinite
  if (any(!is.finite(theta))) {
    errout <- c(errout, "theta contains infinite values. Ensure all(is.finite(theta)) or use default value")
  }

  # beta
  # dimensions don't match num_island num_group num_time
  if (!all(dim(beta) == c(num_island, num_group, num_time))) {
    errout <- c(errout, "beta is not an num_island x num_group x num_time array. Ensure dim(beta) == num_island x num_group x num_time or use default value")
  }
  # values are infinite
  if (any(!is.finite(beta))) {
    errout <- c(errout, "beta contains infinite values. Ensure all(is.finite(beta)) or use default value")
  }
  # sig2
  # is non-positive or infinite
  if (any((sig2 <= 0) | !is.finite(sig2))) {
    errout <- c(errout, "Diagonals of G contain non-positive values. Ensure all diag(G) > 0 and not infinite or use default value")
  }

  # gcor
  # values are infinite
  if (any(!is.finite(gcor))) {
    errout <- c(errout, "Off-diagonals of G contain infinite values. Ensure all(is.finite(G)) or use default value")
  }

  # tau2
  # length not num_group
  if (length(tau2) != num_group) {
    errout <- c(errout, "tau2 is not length num_group. Ensure length(tau2) == num_group or use default value")
  }
  # is non-positive or infinite
  if (any((tau2 <= 0) | !is.finite(tau2))) {
    errout <- c(errout, "tau2 contains non-positive values. Ensure all(tau2 > 0) and not infinite or use default value")
  }

  # rho
  # length not num_group
  if (length(rho) != num_group) {
    errout <- c(errout, "rho is not length num_group. Ensure length(rho) == num_group or use default value")
  }
  # is non-positive or infinite
  if (any((rho <= 0) | !is.finite(rho))) {
    errout <- c(errout, "rho contains non-positive values. Ensure all(rho > 0) and not infinite or use default value")
  }

  # Z
  # dimensions don't match num_region num_group num_time
  if (!all(dim(Z) == c(num_region, num_group, num_time))) {
    errout <- c(errout, "Z is not an num_region x num_group x num_time array. Ensure dim(Z) == dim(Y) or use default value")
  }
  # values are infinite
  if (any(!is.finite(Z))) {
    errout <- c(errout, "Z contains infinite values. Ensure all(is.finite(Z)) or use default value")
  }

  # Ag
  # dimensions don't match num_group num_group
  if (!all(dim(Ag) == c(num_group, num_group))) {
    errout <- c(errout, "Ag is not an num_group x num_group matrix. Ensure dim(Ag) == num_group x num_group or use default value")
  }
  # matrix is not symmetric
  if (!isSymmetric(Ag)) {
    errout <- c(errout, "Ag is not symmetric. Ensure Ag is symmetric or use default value")
  }
  # values are infinite
  if (any(!is.finite(Ag))) {
    errout <- c(errout, "Ag contains infinite values. Ensure Ag is finite or use default value")
  }
  # diagonals are not positive
  if (any(diag(Ag) <= 0)) {
    errout <- c(errout, "diag(Ag) contains non-positive values. Ensure diag(Ag) is positive or use default value")
  }

  if (length(errout) != 0) {
    errout <- paste(seq_along(errout), ":", errout)
    stop(paste(length(errout), "error(s) found in list 'inits':\n", paste(errout, collapse = "\n ")))
  }
}
