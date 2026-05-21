## Utility functions for two-type Galton-Watson calculations.

NM_PGF <- function(z, R, k, p) {
  (1 + (R / k) * (1 - sum(p * z)))^(-k)
}

spectral_radius <- function(M) {
  max(Re(eigen(M, only.values = TRUE)$values))
}

normalize_by_spectral_radius <- function(M) {
  M / spectral_radius(M)
}

make_ngm <- function(rho, r1, p11, p22) {
  ## Columns are parent types; rows are offspring types.
  B <- matrix(
    c(r1 * p11,       (1 - r1) * (1 - p22),
      r1 * (1 - p11), (1 - r1) * p22),
    nrow = 2,
    byrow = TRUE
  )
  
  rho * B / spectral_radius(B)
}

offspring_params <- function(ngm) {
  R <- colSums(ngm)
  p_list <- list(ngm[, 1] / R[1], ngm[, 2] / R[2])
  list(R = R, p_list = p_list)
}

extinction_prob <- function(R, p_list, k, tol = 1e-10) {
  eq_fixed_point <- function(par) {
    G1 <- NM_PGF(par, R = R[1], k = k, p = p_list[[1]])
    G2 <- NM_PGF(par, R = R[2], k = k, p = p_list[[2]])
    c(par[1] - G1, par[2] - G2)
  }
  
  fit <- nleqslv::nleqslv(
    x = c(1e-2, 1e-2),
    fn = eq_fixed_point,
    control = list(ftol = tol, xtol = tol)
  )
  
  pmin(pmax(fit$x, 0), 1)
}

pick_fft_N <- function(dmax) {
  as.integer(max(32, 2^(ceiling(log2(max(1, dmax))) + 1)))
}
