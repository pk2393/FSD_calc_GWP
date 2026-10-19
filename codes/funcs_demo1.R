lg_cache <- function(rad, N, R, p_list, k) {
  theta <- 2 * pi * (0:(N - 1)) / N
  H1 <- rad[1] * exp(1i * theta)
  H2 <- rad[2] * exp(1i * theta)
  
  H1_grid <- matrix(H1, nrow = N, ncol = N, byrow = FALSE)
  H2_grid <- matrix(H2, nrow = N, ncol = N, byrow = TRUE)
  
  denom_1 <- 1 + (R[1] / k) * (1 - p_list[[1]][1] * H1_grid - p_list[[1]][2] * H2_grid)
  denom_2 <- 1 + (R[2] / k) * (1 - p_list[[2]][1] * H1_grid - p_list[[2]][2] * H2_grid)
  
  a11 <- 1 - (R[1] * p_list[[1]][1] * H1_grid) / denom_1
  a12 <-   - (R[1] * p_list[[1]][2] * H1_grid) / denom_1
  a21 <-   - (R[2] * p_list[[2]][1] * H2_grid) / denom_2
  a22 <- 1 - (R[2] * p_list[[2]][2] * H2_grid) / denom_2
  
  list(
    N = N,
    rad = rad,
    H1_grid = H1_grid,
    H2_grid = H2_grid,
    core_detJ = a11 * a22 - a21 * a12,
    G1 = denom_1^(-k),
    G2 = denom_2^(-k)
  )
}


precompute_pows <- function(base_mat, dmax) {
  N <- nrow(base_mat)
  out <- vector("list", dmax + 1L)
  out[[1]] <- matrix(1+0i, nrow = N, ncol = N)

  if (dmax >= 1L) {
    for (e in 1:dmax) {
      out[[e + 1L]] <- out[[e]] * base_mat
    }
  }

  out
}

lg_final_size_table <- function(dmax, index, cache, clamp0 = TRUE) {
  N <- cache$N
  G1_pows <- precompute_pows(cache$G1, dmax)
  G2_pows <- precompute_pows(cache$G2, dmax)
  H1_pows <- precompute_pows(cache$H1_grid, dmax)
  H2_pows <- precompute_pows(cache$H2_grid, dmax)
  
  pmf <- matrix(0, nrow = dmax + 1L, ncol = dmax + 1L)
  rownames(pmf) <- 0:dmax
  colnames(pmf) <- 0:dmax
  
  for (n1 in 0:dmax) {
    for (n2 in 0:dmax) {
      if (index == 1L && n1 < 1L) next
      if (index == 2L && n2 < 1L) next
      
      pow_H1 <- if (index == 1L) n1 - 1L else n1
      pow_H2 <- if (index == 1L) n2 else n2 - 1L
      
      num <- cache$core_detJ * G1_pows[[n1 + 1L]] * G2_pows[[n2 + 1L]]
      den <- (N^2) * H1_pows[[pow_H1 + 1L]] * H2_pows[[pow_H2 + 1L]]
      p <- Re(sum(num / den))
      
      if (clamp0 && p < 0) p <- 0
      pmf[n1 + 1L, n2 + 1L] <- p
    }
  }
  
  pmf
}

run_scenario_LG <- function(rho, r1, p11, p22, k = 1,
                            dmaxs = c(5, 10, 20),
                            eps_rad = 1e-2,
                            rad_cap = 0.95) {
  
  ngm <- make_ngm(rho = rho, r1 = r1, p11 = p11, p22 = p22)
  pars <- offspring_params(ngm)
  R <- pars$R
  p_list <- pars$p_list
  Q <- extinction_prob(R = R, p_list = p_list, k = k)
  
  rad <- pmin(Q * (1 - eps_rad), rad_cap)
  N <- pick_fft_N(max(dmaxs))
  cache <- lg_cache(rad = rad, N = N, R = R, p_list = p_list, k = k)
  
  out_tables <- list()
  for (dmax in dmaxs) {
    tab1 <- lg_final_size_table(dmax = dmax, index = 1, cache = cache)
    tab2 <- lg_final_size_table(dmax = dmax, index = 2, cache = cache)
    
    out_tables[[paste0("dmax", dmax)]] <- list(
      N = N,
      dmax = dmax,
      tab_index1 = tab1,
      tab_index2 = tab2,
      sum_index1 = sum(tab1),
      sum_index2 = sum(tab2)
    )
  }
  
  list(
    params = list(rho = rho, r1 = r1, p11 = p11, p22 = p22, k = k),
    ngm = ngm,
    R = R,
    p_list = p_list,
    Q = Q,
    rad = rad,
    tables = out_tables
  )
}

run_stepA_scenario_solveH <- function(rho, r1, p11, p22, k = 1,
                                      dmaxs = c(5, 10, 20),
                                      rad_cap = 0.95) {
  ngm <- make_ngm(rho = rho, r1 = r1, p11 = p11, p22 = p22)
  pars <- offspring_params(ngm)
  R <- pars$R
  p_list <- pars$p_list
  Q <- extinction_prob(R = R, p_list = p_list, k = k)
  
  rad <- choose_rad_solveH(
    Q = Q,
    R = R,
    p_list = p_list,
    k = k,
    rad_cap = rad_cap,
    target_rho = 0.975,
    par0 = 0.2,
    maxeval = 10000
  )
  
  N <- pick_fft_N(max(dmaxs))
  tables_fft <- precompute_tables(N, R, p_list, k, rad)
  
  out_tables <- list()
  for (dmax in dmaxs) {
    tab1 <- final_size_table_uncond(dmax = dmax, index = 1, tables = tables_fft)
    tab2 <- final_size_table_uncond(dmax = dmax, index = 2, tables = tables_fft)
    
    out_tables[[paste0("dmax", dmax)]] <- list(
      N = N,
      dmax = dmax,
      tab_index1 = tab1,
      tab_index2 = tab2,
      sum_index1 = sum(tab1),
      sum_index2 = sum(tab2)
    )
  }
  
  list(
    params = list(rho = rho, r1 = r1, p11 = p11, p22 = p22, k = k),
    ngm = ngm,
    R = R,
    p_list = p_list,
    Q = Q,
    rad = rad,
    tables = out_tables,
    tables_fft = tables_fft
  )
}
