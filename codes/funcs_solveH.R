fft2 <- function(M) {
  tmp <- mvfft(M)
  t(mvfft(t(tmp)))
}

solve_H_point_newton <- function(z1, z2, R, p_list, k,
                                 h_init = c(0+0i, 0+0i),
                                 tol = 1e-20, maxit = 50) {
  h1 <- h_init[1]
  h2 <- h_init[2]
  R_C <- R[1]
  R_A <- R[2]
  p_C <- p_list[[1]]
  p_A <- p_list[[2]]
  
  h_bound <- c(Mod(z1), Mod(z2))
  inside_region <- function(h) {
    all(Mod(h) <= h_bound * (1 + 1e-12) + 1e-14)
  }
  
  for (it in 1:maxit) {
    Gc <- NM_PGF(c(h1, h2), R_C, k, p_C)
    Ga <- NM_PGF(c(h1, h2), R_A, k, p_A)
    F1 <- h1 - z1 * Gc
    F2 <- h2 - z2 * Ga
    current_ss <- sum(Mod(c(F1, F2))^2)
    
    if (current_ss < tol) break
    
    denomC <- 1 + (R_C / k) * (1 - (p_C[1] * h1 + p_C[2] * h2))
    denomA <- 1 + (R_A / k) * (1 - (p_A[1] * h1 + p_A[2] * h2))
    
    dGc1 <- R_C * p_C[1] * denomC^(-k - 1)
    dGc2 <- R_C * p_C[2] * denomC^(-k - 1)
    dGa1 <- R_A * p_A[1] * denomA^(-k - 1)
    dGa2 <- R_A * p_A[2] * denomA^(-k - 1)
    
    ## Jacobian of F(h) = h - diag(z) G(h).
    J <- rbind(
      c(1 - z1 * dGc1, -z1 * dGc2),
      c(-z2 * dGa1,    1 - z2 * dGa2)
    )
    
    delta <- solve(J, -c(F1, F2))
    step <- 1
    accepted <- FALSE
    
    for (ls in 1:5) {
      h_try <- c(h1, h2) + step * delta
      
      if (inside_region(h_try)) {
        Gc_try <- NM_PGF(h_try, R_C, k, p_C)
        Ga_try <- NM_PGF(h_try, R_A, k, p_A)
        F_try <- c(h_try[1] - z1 * Gc_try,
                   h_try[2] - z2 * Ga_try)
        
        if (sum(Mod(F_try)^2) < current_ss) {
          h1 <- h_try[1]
          h2 <- h_try[2]
          accepted <- TRUE
          break
        }
      }
      
      step <- step / 2
    }
    
    if (!accepted) {
      h1 <- z1 * Gc
      h2 <- z2 * Ga
    }
  }
  
  c(h1, h2)
}

choose_rad_solveH <- function(Q, R, p_list, k,
                              rad_cap = 0.95,
                              target_rho = 0.975,
                              par0 = 0.2,
                              maxeval = 10000) {
  if (all(Q >= 1 - 1e-15)) {
    return(rep(rad_cap, 2))
  }
  
  rho_JG <- function(h) {
    denom1 <- 1 + (R[1] / k) * (1 - (p_list[[1]][1] * h[1] + p_list[[1]][2] * h[2]))
    denom2 <- 1 + (R[2] / k) * (1 - (p_list[[2]][1] * h[1] + p_list[[2]][2] * h[2]))
    
    J <- matrix(
      c(R[1] * p_list[[1]][1] * denom1^(-k - 1),
        R[1] * p_list[[1]][2] * denom1^(-k - 1),
        R[2] * p_list[[2]][1] * denom2^(-k - 1),
        R[2] * p_list[[2]][2] * denom2^(-k - 1)),
      nrow = 2,
      byrow = TRUE
    )
    
    max(Re(eigen(J, only.values = TRUE)$values))
  }
  
  find_rad_upper <- function(t) {
    rad_t <- Q + t * (1 - Q)
    (rho_JG(rad_t) - target_rho)^2
  }
  
  fit <- nloptr::nloptr(
    x0 = par0,
    eval_f = find_rad_upper,
    lb = 0,
    ub = 1,
    opts = list(
      algorithm = "NLOPT_LN_SBPLX",
      xtol_rel = 1e-5,
      xtol_abs = 1e-5,
      maxeval = maxeval,
      print_level = 0
    )
  )
  
  rad <- Q + fit$solution * (1 - Q)
  pmin(rad, rep(rad_cap, 2))
}

solve_H_grid <- function(N, R, p_list, k, rad,
                         tol = 1e-20, maxit = 100) {
  th <- 2 * pi * (0:(N - 1)) / N
  z1 <- rad[1] * exp(1i * th)
  z2 <- rad[2] * exp(1i * th)
  
  H1 <- matrix(0+0i, N, N)
  H2 <- matrix(0+0i, N, N)
  
  for (i in 1:N) {
    for (j in 1:N) {
      hij <- solve_H_point_newton(
        z1 = z1[i],
        z2 = z2[j],
        R = R,
        p_list = p_list,
        k = k,
        h_init = c(0+0i, 0+0i),
        tol = tol,
        maxit = maxit
      )
      
      H1[i, j] <- hij[1]
      H2[i, j] <- hij[2]
    }
  }
  
  list(H1 = H1, H2 = H2, rad = rad, N = N)
}

precompute_tables <- function(N, R, p_list, k, rad) {
  HG <- solve_H_grid(N, R, p_list, k, rad)
  
  list(
    F1 = fft2(HG$H1) / (N^2),
    F2 = fft2(HG$H2) / (N^2),
    N = N,
    rad = rad
  )
}

coef_nm <- function(index, nC, nA, tables, clamp0 = TRUE) {
  r <- tables$rad
  F <- if (index == 1) tables$F1 else tables$F2
  
  p <- Re(F[nC + 1, nA + 1]) * r[1]^(-nC) * r[2]^(-nA)
  if (clamp0 && is.finite(p) && p < 0) p <- 0
  p
}

final_size_table_uncond <- function(dmax, index, tables, clamp0 = TRUE) {
  pmf <- matrix(0, nrow = dmax + 1L, ncol = dmax + 1L)
  rownames(pmf) <- 0:dmax
  colnames(pmf) <- 0:dmax
  
  for (nC in 0:dmax) {
    for (nA in 0:dmax) {
      if (index == 1L && nC < 1L) next
      if (index == 2L && nA < 1L) next
      
      pmf[nC + 1L, nA + 1L] <- coef_nm(index, nC, nA, tables, clamp0 = clamp0)
    }
  }
  
  pmf
}
