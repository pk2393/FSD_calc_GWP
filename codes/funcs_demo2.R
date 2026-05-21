## Demo 2 helper functions: simulation, likelihood, and plotting.

conditional_offspring_params <- function(R, p_list, k, Q) {
  mu_star <- numeric(2)
  p_star <- vector("list", 2)
  s0 <- numeric(2)
  
  for (i in 1:2) {
    s0[i] <- sum(p_list[[i]] * Q)
    mu_star[i] <- R[i] * s0[i] / (1 + (R[i] / k) * (1 - s0[i]))
    p_star[[i]] <- (p_list[[i]] * Q) / s0[i]
  }
  
  list(mu = mu_star, p_list = p_star, s0 = s0)
}

draw_offspring_sums_2G <- function(n_parents, mu, k, p_child) {
  if (n_parents <= 0L) return(c(0L, 0L))
  
  n_total <- rnbinom(n_parents, size = k, mu = mu)
  n_child <- rbinom(n_parents, size = n_total, prob = p_child)
  
  c(
    as.integer(sum(n_child)),
    as.integer(sum(n_total - n_child))
  )
}

simulate_cluster_2G_NM <- function(index,
                                   R, p_list, k, Q,
                                   max_gen = 200L,
                                   max_total = 1000L,
                                   condition_on_extinction = TRUE) {
  gen <- 0L
  g_now <- if (index == 1L) c(1L, 0L) else c(0L, 1L)
  total <- g_now
  truncated <- FALSE
  
  ## For supercritical settings, this option simulates finite clusters
  ## from the offspring law conditioned on eventual extinction.
  if (condition_on_extinction && any(Q < 1 - 1e-12)) {
    cond <- conditional_offspring_params(R, p_list, k, Q)
    R_use <- cond$mu
    p_use <- cond$p_list
  } else {
    R_use <- R
    p_use <- p_list
  }
  
  while (sum(g_now) > 0L) {
    gen <- gen + 1L
    if (gen > max_gen) {
      truncated <- TRUE
      break
    }
    
    off1 <- draw_offspring_sums_2G(g_now[1], R_use[1], k, p_use[[1]][1])
    off2 <- draw_offspring_sums_2G(g_now[2], R_use[2], k, p_use[[2]][1])
    
    g_now <- off1 + off2
    total <- total + g_now
    
    if (sum(total) > max_total) {
      truncated <- TRUE
      break
    }
  }
  
  list(
    index = as.integer(index),
    n_child = as.integer(total[1]),
    n_adult = as.integer(total[2]),
    n_total = as.integer(sum(total)),
    extinct = (!truncated && sum(g_now) == 0L),
    truncated = truncated,
    n_gen = as.integer(gen)
  )
}

simulate_dataset_2G_NM <- function(n_index = c(20L, 20L),
                                   index_vec = NULL,
                                   R, p_list, k, Q,
                                   seed = NULL,
                                   shuffle = TRUE,
                                   max_gen = 200L,
                                   max_total = 1000L,
                                   condition_on_extinction = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(index_vec)) {
    index_vec <- c(rep(1L, n_index[1]), rep(2L, n_index[2]))
    if (shuffle) index_vec <- sample(index_vec)
  } else {
    index_vec <- as.integer(index_vec)
  }
  
  out <- vector("list", length(index_vec))
  for (i in seq_along(index_vec)) {
    out[[i]] <- simulate_cluster_2G_NM(
      index = index_vec[i],
      R = R,
      p_list = p_list,
      k = k,
      Q = Q,
      max_gen = max_gen,
      max_total = max_total,
      condition_on_extinction = condition_on_extinction
    )
  }
  
  df <- data.frame(
    cluster_id = seq_along(index_vec),
    index = vapply(out, `[[`, integer(1), "index"),
    n_child = vapply(out, `[[`, integer(1), "n_child"),
    n_adult = vapply(out, `[[`, integer(1), "n_adult"),
    n_total = vapply(out, `[[`, integer(1), "n_total"),
    extinct = vapply(out, `[[`, logical(1), "extinct"),
    truncated = vapply(out, `[[`, logical(1), "truncated"),
    n_gen = vapply(out, `[[`, integer(1), "n_gen")
  )
  
  attr(df, "Q") <- Q
  df
}

sim_df_to_complete_lists <- function(sim_df) {
  make_one <- function(sub) {
    if (nrow(sub) == 0L) {
      return(data.frame(n_child = integer(0), n_adult = integer(0), count = integer(0)))
    }
    
    tab <- as.data.frame(
      table(n_child = sub$n_child, n_adult = sub$n_adult),
      stringsAsFactors = FALSE
    )
    tab$count <- as.integer(tab$Freq)
    tab$Freq <- NULL
    tab$n_child <- as.integer(tab$n_child)
    tab$n_adult <- as.integer(tab$n_adult)
    tab[tab$count > 0L, ]
  }
  
  list(
    list_complete_C = make_one(sim_df[sim_df$index == 1L & !sim_df$truncated, ]),
    list_complete_A = make_one(sim_df[sim_df$index == 2L & !sim_df$truncated, ])
  )
}

simulate_scenario_2G_NM <- function(rho, ngm_base, k,
                                    n_index = c(20L, 20L),
                                    seed = 1L,
                                    max_gen = 200L,
                                    max_total = 1000L,
                                    condition_on_extinction = TRUE) {
  ngm <- rho * ngm_base
  pars <- offspring_params(ngm)
  R <- pars$R
  p_list <- pars$p_list
  Q <- if (rho <= 1) c(1, 1) else extinction_prob(R, p_list, k)
  
  sim_df <- simulate_dataset_2G_NM(
    n_index = n_index,
    R = R,
    p_list = p_list,
    k = k,
    Q = Q,
    seed = seed,
    shuffle = TRUE,
    max_gen = max_gen,
    max_total = max_total,
    condition_on_extinction = condition_on_extinction
  )
  
  list(
    params = list(rho = rho, k = k),
    ngm = ngm,
    R = R,
    p_list = p_list,
    Q = Q,
    sim = list(
      n_index = n_index,
      seed = seed,
      max_gen = max_gen,
      max_total = max_total,
      df = sim_df
    )
  )
}

simulate_no_trunc <- function(rho_true, k_true, seed0, n_index,
                              ngm_base,
                              max_gen, max_total, max_resim,
                              condition_on_extinction = TRUE,
                              verbose = TRUE) {
  seed_use <- as.integer(seed0)
  
  for (attempt in 1:max_resim) {
    sim_obj <- simulate_scenario_2G_NM(
      rho = rho_true,
      ngm_base = ngm_base,
      k = k_true,
      n_index = n_index,
      seed = seed_use,
      max_gen = max_gen,
      max_total = max_total,
      condition_on_extinction = condition_on_extinction
    )
    
    sim_df <- sim_obj$sim$df
    n_trunc <- sum(sim_df$truncated)
    
    if (n_trunc == 0L) {
      return(list(
        sim_obj = sim_obj,
        sim_df = sim_df,
        n_trunc = n_trunc,
        seed = seed_use,
        attempt = attempt
      ))
    }
    
    if (verbose) {
      cat(sprintf(
        "[resim] attempt %d/%d: truncated=%d -> seed=%d\n",
        attempt, max_resim, n_trunc, seed_use + 1L
      ))
    }
    seed_use <- seed_use + 1L
  }
  
  warning(sprintf(
    "Still truncated after %d attempts (max_total=%d, max_gen=%d).",
    max_resim, max_total, max_gen
  ))
  
  list(
    sim_obj = sim_obj,
    sim_df = sim_df,
    n_trunc = n_trunc,
    seed = seed_use,
    attempt = max_resim
  )
}

choose_dmax_from_lists <- function(listC, listA) {
  vals <- c(
    if (nrow(listC) > 0L) c(listC$n_child, listC$n_adult) else integer(0),
    if (nrow(listA) > 0L) c(listA$n_child, listA$n_adult) else integer(0)
  )
  
  dmax <- if (length(vals) == 0L) 1L else max(vals, na.rm = TRUE)
  dmax <- max(1L, as.integer(dmax))
  list(dmax = dmax, N = pick_fft_N(dmax))
}

make_loglik_uncond <- function(listC, listA, N, ngm_base,
                               rad_cap = 0.95,
                               target_rho = 0.975) {
  cache_env <- new.env(parent = emptyenv())
  
  function(par) {
    rho <- exp(par[1])
    k <- exp(par[2])
    key <- paste(sprintf("%.8f", par), collapse = "|")
    
    if (exists(key, envir = cache_env, inherits = FALSE)) {
      return(get(key, envir = cache_env, inherits = FALSE))
    }
    
    ngm <- rho * ngm_base
    pars <- offspring_params(ngm)
    R <- pars$R
    p_list <- pars$p_list
    Q <- if (rho <= 1) c(1, 1) else extinction_prob(R, p_list, k)
    rad <- choose_rad_solveH(Q, R, p_list, k,
                             rad_cap = rad_cap,
                             target_rho = target_rho)
    tables <- precompute_tables(N, R, p_list, k, rad)
    ll <- 0
    
    if (nrow(listC) > 0L) {
      for (ii in seq_len(nrow(listC))) {
        p <- coef_nm(1, listC$n_child[ii], listC$n_adult[ii], tables)
        if (!is.finite(p) || p <= 0) {
          ll <- -Inf
          break
        }
        ll <- ll + log(p) * listC$count[ii]
      }
    }
    
    if (is.finite(ll) && nrow(listA) > 0L) {
      for (ii in seq_len(nrow(listA))) {
        p <- coef_nm(2, listA$n_child[ii], listA$n_adult[ii], tables)
        if (!is.finite(p) || p <= 0) {
          ll <- -Inf
          break
        }
        ll <- ll + log(p) * listA$count[ii]
      }
    }
    
    assign(key, ll, envir = cache_env)
    ll
  }
}

loglik_safe <- function(loglik_fun, par) {
  ll <- tryCatch(loglik_fun(par), error = function(e) NA_real_)
  if (!is.finite(ll)) return(-Inf)
  ll
}

fit_mle_nloptr <- function(loglik_fun, starts, lb, ub, opts) {
  obj <- function(x) {
    ll <- loglik_safe(loglik_fun, x)
    if (!is.finite(ll)) return(1e100)
    -ll
  }
  
  alg <- if (!is.null(opts$algorithm)) as.character(opts$algorithm) else ""
  if (grepl("^NLOPT_GN_DIRECT", alg)) {
    starts <- starts[1, , drop = FALSE]
  }
  
  best_obj <- Inf
  best_sol <- NULL
  best_status <- NA_integer_
  last_err <- NULL
  
  for (i in seq_len(nrow(starts))) {
    x0 <- pmin(pmax(starts[i, ], lb), ub)
    
    res <- tryCatch(
      nloptr::nloptr(
        x0 = x0,
        eval_f = obj,
        lb = lb,
        ub = ub,
        opts = opts
      ),
      error = function(e) {
        last_err <<- e$message
        NULL
      }
    )
    
    if (is.null(res)) next
    
    if (is.finite(res$objective) && res$objective < best_obj) {
      best_obj <- res$objective
      best_sol <- res$solution
      best_status <- res$status
    }
  }
  
  if (!is.finite(best_obj) || is.null(best_sol)) {
    return(list(ok = FALSE, err = last_err))
  }
  
  list(
    ok = TRUE,
    par = best_sol,
    logrho = best_sol[1],
    logk = best_sol[2],
    rho = exp(best_sol[1]),
    k = exp(best_sol[2]),
    logLik = -best_obj,
    status = best_status,
    opts = opts
  )
}

eval_ll_grid <- function(logrho_grid, logk_grid, loglik_fun) {
  df <- expand.grid(log_rho = logrho_grid, log_k = logk_grid)
  df$ll <- NA_real_
  df$ok <- FALSE
  
  for (i in seq_len(nrow(df))) {
    ll <- tryCatch(
      loglik_fun(c(df$log_rho[i], df$log_k[i])),
      error = function(e) NA_real_
    )
    df$ll[i] <- ll
    df$ok[i] <- is.finite(ll)
  }
  
  df
}

hpd_thresh_unnorm <- function(pmf, alpha = 0.95) {
  v <- as.vector(pmf)
  v <- v[is.finite(v) & v > 0]
  if (length(v) == 0L) return(NA_real_)
  
  o <- order(v, decreasing = TRUE)
  cs <- cumsum(v[o])
  target <- alpha * sum(v[o])
  v[o][which(cs >= target)[1]]
}

plot_ll_surface_2d <- function(df_grid, mle, lrt_levels, fill_cols,
                               clip_pad = 2) {
  ll_ref <- max(mle$logLik, max(df_grid$ll[df_grid$ok], na.rm = TRUE))
  ll_floor <- min(df_grid$ll[df_grid$ok], na.rm = TRUE) - clip_pad
  
  df_grid$ll_plot <- ifelse(df_grid$ok, df_grid$ll, ll_floor)
  df_grid$ll_contour <- ifelse(df_grid$ok, df_grid$ll, NA_real_)
  df_grid$ll_rel_plot <- df_grid$ll_plot - ll_ref
  lrt_breaks <- ll_ref - stats::qchisq(lrt_levels, df = 2) / 2
  
  ggplot(df_grid, aes(x = log_rho, y = log_k)) +
    geom_raster(aes(fill = ll_rel_plot)) +
    geom_contour(aes(z = ll_contour), breaks = lrt_breaks,
                 colour = "black", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = 3, linewidth = 0.35) +
    geom_point(aes(x = mle$logrho, y = mle$logk),
               inherit.aes = FALSE,
               shape = 4, size = 3.0, stroke = 1.1,
               colour = "black") +
    coord_cartesian(expand = FALSE) +
    scale_fill_gradientn(
      colours = fill_cols,
      name = expression(llk(R, k) - llk(hat(R), hat(k))),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = grid::unit(34, "mm"),
        barheight = grid::unit(3.2, "mm")
      )
    ) +
    labs(
      title = "Likelihood Surface at MLE",
      x = expression(log(R)),
      y = expression(log(k))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      aspect.ratio = 1,
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

plot_pmf_overlay <- function(pmf, pts, index_label, mass_levels, fill_cols) {
  dmax <- nrow(pmf) - 1L
  dfp <- expand.grid(n_child = 0:dmax, n_adult = 0:dmax)
  dfp$p <- as.vector(pmf)
  dfp$p[dfp$p < 1e-8] <- NA_real_
  
  contour_breaks <- sapply(mass_levels, function(a) hpd_thresh_unnorm(pmf, a))
  contour_breaks <- contour_breaks[is.finite(contour_breaks) & contour_breaks > 0]
  int_breaks <- function(x) unique(as.integer(pretty(x)))
  
  g <- ggplot(dfp, aes(x = n_child, y = n_adult)) +
    geom_raster(aes(fill = p)) +
    coord_equal(expand = FALSE) +
    scale_x_continuous(breaks = int_breaks) +
    scale_y_continuous(breaks = int_breaks) +
    scale_fill_gradientn(
      colours = fill_cols,
      trans = "log10",
      breaks = c(1e-5, 1e-3, 1e-1),
      labels = c("1e-5", "1e-3", "1e-1"),
      na.value = fill_cols[1],
      name = "Probability",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = grid::unit(34, "mm"),
        barheight = grid::unit(3.2, "mm")
      )
    ) +
    labs(
      title = sprintf("PMF at MLE: index = %s", index_label),
      x = expression(n[C]),
      y = expression(n[A])
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
  
  if (length(contour_breaks) > 0L) {
    g <- g + geom_contour(aes(z = p), breaks = contour_breaks,
                          colour = "black", linewidth = 0.35)
  }
  
  if (!is.null(pts) && nrow(pts) > 0L) {
    pts <- pts
    pts$pt_size <- sqrt(pts$count)
    g <- g +
      geom_point(
        data = pts,
        aes(x = n_child, y = n_adult, size = pt_size),
        colour = "black",
        shape = 4,
        stroke = 1.0,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      scale_size_continuous(range = c(1.5, 4.5))
  }
  
  g
}

make_pmf_tables_at_fit <- function(fit, ngm_base, N, dmax,
                                   rad_cap = 0.95,
                                   target_rho = 0.975) {
  ngm <- fit$rho * ngm_base
  pars <- offspring_params(ngm)
  R <- pars$R
  p_list <- pars$p_list
  Q <- if (fit$rho <= 1) c(1, 1) else extinction_prob(R, p_list, fit$k)
  rad <- choose_rad_solveH(Q, R, p_list, fit$k,
                           rad_cap = rad_cap,
                           target_rho = target_rho)
  tables <- precompute_tables(N, R, p_list, fit$k, rad)
  
  list(
    ngm = ngm,
    R = R,
    p_list = p_list,
    Q = Q,
    rad = rad,
    tables = tables,
    pmf_child = final_size_table_uncond(dmax = dmax, index = 1, tables = tables),
    pmf_adult = final_size_table_uncond(dmax = dmax, index = 2, tables = tables)
  )
}

save_three_panel_pdf <- function(p1, p2, p3, file,
                                 widths = c(1, 1, 1),
                                 width = 15.0,
                                 height = 5.0) {
  grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 1,
    ncol = 3,
    widths = grid::unit(widths, "null")
  )))
  
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p3, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
}
